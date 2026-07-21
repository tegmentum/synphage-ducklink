//! DuckLink dispatch surface for the BLAST component.
//!
//! Registers two table functions via `runtime::TableRegistry` and dispatches
//! them row-major through `callback_dispatch::call_table`:
//!
//! - `blastn(queries VARCHAR, subjects VARCHAR, options VARCHAR)` — JSON args.
//! - `blastp(...)` — same signature, BLASTP-default scoring.
//!
//! Both delegate to `crate::run_search` with a fixed scoring preset
//! (`BlastnDefault` / `BlastpDefault`). Tunable scoring is future work — the
//! SQL surface stays clean for the common case first.
//!
//! ## Why row-major dispatch (call_table) instead of streaming
//!
//! ducklink-extension v5.0.0's native DuckDB path only expands table functions
//! registered on the `runtime::TableRegistry` capability — the alternative
//! `table_stream::register_filterable_table` route (which we previously used
//! for optional filter pushdown and streaming) is not surfaced. Both paths are
//! defined in the WIT, but only the runtime path is dispatched into DuckDB by
//! the current native extension.
//!
//! We eagerly materialise all hits in memory anyway (Synphage-scale inputs
//! are dozens of genomes × thousands of features), so giving up the streaming
//! cursor costs nothing. The `pushdown` module and its tests stay in the tree
//! as helpers ready to hook back in when the streaming path returns.
//!
//! ## Sample surface
//!
//! ```sql
//! LOAD 'ducklink.duckdb_extension';
//! FROM ducklink_load('blast.wasm');
//! SELECT * FROM blastn(
//!     '[{"key":"q1","data":"ACGT..."}]',
//!     '[{"key":"s1","data":"ACGT..."}]',
//!     '{"evalue_max": 1e-5}'
//! );
//! ```

use std::collections::HashMap;
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::{Mutex, OnceLock};

use serde::Deserialize;

use crate::bindings::duckdb::extension::runtime::{
    self, Capability, Extopts, TableCallback,
};
use crate::bindings::duckdb::extension::types::{
    Capabilitykind, Columndef, Complexvalue, Duckerror, Duckvalue, Funcarg, Loadresult,
    Logicaltype, Resultset,
};
use crate::bindings::exports::duckdb::extension::{
    callback_dispatch, guest, table_stream_dispatch,
};
use crate::{Component, Hit, Scoring, SearchOptions, Sequence, Strand};

/// Discriminator stored in `TABLE_HANDLERS` keyed by the callback handle we
/// mint in `load()` and pass to `TableCallback::new`. The runtime threads the
/// same handle back into `callback_dispatch::call_table`, and we use it to
/// pick which scoring preset to run.
#[derive(Copy, Clone)]
enum TableHandler {
    Blastn,
    Blastp,
}

fn table_handlers() -> &'static Mutex<HashMap<u32, TableHandler>> {
    static M: OnceLock<Mutex<HashMap<u32, TableHandler>>> = OnceLock::new();
    M.get_or_init(|| Mutex::new(HashMap::new()))
}

fn next_handle() -> u32 {
    static N: AtomicU32 = AtomicU32::new(1);
    N.fetch_add(1, Ordering::Relaxed)
}

// ---- guest lifecycle ---------------------------------------------------

impl guest::Guest for Component {
    fn load() -> Result<Loadresult, Duckerror> {
        let capability = runtime::get_capability(Capabilitykind::Table).ok_or_else(|| {
            Duckerror::Internal("blast: host did not expose table capability".into())
        })?;
        let registry = match capability {
            Capability::Table(r) => r,
            _ => {
                return Err(Duckerror::Internal(
                    "blast: table capability returned unexpected variant".into(),
                ))
            }
        };

        let cols = hit_columns();
        let args = seq_search_arg_types();

        let blastn_handle = next_handle();
        table_handlers()
            .lock()
            .unwrap()
            .insert(blastn_handle, TableHandler::Blastn);
        registry.register(
            "blastn",
            &args,
            &cols,
            TableCallback::new(blastn_handle),
            Some(&Extopts {
                description: Some(
                    "Nucleotide BLAST (rust-bio Smith-Waterman + Karlin-Altschul stats)".into(),
                ),
                tags: vec!["bioinformatics".into(), "alignment".into()],
            }),
        )?;

        let blastp_handle = next_handle();
        table_handlers()
            .lock()
            .unwrap()
            .insert(blastp_handle, TableHandler::Blastp);
        registry.register(
            "blastp",
            &args,
            &cols,
            TableCallback::new(blastp_handle),
            Some(&Extopts {
                description: Some(
                    "Protein BLAST (rust-bio Smith-Waterman + BLOSUM62 + Karlin-Altschul)".into(),
                ),
                tags: vec!["bioinformatics".into(), "alignment".into()],
            }),
        )?;

        Ok(Loadresult {
            name: "blast".into(),
            version: Some(env!("CARGO_PKG_VERSION").into()),
            requires: vec![Capabilitykind::Table],
        })
    }

    fn reconfigure(_keys: Vec<String>) -> Result<bool, Duckerror> {
        Ok(false)
    }

    fn shutdown() -> Result<bool, Duckerror> {
        Ok(false)
    }
}

fn seq_search_arg_types() -> Vec<Funcarg> {
    vec![
        Funcarg {
            name: Some("queries".into()),
            logical: Logicaltype::Text,
        },
        Funcarg {
            name: Some("subjects".into()),
            logical: Logicaltype::Text,
        },
        Funcarg {
            name: Some("options".into()),
            logical: Logicaltype::Text,
        },
    ]
}

fn hit_columns() -> Vec<Columndef> {
    let text = |name: &str| Columndef {
        name: name.into(),
        logical: Logicaltype::Text,
    };
    let u32c = |name: &str| Columndef {
        name: name.into(),
        logical: Logicaltype::Uint32,
    };
    let f64c = |name: &str| Columndef {
        name: name.into(),
        logical: Logicaltype::Float64,
    };
    vec![
        text("query_key"),
        text("subject_key"),
        u32c("query_start"),
        u32c("query_end"),
        u32c("subject_start"),
        u32c("subject_end"),
        text("strand"),
        u32c("identity_count"),
        u32c("alignment_length"),
        f64c("percent_identity"),
        f64c("bit_score"),
        f64c("raw_score"),
        f64c("evalue"),
    ]
}

// ---- callback-dispatch: real call_table, stubs for the rest -----------

impl callback_dispatch::Guest for Component {
    fn call_scalar_batch_col(
        _handle: u32,
        _args: Vec<callback_dispatch::Colvec>,
        _ctx: callback_dispatch::Invokeinfo,
    ) -> Result<callback_dispatch::Colvec, Duckerror> {
        Err(unsupported("no scalar functions"))
    }
    fn call_aggregate_col(
        _handle: u32,
        _args: Vec<callback_dispatch::Colvec>,
    ) -> Result<Duckvalue, Duckerror> {
        Err(unsupported("no aggregates"))
    }
    fn call_cast_col(
        _handle: u32,
        _arg: callback_dispatch::Colvec,
    ) -> Result<callback_dispatch::Colvec, Duckerror> {
        Err(unsupported("no casts"))
    }
    fn call_scalar(
        _handle: u32,
        _args: Vec<Duckvalue>,
        _ctx: callback_dispatch::Invokeinfo,
    ) -> Result<Duckvalue, Duckerror> {
        Err(unsupported("no scalar functions"))
    }
    fn call_table(
        handle: u32,
        args: Vec<Duckvalue>,
    ) -> Result<Resultset, Duckerror> {
        let handler = table_handlers()
            .lock()
            .unwrap()
            .get(&handle)
            .copied()
            .ok_or_else(|| {
                Duckerror::Internal(format!("blast: unknown table handle {handle}"))
            })?;

        let queries = parse_sequence_list(args.get(0), "queries")?;
        let subjects = parse_sequence_list(args.get(1), "subjects")?;
        let options = parse_options(args.get(2))?;

        let scoring = match handler {
            TableHandler::Blastn => Scoring::BlastnDefault,
            TableHandler::Blastp => Scoring::BlastpDefault,
        };

        // Pushdown fields on `options` (query_keys, subject_keys, strand)
        // are honored inside `run_search` — both the DuckLink dispatch
        // (here) and the biology-only `sequence-search::search` Guest see
        // the same filter semantics.
        let hits = crate::run_search(&queries, &subjects, &scoring, &options)
            .map_err(search_error_to_duck)?;

        Ok(hits.iter().map(hit_to_row).collect())
    }
    fn call_pragma(
        _handle: u32,
        _args: Vec<Duckvalue>,
    ) -> Result<Option<Duckvalue>, Duckerror> {
        Err(unsupported("no pragmas"))
    }
    fn call_cast(_handle: u32, _value: Duckvalue) -> Result<Duckvalue, Duckerror> {
        Err(unsupported("no casts"))
    }
}

fn unsupported(msg: &str) -> Duckerror {
    Duckerror::Unsupported(format!("blast: {msg}"))
}

// ---- streaming-table stubs: world exports the interface but no functions
// are registered on it, so every dispatch returns "not supported". Ready to
// hook back in when v5.x learns to dispatch streaming tables.

impl table_stream_dispatch::Guest for Component {
    fn call_table_open(
        _handle: u32,
        _args: Vec<Duckvalue>,
        _projection: Vec<u32>,
    ) -> Result<table_stream_dispatch::TableOpenResult, Duckerror> {
        Err(unsupported(
            "no streaming table functions registered; dispatch uses call_table",
        ))
    }
    fn call_table_open_filtered(
        _handle: u32,
        _args: Vec<Duckvalue>,
        _projection: Vec<u32>,
        _filters: Vec<table_stream_dispatch::TableFilter>,
    ) -> Result<table_stream_dispatch::TableOpenResult, Duckerror> {
        Err(unsupported(
            "no streaming table functions registered; dispatch uses call_table",
        ))
    }
    fn call_table_next(
        _handle: u32,
        _cursor: u32,
        _max_rows: u32,
    ) -> Result<table_stream_dispatch::Resultset, Duckerror> {
        Err(unsupported("no streaming cursors"))
    }
    fn call_table_close(_handle: u32, _cursor: u32) -> Result<bool, Duckerror> {
        Ok(false)
    }
}

// ---- row emission ------------------------------------------------------

fn hit_to_row(hit: &Hit) -> Vec<Duckvalue> {
    let strand_str = match hit.strand {
        Strand::Plus => "plus",
        Strand::Minus => "minus",
    };
    vec![
        Duckvalue::Text(hit.query_key.clone()),
        Duckvalue::Text(hit.subject_key.clone()),
        Duckvalue::Uint32(hit.query_start),
        Duckvalue::Uint32(hit.query_end),
        Duckvalue::Uint32(hit.subject_start),
        Duckvalue::Uint32(hit.subject_end),
        Duckvalue::Text(strand_str.into()),
        Duckvalue::Uint32(hit.identity_count),
        Duckvalue::Uint32(hit.alignment_length),
        Duckvalue::Float64(hit.percent_identity),
        Duckvalue::Float64(hit.bit_score),
        Duckvalue::Float64(hit.raw_score),
        Duckvalue::Float64(hit.evalue),
    ]
}

// ---- arg parsing ------------------------------------------------------

#[derive(Deserialize)]
struct RawSequence {
    key: String,
    data: String,
}

#[derive(Deserialize, Default)]
struct RawOptions {
    #[serde(default)]
    evalue_max: Option<f64>,
    #[serde(default)]
    max_target_seqs: Option<u32>,
    #[serde(default)]
    min_identity: Option<f64>,
    /// Whitelist of query keys — see the WIT `search-options.query-keys`
    /// docstring. `None` = no restriction; empty vec short-circuits the
    /// scan.
    #[serde(default)]
    query_keys: Option<Vec<String>>,
    /// Same shape and semantics as `query_keys`, applied to subjects.
    #[serde(default)]
    subject_keys: Option<Vec<String>>,
    /// One of `"plus"` or `"minus"` to restrict emitted hits by strand.
    /// Unknown labels are silently ignored (no filter applied).
    #[serde(default)]
    strand: Option<String>,
}

fn parse_sequence_list(v: Option<&Duckvalue>, name: &str) -> Result<Vec<Sequence>, Duckerror> {
    let json = match v {
        Some(Duckvalue::Text(s)) => s.as_str(),
        // Kept in case DuckLink starts preserving complex type-expressions.
        Some(Duckvalue::Complex(Complexvalue { json, .. })) => json.as_str(),
        _ => {
            return Err(Duckerror::Invalidargument(format!(
                "blast: '{name}' must be a JSON VARCHAR (or LIST(STRUCT(key,data)))"
            )));
        }
    };
    let raws: Vec<RawSequence> = serde_json::from_str(json).map_err(|e| {
        Duckerror::Invalidargument(format!("blast: cannot parse '{name}' JSON: {e}"))
    })?;
    Ok(raws
        .into_iter()
        .map(|r| Sequence {
            key: r.key,
            data: r.data,
        })
        .collect())
}

fn parse_options(v: Option<&Duckvalue>) -> Result<SearchOptions, Duckerror> {
    let raw: RawOptions = match v {
        None | Some(Duckvalue::Null) => RawOptions::default(),
        Some(Duckvalue::Text(s)) if s.is_empty() => RawOptions::default(),
        Some(Duckvalue::Text(s)) => serde_json::from_str(s).map_err(|e| {
            Duckerror::Invalidargument(format!("blast: cannot parse 'options' JSON: {e}"))
        })?,
        Some(Duckvalue::Complex(Complexvalue { json, .. })) => serde_json::from_str(json)
            .map_err(|e| {
                Duckerror::Invalidargument(format!("blast: cannot parse 'options' JSON: {e}"))
            })?,
        _ => {
            return Err(Duckerror::Invalidargument(
                "blast: 'options' must be a JSON VARCHAR (or NULL)".into(),
            ));
        }
    };
    Ok(SearchOptions {
        evalue_max: raw.evalue_max,
        max_target_seqs: raw.max_target_seqs,
        min_identity: raw.min_identity,
        query_keys: raw.query_keys,
        subject_keys: raw.subject_keys,
        strand: raw.strand,
    })
}

fn search_error_to_duck(e: crate::SearchError) -> Duckerror {
    use crate::SearchError::*;
    match e {
        EmptyQueries => Duckerror::Invalidargument("blast: queries is empty".into()),
        EmptySubjects => Duckerror::Invalidargument("blast: subjects is empty".into()),
        InvalidSequence(m) => Duckerror::Invalidargument(format!("blast: invalid sequence: {m}")),
        InvalidScoring(m) => Duckerror::Invalidargument(format!("blast: invalid scoring: {m}")),
        AlignmentFailed(m) => Duckerror::Internal(format!("blast: alignment failed: {m}")),
    }
}
