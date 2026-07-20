//! DuckLink dispatch surface for the BLAST component.
//!
//! Implements the three Guest traits of the `duckdb-extension-table-stream`
//! world: `guest` (lifecycle), `callback-dispatch` (scalar/aggregate/cast
//! stubs), `table-stream-dispatch` (the streaming table cursor protocol).
//!
//! Two table functions are registered during `load()`:
//! - `blastn(queries LIST(STRUCT(key VARCHAR, data VARCHAR)),
//!            subjects LIST(STRUCT(key VARCHAR, data VARCHAR)),
//!            options STRUCT(evalue_max DOUBLE, max_target_seqs INTEGER,
//!                           min_identity DOUBLE))`
//! - `blastp(...)` with the same signature.
//!
//! Both delegate to `crate::run_search` with a fixed scoring preset
//! (`BlastnDefault` / `BlastpDefault`). Tunable scoring is future work —
//! the SQL surface stays clean for the common case first.

use std::collections::HashMap;
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::{Mutex, OnceLock};

use serde::Deserialize;

use crate::bindings::duckdb::extension::table_stream;
use crate::bindings::duckdb::extension::types::{
    Capabilitykind, Columndef, Complexvalue, Duckerror, Duckvalue, Funcarg, Loadresult, Logicaltype,
};
use crate::bindings::exports::duckdb::extension::{
    callback_dispatch, guest, table_stream_dispatch,
};
use crate::pushdown::{self, PushdownPlan};
use crate::{Component, Hit, Scoring, SearchOptions, Sequence, Strand};

const HANDLE_BLASTN: u32 = 1;
const HANDLE_BLASTP: u32 = 2;

/// Type-expression string carried in the `complex(...)` logicaltype for the
/// two sequence-list args. DuckLink echoes it back verbatim as the
/// `type-expr` field of the incoming `complexvalue`; we ignore that echo
/// and parse the `json` payload directly.
const TYPE_EXPR_SEQ_LIST: &str = "LIST(STRUCT(key VARCHAR, data VARCHAR))";
const TYPE_EXPR_OPTIONS: &str =
    "STRUCT(evalue_max DOUBLE, max_target_seqs INTEGER, min_identity DOUBLE)";

/// State per open scan. `hits` is materialised eagerly by `run_search`
/// during `call_table_open_filtered`; the cursor merely paginates it in
/// `max-rows`-sized batches through `call_table_next`.
struct Cursor {
    hits: Vec<Hit>,
    next: usize,
    projection: Vec<u32>,
}

fn cursors() -> &'static Mutex<HashMap<u32, Cursor>> {
    static C: OnceLock<Mutex<HashMap<u32, Cursor>>> = OnceLock::new();
    C.get_or_init(|| Mutex::new(HashMap::new()))
}

fn next_cursor_id() -> u32 {
    static N: AtomicU32 = AtomicU32::new(1);
    N.fetch_add(1, Ordering::Relaxed)
}

// ---- guest lifecycle ---------------------------------------------------

impl guest::Guest for Component {
    fn load() -> Result<Loadresult, Duckerror> {
        let args = seq_search_arg_types();
        let cols = hit_columns();

        table_stream::register_filterable_table("blastn", &args, &cols, HANDLE_BLASTN)?;
        table_stream::register_filterable_table("blastp", &args, &cols, HANDLE_BLASTP)?;

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
    // NOTE: The natural SQL surface would be LIST(STRUCT(key, data)) via
    // `Logicaltype::Complex(TYPE_EXPR_SEQ_LIST)`, but DuckLink 4.0.0 flattens
    // the type-expression string down to `VARCHAR[]` at the DuckDB binder so
    // the STRUCT payload is lost and duckvalues arrive as a plain array of
    // strings instead of `complex(json)`. Until DuckLink learns to preserve
    // arbitrary type-expressions, we take three VARCHAR args carrying JSON
    // payloads. Consumers can wrap with a DuckDB macro that hides the
    // to_json() calls; see sql/blast_macros.sql (future).
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

fn hit_columns_projected(projection: &[u32]) -> Vec<Columndef> {
    let full = hit_columns();
    if projection.is_empty() {
        full
    } else {
        projection
            .iter()
            .filter_map(|&i| full.get(i as usize).cloned())
            .collect()
    }
}

// ---- callback-dispatch stubs (table-only extension) -------------------

impl callback_dispatch::Guest for Component {
    fn call_scalar_batch_col(
        _handle: u32,
        _args: Vec<callback_dispatch::Colvec>,
        _ctx: callback_dispatch::Invokeinfo,
    ) -> Result<callback_dispatch::Colvec, Duckerror> {
        Err(unsupported("blast exports no scalar functions"))
    }
    fn call_aggregate_col(
        _handle: u32,
        _args: Vec<callback_dispatch::Colvec>,
    ) -> Result<Duckvalue, Duckerror> {
        Err(unsupported("blast exports no aggregates"))
    }
    fn call_cast_col(
        _handle: u32,
        _arg: callback_dispatch::Colvec,
    ) -> Result<callback_dispatch::Colvec, Duckerror> {
        Err(unsupported("blast exports no casts"))
    }
    fn call_scalar(
        _handle: u32,
        _args: Vec<Duckvalue>,
        _ctx: callback_dispatch::Invokeinfo,
    ) -> Result<Duckvalue, Duckerror> {
        Err(unsupported("blast exports no scalar functions"))
    }
    fn call_table(
        _handle: u32,
        _args: Vec<Duckvalue>,
    ) -> Result<callback_dispatch::Resultset, Duckerror> {
        Err(unsupported(
            "blast uses the streaming table dispatch, not the row-major one",
        ))
    }
    fn call_pragma(
        _handle: u32,
        _args: Vec<Duckvalue>,
    ) -> Result<Option<Duckvalue>, Duckerror> {
        Err(unsupported("blast exports no pragmas"))
    }
    fn call_cast(
        _handle: u32,
        _value: Duckvalue,
    ) -> Result<Duckvalue, Duckerror> {
        Err(unsupported("blast exports no casts"))
    }
}

fn unsupported(msg: &str) -> Duckerror {
    Duckerror::Unsupported(msg.into())
}

// ---- streaming table dispatch -----------------------------------------

impl table_stream_dispatch::Guest for Component {
    fn call_table_open(
        handle: u32,
        args: Vec<Duckvalue>,
        projection: Vec<u32>,
    ) -> Result<table_stream_dispatch::TableOpenResult, Duckerror> {
        open_inner(handle, args, projection, PushdownPlan::default())
    }

    fn call_table_open_filtered(
        handle: u32,
        args: Vec<Duckvalue>,
        projection: Vec<u32>,
        filters: Vec<table_stream_dispatch::TableFilter>,
    ) -> Result<table_stream_dispatch::TableOpenResult, Duckerror> {
        // Filter pushdown: recognise the handful of clauses we can profitably
        // absorb (evalue ceiling, identity floor, key restrictions, strand)
        // and let the DuckLink core re-check anything unrecognised above the
        // scan — that's still correct per the freeze policy, just less
        // efficient. See `pushdown::plan` for the full recognition set.
        let plan = pushdown::plan(&filters);
        open_inner(handle, args, projection, plan)
    }

    fn call_table_next(
        _handle: u32,
        cursor: u32,
        max_rows: u32,
    ) -> Result<table_stream_dispatch::Resultset, Duckerror> {
        let mut guard = cursors().lock().unwrap();
        let cur = guard.get_mut(&cursor).ok_or_else(|| {
            Duckerror::Invalidstate(format!("blast: unknown cursor {cursor}"))
        })?;

        let mut rows: Vec<Vec<Duckvalue>> = Vec::new();
        while cur.next < cur.hits.len() && (rows.len() as u32) < max_rows {
            let hit = &cur.hits[cur.next];
            cur.next += 1;
            rows.push(project(hit_to_row(hit), &cur.projection));
        }
        Ok(rows)
    }

    fn call_table_close(_handle: u32, cursor: u32) -> Result<bool, Duckerror> {
        Ok(cursors().lock().unwrap().remove(&cursor).is_some())
    }
}

fn open_inner(
    handle: u32,
    args: Vec<Duckvalue>,
    projection: Vec<u32>,
    plan: PushdownPlan,
) -> Result<table_stream_dispatch::TableOpenResult, Duckerror> {
    let mut queries = parse_sequence_list(args.get(0), "queries")?;
    let mut subjects = parse_sequence_list(args.get(1), "subjects")?;
    let mut options = parse_options(args.get(2))?;

    let scoring = match handle {
        HANDLE_BLASTN => Scoring::BlastnDefault,
        HANDLE_BLASTP => Scoring::BlastpDefault,
        other => {
            return Err(Duckerror::Invalidstate(format!(
                "blast: unknown callback handle {other}"
            )));
        }
    };

    let hits = if plan.short_circuit {
        // The pushed filter set is unsatisfiable (e.g. contradictory key
        // clauses). Open an empty cursor and skip alignment entirely.
        Vec::new()
    } else {
        pushdown::tighten_options(&plan, &mut options);
        let batches_empty = pushdown::prune_batches(&plan, &mut queries, &mut subjects);
        if batches_empty {
            Vec::new()
        } else {
            let raw = crate::run_search(&queries, &subjects, &scoring, &options)
                .map_err(search_error_to_duck)?;
            // The scoring params ran both strands for BLASTN; a
            // `WHERE strand = 'plus'|'minus'` clause is honoured by
            // dropping the non-matching orientation post-alignment. This
            // wastes minus-strand work when strand='plus' is pinned, but
            // stays contained and correct — the alternative would need to
            // reach into `scoring::resolve` from the dispatch layer.
            if let Some(sf) = plan.strand_keep {
                raw.into_iter().filter(|h| sf.matches(h.strand)).collect()
            } else {
                raw
            }
        }
    };

    let id = next_cursor_id();
    let projected_cols = hit_columns_projected(&projection);
    cursors().lock().unwrap().insert(
        id,
        Cursor {
            hits,
            next: 0,
            projection,
        },
    );

    Ok(table_stream_dispatch::TableOpenResult {
        cursor: id,
        columns: projected_cols,
    })
}

fn project(row: Vec<Duckvalue>, projection: &[u32]) -> Vec<Duckvalue> {
    if projection.is_empty() {
        row
    } else {
        projection
            .iter()
            .map(|&i| row.get(i as usize).cloned().unwrap_or(Duckvalue::Null))
            .collect()
    }
}

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
