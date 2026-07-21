//! DuckLink dispatch surface for the genome-format component.
//!
//! Registers two table functions plus a replacement scan, all routed through
//! `runtime::TableRegistry` (the row-major `callback_dispatch::call_table`
//! path — the only registration path the `.gb` / `.gbk` replacement scan can
//! reference, since the replacement-scan handle registry looks up names via
//! the runtime-side `table_handle_names` populated only by that call).
//!
//! - `genbank_scan(contents VARCHAR)` — parses raw GenBank text passed as a
//!   VARCHAR. Multi-record files (multiple LOCUS blocks separated by `//`)
//!   are handled by the parser.
//! - `genbank_read_path(path VARCHAR)` — reads the file at `path` from
//!   inside the extension via `std::fs::read`, then hands the bytes to the
//!   same parser.
//! - **Replacement scan** on file extensions `.gb` / `.gbk` (mode:
//!   ExtensionOnly): `SELECT * FROM 'lambda.gb'` rewrites to
//!   `genbank_read_path('lambda.gb')` at parse time — sidesteps the
//!   DuckDB TVF-subquery binder rule entirely.
//!
//! ## Why row-major dispatch (call_table) instead of streaming
//!
//! The replacement-scan mechanism looks up its target table function via
//! `runtime`'s `table_handle_names` registry, which is populated only by
//! `runtime::TableRegistry::register`. Table functions registered through
//! the alternative `table_stream::register_filterable_table` path (which
//! we previously used) land in a different registry entirely — the
//! replacement scan can't reach them. So both `genbank_scan` and
//! `genbank_read_path` are registered on the runtime path and dispatched
//! row-major via `callback_dispatch::call_table`. Filter pushdown and
//! streaming are unavailable on this path; for GenBank workloads (bounded
//! by file size, always materialised in one pass anyway) that is fine.
//!
//! ## Composing with DuckDB's read_text
//!
//! With the replacement scan wired, the daily-driver invocation is:
//!
//! ```sql
//! SELECT * FROM 'lambda.gb';
//! SELECT * FROM 'phages/*.gb';   -- DuckDB expands the glob before matching
//! ```
//!
//! When file access is unavailable (e.g. the standalone DuckLink CLI does
//! not thread `--dir` preopens through to loaded extensions), users fall
//! back to `genbank_scan(<inlined-content>)` or, in DuckDB builds that
//! permit TVF-arg subqueries, `genbank_scan((SELECT content FROM
//! read_text('lambda.gb')))`.

use std::collections::HashMap;
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::{Mutex, OnceLock};

use crate::bindings::duckdb::extension::files::{
    self, DetectionMode, ReplacementScan,
};
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
use crate::model::{Parsed, Qualifier};
use crate::{parser, Component};

/// Discriminator stored in `TABLE_HANDLERS` keyed by the callback handle we
/// mint in `load()` and pass to `TableCallback::new`. The runtime threads
/// that same handle back into `callback_dispatch::call_table`, and we use it
/// to pick which parser front-end to invoke.
#[derive(Copy, Clone)]
enum TableHandler {
    /// `genbank_scan(contents VARCHAR)` — parse a VARCHAR of GenBank text.
    Scan,
    /// `genbank_read_path(path VARCHAR)` — read the file at `path`, parse.
    ReadPath,
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
            Duckerror::Internal("genome-format: host did not expose table capability".into())
        })?;
        let registry = match capability {
            Capability::Table(r) => r,
            _ => {
                return Err(Duckerror::Internal(
                    "genome-format: table capability returned unexpected variant".into(),
                ))
            }
        };
        let cols = feature_columns();

        // 1) genbank_scan(contents VARCHAR).
        let scan_handle = next_handle();
        table_handlers()
            .lock()
            .unwrap()
            .insert(scan_handle, TableHandler::Scan);
        registry.register(
            "genbank_scan",
            &contents_arg_types(),
            &cols,
            TableCallback::new(scan_handle),
            Some(&Extopts {
                description: Some(
                    "Parse a VARCHAR of raw GenBank text into wide feature rows".into(),
                ),
                tags: vec!["genomics".into(), "genbank".into()],
            }),
        )?;

        // 2) genbank_read_path(path VARCHAR) — captured for the replacement scan.
        let read_path_handle = next_handle();
        table_handlers()
            .lock()
            .unwrap()
            .insert(read_path_handle, TableHandler::ReadPath);
        let read_path_reg_id = registry.register(
            "genbank_read_path",
            &path_arg_types(),
            &cols,
            TableCallback::new(read_path_handle),
            Some(&Extopts {
                description: Some(
                    "Read a GenBank file at the given path and emit wide feature rows".into(),
                ),
                tags: vec!["genomics".into(), "genbank".into()],
            }),
        )?;

        // 3) Wire the replacement scan: `SELECT * FROM 'lambda.gb'` -> the
        //    read-path function. ExtensionOnly matches on the string suffix;
        //    Signature would peek at the file bytes for detection.
        files::register_replacement_scan(&ReplacementScan {
            extensions: vec!["gb".into(), "gbk".into()],
            table_function: read_path_reg_id,
            mode: DetectionMode::ExtensionOnly,
        })
        .map_err(|e| {
            Duckerror::Internal(format!(
                "genome-format: register_replacement_scan(gb, gbk): {e}"
            ))
        })?;

        Ok(Loadresult {
            name: "genome-format".into(),
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

fn contents_arg_types() -> Vec<Funcarg> {
    vec![Funcarg {
        name: Some("contents".into()),
        logical: Logicaltype::Text,
    }]
}

fn path_arg_types() -> Vec<Funcarg> {
    vec![Funcarg {
        name: Some("path".into()),
        logical: Logicaltype::Text,
    }]
}

fn feature_columns() -> Vec<Columndef> {
    let text = |name: &str| Columndef {
        name: name.into(),
        logical: Logicaltype::Text,
    };
    let u32c = |name: &str| Columndef {
        name: name.into(),
        logical: Logicaltype::Uint32,
    };
    let i32c = |name: &str| Columndef {
        name: name.into(),
        logical: Logicaltype::Int32,
    };
    vec![
        text("record_id"),
        text("accession"),
        text("version"),
        text("organism"),
        u32c("record_length"),
        u32c("feature_index"),
        text("feature_type"),
        u32c("start_position"),
        u32c("end_position"),
        i32c("strand"),
        text("qualifiers_json"),
        text("sequence"),
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
                Duckerror::Internal(format!("genome-format: unknown table handle {handle}"))
            })?;

        let bytes: Vec<u8> = match handler {
            TableHandler::Scan => parse_contents_arg(args.first())?.into_bytes(),
            TableHandler::ReadPath => {
                let path = parse_path_arg(args.first())?;
                std::fs::read(&path).map_err(|e| {
                    Duckerror::Io(format!(
                        "genbank_read_path: cannot read '{path}': {e}. \
                         If this DuckLink build does not grant filesystem access to \
                         extensions, read the file with DuckDB's read_text and call \
                         genbank_scan(<contents>) instead."
                    ))
                })?
            }
        };

        let parsed = parser::parse_genbank(&bytes).map_err(|e| {
            let reason = match e {
                crate::model::ParseError::Malformed(m) => m,
                crate::model::ParseError::UnsupportedVersion(m) => {
                    format!("unsupported version: {m}")
                }
            };
            Duckerror::Invalidargument(format!("genome-format: cannot parse GenBank: {reason}"))
        })?;

        let mut rows: Vec<Vec<Duckvalue>> = Vec::new();
        emit_rows(&parsed, &mut rows);
        Ok(rows)
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
    Duckerror::Unsupported(format!("genome-format: {msg}"))
}

// ---- streaming-table stubs: world exports this interface but we register
// nothing here, so every dispatch is a "no such handle" error.

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

// ---- argument parsing --------------------------------------------------

fn parse_contents_arg(v: Option<&Duckvalue>) -> Result<String, Duckerror> {
    match v {
        Some(Duckvalue::Text(s)) => Ok(s.clone()),
        // Kept in case DuckLink starts preserving complex type-expressions —
        // parallels the escape hatch in components/blast/src/ducklink.rs.
        Some(Duckvalue::Complex(Complexvalue { json, .. })) => Ok(json.clone()),
        _ => Err(Duckerror::Invalidargument(
            "genbank_scan: 'contents' must be a VARCHAR carrying raw GenBank text".into(),
        )),
    }
}

fn parse_path_arg(v: Option<&Duckvalue>) -> Result<String, Duckerror> {
    match v {
        Some(Duckvalue::Text(s)) => Ok(s.clone()),
        _ => Err(Duckerror::Invalidargument(
            "genbank_read_path: 'path' must be a VARCHAR filesystem path".into(),
        )),
    }
}

// ---- row emission ------------------------------------------------------

/// Group qualifiers by (record_id, feature_index) so `emit_rows` can build
/// each feature's JSON object in one pass instead of scanning the full
/// qualifier vector per feature.
fn index_qualifiers(qs: &[Qualifier]) -> HashMap<(String, u32), Vec<&Qualifier>> {
    let mut m: HashMap<(String, u32), Vec<&Qualifier>> = HashMap::new();
    for q in qs {
        m.entry((q.record_id.clone(), q.feature_index))
            .or_default()
            .push(q);
    }
    m
}

/// Build the JSON qualifier object for one feature. Multiple values under
/// the same name collapse to the FIRST occurrence — full multi-value
/// support is a later slice. Empty object `{}` if no qualifiers.
///
/// Serialised by hand so we don't have to spin up a serde_json::Map + Value
/// tree for a small object we're only going to serialise once. The key /
/// value escaping still goes through serde_json to keep quoting correct.
fn qualifiers_to_json(qs: &[&Qualifier]) -> String {
    if qs.is_empty() {
        return "{}".to_string();
    }
    let mut seen: HashMap<&str, ()> = HashMap::new();
    let mut buf = String::from("{");
    let mut first = true;
    for q in qs {
        if seen.contains_key(q.name.as_str()) {
            continue;
        }
        seen.insert(q.name.as_str(), ());
        if !first {
            buf.push(',');
        }
        first = false;
        buf.push_str(&serde_json::to_string(&q.name).unwrap_or_else(|_| "\"\"".into()));
        buf.push(':');
        buf.push_str(&serde_json::to_string(&q.value).unwrap_or_else(|_| "\"\"".into()));
    }
    buf.push('}');
    buf
}

/// Turn one `Parsed` bundle into wide feature rows and append them to
/// `rows`. Every feature contributes one row; the enclosing record's
/// metadata (`accession`, `organism`, `record_length`) is denormalised on.
fn emit_rows(parsed: &Parsed, rows: &mut Vec<Vec<Duckvalue>>) {
    // Records may share record_ids across a multi-record file only if the
    // producer duplicated them — treat the first occurrence as authoritative
    // and let SQL callers deduplicate downstream if needed.
    let mut record_by_id: HashMap<&str, &crate::model::Record> = HashMap::new();
    for r in &parsed.records {
        record_by_id.entry(r.record_id.as_str()).or_insert(r);
    }

    let mut sequence_by_id: HashMap<&str, &str> = HashMap::new();
    for s in &parsed.sequences {
        sequence_by_id
            .entry(s.record_id.as_str())
            .or_insert(s.data.as_str());
    }

    let qual_ix = index_qualifiers(&parsed.qualifiers);

    for feat in &parsed.features {
        let rec = record_by_id.get(feat.record_id.as_str());
        let (accession, version, organism, record_length) = match rec {
            Some(r) => (
                r.accession.clone(),
                r.version.clone(),
                r.organism.clone(),
                r.length,
            ),
            None => (String::new(), String::new(), String::new(), 0u32),
        };

        let sequence = sequence_by_id
            .get(feat.record_id.as_str())
            .map(|s| {
                extract_feature_sequence(s, feat.start_position, feat.end_position, feat.strand)
            })
            .unwrap_or_default();

        // Skip our synthetic `_location` qualifier — it was emitted by the
        // parser as an escape hatch for SQL callers that want the raw
        // location string, but it isn't a real GenBank qualifier and would
        // pollute every JSON object here. Users who need it can query the
        // biology-only relations directly.
        let quals: Vec<&Qualifier> = qual_ix
            .get(&(feat.record_id.clone(), feat.feature_index))
            .map(|v| v.iter().copied().filter(|q| q.name != "_location").collect())
            .unwrap_or_default();
        let qualifiers_json = qualifiers_to_json(&quals);

        rows.push(vec![
            Duckvalue::Text(feat.record_id.clone()),
            Duckvalue::Text(accession),
            Duckvalue::Text(version),
            Duckvalue::Text(organism),
            Duckvalue::Uint32(record_length),
            Duckvalue::Uint32(feat.feature_index),
            Duckvalue::Text(feat.feature_type.clone()),
            Duckvalue::Uint32(feat.start_position),
            Duckvalue::Uint32(feat.end_position),
            Duckvalue::Int32(feat.strand as i32),
            Duckvalue::Text(qualifiers_json),
            Duckvalue::Text(sequence),
        ]);
    }
}

// ---- sequence extraction ----------------------------------------------

/// Extract the DNA slice for one feature. 1-indexed inclusive positions
/// (matching NCBI convention); returns the empty string when the interval
/// is out of bounds or when a fuzzy / joined location collapsed to (0, 0)
/// during parsing — the empty string is a signal to the SQL caller that
/// they should reach for the raw `/_location` qualifier if they need the
/// exact structure.
fn extract_feature_sequence(seq: &str, start: u32, end: u32, strand: i8) -> String {
    let record_len = seq.len() as u32;
    if start < 1 || end < start || end > record_len {
        return String::new();
    }
    let s = (start - 1) as usize;
    let e = end as usize;
    let slice = &seq.as_bytes()[s..e];
    if strand == -1 {
        String::from_utf8(revcomp(slice)).unwrap_or_default()
    } else {
        String::from_utf8(slice.to_vec()).unwrap_or_default()
    }
}

/// Nucleotide reverse-complement over the ambiguity-tolerant alphabet
/// (matches `components/blast/src/strand.rs` — kept local so this crate
/// doesn't depend on the sibling crate for eight lines of biology).
fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'a' => b't',
        b'T' => b'A',
        b't' => b'a',
        b'C' => b'G',
        b'c' => b'g',
        b'G' => b'C',
        b'g' => b'c',
        b'U' => b'A',
        b'u' => b'a',
        b'N' | b'n' => b,
        _ => b'N',
    }
}

// ---- tests --------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extract_forward_strand() {
        // 1-indexed inclusive: positions 2..5 of ACGTAC -> "CGTA".
        assert_eq!(extract_feature_sequence("ACGTAC", 2, 5, 1), "CGTA");
    }

    #[test]
    fn extract_reverse_strand_revcomps() {
        // Slice CGTA -> revcomp TACG.
        assert_eq!(extract_feature_sequence("ACGTAC", 2, 5, -1), "TACG");
    }

    #[test]
    fn extract_out_of_bounds_is_empty() {
        assert_eq!(extract_feature_sequence("ACGT", 1, 10, 1), "");
        assert_eq!(extract_feature_sequence("ACGT", 0, 3, 1), "");
        assert_eq!(extract_feature_sequence("ACGT", 3, 2, 1), "");
    }

    #[test]
    fn qualifiers_json_first_value_wins() {
        // gene appears twice; the first value wins per this component's
        // first-slice rule.
        let a = Qualifier {
            record_id: "r".into(),
            feature_index: 0,
            name: "gene".into(),
            value: "thrL".into(),
        };
        let b = Qualifier {
            record_id: "r".into(),
            feature_index: 0,
            name: "gene".into(),
            value: "thrL_alt".into(),
        };
        let c = Qualifier {
            record_id: "r".into(),
            feature_index: 0,
            name: "product".into(),
            value: "thr operon leader".into(),
        };
        let json = qualifiers_to_json(&[&a, &b, &c]);
        assert!(json.contains("\"gene\":\"thrL\""));
        assert!(!json.contains("thrL_alt"));
        assert!(json.contains("\"product\":\"thr operon leader\""));
    }

    #[test]
    fn qualifiers_json_empty() {
        assert_eq!(qualifiers_to_json(&[]), "{}");
    }

    #[test]
    fn parse_contents_arg_accepts_text() {
        let v = Duckvalue::Text("LOCUS x".into());
        assert_eq!(parse_contents_arg(Some(&v)).unwrap(), "LOCUS x");
    }

    #[test]
    fn parse_contents_arg_rejects_null() {
        assert!(parse_contents_arg(Some(&Duckvalue::Null)).is_err());
    }
}
