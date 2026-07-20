//! DuckLink dispatch surface for the genome-format component.
//!
//! Implements the three Guest traits of the `duckdb-extension-table-stream`
//! world: `guest` (lifecycle), `callback-dispatch` (scalar/aggregate/cast
//! stubs), `table-stream-dispatch` (the streaming table cursor protocol).
//!
//! One table function is registered during `load()`:
//! - `genbank_scan(contents VARCHAR)` — takes a single VARCHAR carrying
//!   the raw GenBank text (may contain multiple LOCUS records concatenated
//!   together — the parser already handles multi-record input), runs it
//!   through the same parser this component exports as
//!   `tegmentum:bio/genome-format`, and emits one row per feature with a
//!   "wide" projection that includes the qualifier map (as a JSON string)
//!   and the feature's extracted DNA sequence.
//!
//! ## Why `contents VARCHAR` instead of `paths VARCHAR`
//!
//! Wasm extensions run in the DuckLink loader's WASI sandbox, and the
//! `--dir HOST::GUEST` preopens registered on the CLI do NOT thread through
//! to loaded extension instances — `std::fs::read` from inside an extension
//! always returns `No such file or directory`, even for paths that DuckDB
//! itself can open. Rather than build a bespoke file-reading host import,
//! we let DuckDB do what it's already good at (I/O + globbing via
//! `read_text`) and take just the bytes.
//!
//! ## Intended composition
//!
//! Single file (composed via `read_text`, in a DuckDB build that permits
//! subqueries as TVF arguments):
//!
//! ```sql
//! SELECT *
//! FROM genbank_scan((SELECT content FROM read_text('lambda.gb')));
//! ```
//!
//! Multi-file glob (concatenate with a newline so the parser sees a
//! properly separated stream of LOCUS records):
//!
//! ```sql
//! SELECT *
//! FROM genbank_scan((
//!     SELECT string_agg(content, chr(10))
//!     FROM read_text('data/*.gb')
//! ));
//! ```
//!
//! A SQL macro is the natural next sugar to hide the sub-select, e.g.
//!
//! ```sql
//! CREATE OR REPLACE MACRO genbank_files(glob) AS TABLE
//!     SELECT * FROM genbank_scan((
//!         SELECT string_agg(content, chr(10)) FROM read_text(glob)
//!     ));
//! ```
//!
//! Note: DuckLink 4.0.0's vendored DuckDB binder currently rejects both
//! scalar subqueries (`Table function cannot contain subqueries`) and
//! LATERAL column references (`only supports literals as parameters`) in
//! TVF argument position. Until that lifts, callers who want to skip the
//! macro have to inline the file contents as a SQL string literal — which
//! is exactly what the acceptance harness does. `read_text` itself works
//! fine, so the moment the binder learns TVF subqueries the sub-select
//! form above becomes the daily-driver.
//!
//! Filter pushdown is not implemented on this first pass: the DuckLink core
//! re-checks every pushed filter above the scan, so correctness holds — we
//! simply forgo the optimisation. `call_table_open_filtered` reduces to
//! `call_table_open`.

use std::collections::HashMap;
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::{Mutex, OnceLock};

use crate::bindings::duckdb::extension::table_stream;
use crate::bindings::duckdb::extension::types::{
    Capabilitykind, Columndef, Complexvalue, Duckerror, Duckvalue, Funcarg, Loadresult,
    Logicaltype,
};
use crate::bindings::exports::duckdb::extension::{
    callback_dispatch, guest, table_stream_dispatch,
};
use crate::model::{Parsed, Qualifier};
use crate::{parser, Component};

const HANDLE_GENBANK_SCAN: u32 = 1;

/// One row of the wide feature projection `genbank_scan` emits. Kept in
/// memory between `call_table_open_filtered` and the paginating
/// `call_table_next` calls — Synphage-scale inputs are dozens of genomes
/// with thousands of features apiece, comfortably below the point where we
/// need a lazy per-file iterator.
struct Cursor {
    rows: Vec<Vec<Duckvalue>>,
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
        let args = genbank_scan_arg_types();
        let cols = feature_columns();

        table_stream::register_filterable_table(
            "genbank_scan",
            &args,
            &cols,
            HANDLE_GENBANK_SCAN,
        )?;

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

fn genbank_scan_arg_types() -> Vec<Funcarg> {
    // Single VARCHAR carrying the raw GenBank text. See the module doc for
    // why this isn't a path or a LIST(STRUCT(...)): WASI preopens don't
    // thread through to extensions, and DuckLink 4.0.0's Complex(...) type
    // erases to VARCHAR[] at the binder anyway.
    vec![Funcarg {
        name: Some("contents".into()),
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

fn feature_columns_projected(projection: &[u32]) -> Vec<Columndef> {
    let full = feature_columns();
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
        Err(unsupported("genome-format exports no scalar functions"))
    }
    fn call_aggregate_col(
        _handle: u32,
        _args: Vec<callback_dispatch::Colvec>,
    ) -> Result<Duckvalue, Duckerror> {
        Err(unsupported("genome-format exports no aggregates"))
    }
    fn call_cast_col(
        _handle: u32,
        _arg: callback_dispatch::Colvec,
    ) -> Result<callback_dispatch::Colvec, Duckerror> {
        Err(unsupported("genome-format exports no casts"))
    }
    fn call_scalar(
        _handle: u32,
        _args: Vec<Duckvalue>,
        _ctx: callback_dispatch::Invokeinfo,
    ) -> Result<Duckvalue, Duckerror> {
        Err(unsupported("genome-format exports no scalar functions"))
    }
    fn call_table(
        _handle: u32,
        _args: Vec<Duckvalue>,
    ) -> Result<callback_dispatch::Resultset, Duckerror> {
        Err(unsupported(
            "genome-format uses the streaming table dispatch, not the row-major one",
        ))
    }
    fn call_pragma(
        _handle: u32,
        _args: Vec<Duckvalue>,
    ) -> Result<Option<Duckvalue>, Duckerror> {
        Err(unsupported("genome-format exports no pragmas"))
    }
    fn call_cast(_handle: u32, _value: Duckvalue) -> Result<Duckvalue, Duckerror> {
        Err(unsupported("genome-format exports no casts"))
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
        open_inner(handle, args, projection)
    }

    fn call_table_open_filtered(
        handle: u32,
        args: Vec<Duckvalue>,
        projection: Vec<u32>,
        _filters: Vec<table_stream_dispatch::TableFilter>,
    ) -> Result<table_stream_dispatch::TableOpenResult, Duckerror> {
        // Filter pushdown is ignored on this first pass. Per the DuckLink
        // freeze policy the core re-checks every pushed filter above the
        // scan, so correctness holds — we just forgo the optimisation.
        open_inner(handle, args, projection)
    }

    fn call_table_next(
        _handle: u32,
        cursor: u32,
        max_rows: u32,
    ) -> Result<table_stream_dispatch::Resultset, Duckerror> {
        let mut guard = cursors().lock().unwrap();
        let cur = guard.get_mut(&cursor).ok_or_else(|| {
            Duckerror::Invalidstate(format!("genome-format: unknown cursor {cursor}"))
        })?;

        let mut out: Vec<Vec<Duckvalue>> = Vec::new();
        while cur.next < cur.rows.len() && (out.len() as u32) < max_rows {
            let row = cur.rows[cur.next].clone();
            cur.next += 1;
            out.push(project(row, &cur.projection));
        }
        Ok(out)
    }

    fn call_table_close(_handle: u32, cursor: u32) -> Result<bool, Duckerror> {
        Ok(cursors().lock().unwrap().remove(&cursor).is_some())
    }
}

fn open_inner(
    handle: u32,
    args: Vec<Duckvalue>,
    projection: Vec<u32>,
) -> Result<table_stream_dispatch::TableOpenResult, Duckerror> {
    if handle != HANDLE_GENBANK_SCAN {
        return Err(Duckerror::Invalidstate(format!(
            "genome-format: unknown callback handle {handle}"
        )));
    }

    let contents = parse_contents_arg(args.first())?;
    let parsed = parser::parse_genbank(contents.as_bytes()).map_err(|e| {
        let reason = match e {
            crate::model::ParseError::Malformed(m) => m,
            crate::model::ParseError::UnsupportedVersion(m) => {
                format!("unsupported version: {m}")
            }
        };
        Duckerror::Invalidargument(format!("genbank_scan: cannot parse contents: {reason}"))
    })?;

    let mut rows: Vec<Vec<Duckvalue>> = Vec::new();
    emit_rows(&parsed, &mut rows);

    let id = next_cursor_id();
    let projected_cols = feature_columns_projected(&projection);
    cursors().lock().unwrap().insert(
        id,
        Cursor {
            rows,
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
