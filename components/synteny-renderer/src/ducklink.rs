//! DuckLink dispatch surface for the synteny-renderer component.
//!
//! Implements the three Guest traits of the compound `synteny-renderer`
//! world: `guest` (lifecycle), `callback-dispatch` (scalar/aggregate/cast
//! stubs -- all Unsupported), `table-stream-dispatch` (the streaming table
//! cursor protocol).
//!
//! One table function is registered during `load()`:
//!
//! ```text
//! render_synteny_svg(
//!     tracks   VARCHAR,   -- JSON: [{"track_id":..., "label":..., "length":...}, ...]
//!     features VARCHAR,   -- JSON: [{"track_id":..., "feature_id":...,
//!                         --         "start_position":..., "end_position":...,
//!                         --         "strand":..., "colour":..., "label":...}, ...]
//!     links    VARCHAR    -- JSON: [{"query_track":..., "query_feature":...,
//!                         --         "subject_track":..., "subject_feature":...,
//!                         --         "identity":..., "colour":...}, ...]
//! ) -> TABLE(svg BLOB, bytes_len UINTEGER)
//! ```
//!
//! Example:
//!
//! ```sql
//! SELECT bytes_len FROM render_synteny_svg(
//!     '[{"track_id":"t1","label":"genome A","length":1000}]',
//!     '[{"track_id":"t1","feature_id":"g1","start_position":100,
//!        "end_position":300,"strand":1,"colour":null,"label":"geneA"}]',
//!     '[]');
//! ```
//!
//! Args are JSON VARCHAR strings rather than the natural LIST(STRUCT(...))
//! shape because DuckLink 4.0.0 flattens type-expressions down to VARCHAR[]
//! at the DuckDB binder -- the STRUCT payload is lost. Once DuckLink learns
//! to preserve arbitrary type-expressions, the parsers already accept the
//! `complex(json)` fallback path and only the Funcarg types need lifting.
//!
//! The cursor emits exactly one row (the rendered SVG + its byte length)
//! on the first `call_table_next`, then empty on all subsequent calls.
//! If the underlying renderer reports `EmptyInput`, the cursor opens with
//! no payload and the very first `call_table_next` returns zero rows -- so
//! `SELECT * FROM render_synteny_svg(NULL, NULL, NULL)` returns an empty
//! resultset instead of a SQL exception.

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
use crate::render::{self, Feature, Link, RenderError, Track};
use crate::Component;

const HANDLE_RENDER: u32 = 1;

/// State per open scan. Rendering runs up-front in `call_table_open_*`; the
/// cursor just holds the resulting bytes and hands them out once.
///
/// `payload = None` marks an empty-input open (RenderError::EmptyInput). In
/// that case `already_emitted` stays true throughout so `call_table_next`
/// always returns zero rows -- the empty resultset the task requires.
struct Cursor {
    payload: Option<Vec<u8>>,
    already_emitted: bool,
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
        let args = render_arg_types();
        let cols = svg_columns();

        table_stream::register_filterable_table(
            "render_synteny_svg",
            &args,
            &cols,
            HANDLE_RENDER,
        )?;

        Ok(Loadresult {
            name: "synteny-renderer".into(),
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

fn render_arg_types() -> Vec<Funcarg> {
    // NOTE: The natural SQL surface would be LIST(STRUCT(...)) via
    // `Logicaltype::Complex(...)`, but DuckLink 4.0.0 flattens the
    // type-expression string down to `VARCHAR[]` at the DuckDB binder so the
    // STRUCT payload is lost and duckvalues arrive as a plain array of
    // strings instead of `complex(json)`. Until DuckLink learns to preserve
    // arbitrary type-expressions, we take three VARCHAR args carrying JSON
    // payloads. Same workaround as the blast component.
    vec![
        Funcarg {
            name: Some("tracks".into()),
            logical: Logicaltype::Text,
        },
        Funcarg {
            name: Some("features".into()),
            logical: Logicaltype::Text,
        },
        Funcarg {
            name: Some("links".into()),
            logical: Logicaltype::Text,
        },
    ]
}

fn svg_columns() -> Vec<Columndef> {
    vec![
        Columndef {
            name: "svg".into(),
            logical: Logicaltype::Blob,
        },
        Columndef {
            name: "bytes_len".into(),
            logical: Logicaltype::Uint32,
        },
    ]
}

fn svg_columns_projected(projection: &[u32]) -> Vec<Columndef> {
    let full = svg_columns();
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
        Err(unsupported("synteny-renderer exports no scalar functions"))
    }
    fn call_aggregate_col(
        _handle: u32,
        _args: Vec<callback_dispatch::Colvec>,
    ) -> Result<Duckvalue, Duckerror> {
        Err(unsupported("synteny-renderer exports no aggregates"))
    }
    fn call_cast_col(
        _handle: u32,
        _arg: callback_dispatch::Colvec,
    ) -> Result<callback_dispatch::Colvec, Duckerror> {
        Err(unsupported("synteny-renderer exports no casts"))
    }
    fn call_scalar(
        _handle: u32,
        _args: Vec<Duckvalue>,
        _ctx: callback_dispatch::Invokeinfo,
    ) -> Result<Duckvalue, Duckerror> {
        Err(unsupported("synteny-renderer exports no scalar functions"))
    }
    fn call_table(
        _handle: u32,
        _args: Vec<Duckvalue>,
    ) -> Result<callback_dispatch::Resultset, Duckerror> {
        Err(unsupported(
            "synteny-renderer uses the streaming table dispatch, not the row-major one",
        ))
    }
    fn call_pragma(
        _handle: u32,
        _args: Vec<Duckvalue>,
    ) -> Result<Option<Duckvalue>, Duckerror> {
        Err(unsupported("synteny-renderer exports no pragmas"))
    }
    fn call_cast(_handle: u32, _value: Duckvalue) -> Result<Duckvalue, Duckerror> {
        Err(unsupported("synteny-renderer exports no casts"))
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
        // Filter pushdown ignored -- correctness holds because the DuckLink
        // core re-checks any pushed filters above the scan. Only two
        // columns (svg, bytes_len) come out per open anyway, so there is
        // nothing meaningful to prune at this layer.
        open_inner(handle, args, projection)
    }

    fn call_table_next(
        _handle: u32,
        cursor: u32,
        max_rows: u32,
    ) -> Result<table_stream_dispatch::Resultset, Duckerror> {
        let mut guard = cursors().lock().unwrap();
        let cur = guard.get_mut(&cursor).ok_or_else(|| {
            Duckerror::Invalidstate(format!("synteny-renderer: unknown cursor {cursor}"))
        })?;

        // Zero-row result on empty input, or on any subsequent poll after
        // the single payload row has been emitted -- both signal EOF.
        if max_rows == 0 || cur.already_emitted {
            return Ok(Vec::new());
        }
        let Some(payload) = cur.payload.take() else {
            cur.already_emitted = true;
            return Ok(Vec::new());
        };
        cur.already_emitted = true;

        let bytes_len = payload.len() as u32;
        let row = vec![Duckvalue::Blob(payload), Duckvalue::Uint32(bytes_len)];
        Ok(vec![project(row, &cur.projection)])
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
    if handle != HANDLE_RENDER {
        return Err(Duckerror::Invalidstate(format!(
            "synteny-renderer: unknown callback handle {handle}"
        )));
    }

    let tracks = parse_tracks(args.get(0))?;
    let features = parse_features(args.get(1))?;
    let links = parse_links(args.get(2))?;

    // EmptyInput -> zero-row cursor; other errors -> Invalidargument. That
    // way callers can `SELECT * FROM render_synteny_svg(NULL, NULL, NULL)`
    // and get an empty resultset rather than a SQL exception.
    let (payload, already_emitted) = match render::render_svg(&tracks, &features, &links) {
        Ok(bytes) => (Some(bytes), false),
        Err(RenderError::EmptyInput) => (None, true),
        Err(RenderError::InvalidModel(msg)) => {
            return Err(Duckerror::Invalidargument(format!(
                "synteny-renderer: {msg}"
            )));
        }
    };

    let id = next_cursor_id();
    let projected_cols = svg_columns_projected(&projection);
    cursors().lock().unwrap().insert(
        id,
        Cursor {
            payload,
            already_emitted,
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

// ---- arg parsing ------------------------------------------------------

#[derive(Deserialize)]
struct RawTrack {
    track_id: String,
    label: String,
    length: u32,
}

#[derive(Deserialize)]
struct RawFeature {
    track_id: String,
    feature_id: String,
    start_position: u32,
    end_position: u32,
    strand: i8,
    #[serde(default)]
    colour: Option<String>,
    #[serde(default)]
    label: Option<String>,
}

#[derive(Deserialize)]
struct RawLink {
    query_track: String,
    query_feature: String,
    subject_track: String,
    subject_feature: String,
    identity: f64,
    #[serde(default)]
    colour: Option<String>,
}

/// Common shape for the three list-of-struct args: NULL / missing / empty
/// string all yield an empty vec (so the renderer's own EmptyInput handling
/// kicks in); a `Duckvalue::Text(json)` -- the DuckLink 4.0.0 path -- or a
/// `Duckvalue::Complex(json)` -- kept as a fallback for the day DuckLink
/// preserves type-expressions -- both hand back the JSON string to parse.
fn expect_list_json<'a>(
    v: Option<&'a Duckvalue>,
    name: &str,
) -> Result<Option<&'a str>, Duckerror> {
    match v {
        None | Some(Duckvalue::Null) => Ok(None),
        Some(Duckvalue::Text(s)) if s.is_empty() => Ok(None),
        Some(Duckvalue::Text(s)) => Ok(Some(s.as_str())),
        // Kept in case DuckLink starts preserving complex type-expressions.
        Some(Duckvalue::Complex(Complexvalue { json, .. })) => Ok(Some(json.as_str())),
        _ => Err(Duckerror::Invalidargument(format!(
            "synteny-renderer: '{name}' must be a JSON VARCHAR (or NULL)"
        ))),
    }
}

fn parse_tracks(v: Option<&Duckvalue>) -> Result<Vec<Track>, Duckerror> {
    let Some(json) = expect_list_json(v, "tracks")? else {
        return Ok(Vec::new());
    };
    let raws: Vec<RawTrack> = serde_json::from_str(json).map_err(|e| {
        Duckerror::Invalidargument(format!(
            "synteny-renderer: cannot parse 'tracks' JSON: {e}"
        ))
    })?;
    Ok(raws
        .into_iter()
        .map(|r| Track {
            track_id: r.track_id,
            label: r.label,
            length: r.length,
        })
        .collect())
}

fn parse_features(v: Option<&Duckvalue>) -> Result<Vec<Feature>, Duckerror> {
    let Some(json) = expect_list_json(v, "features")? else {
        return Ok(Vec::new());
    };
    let raws: Vec<RawFeature> = serde_json::from_str(json).map_err(|e| {
        Duckerror::Invalidargument(format!(
            "synteny-renderer: cannot parse 'features' JSON: {e}"
        ))
    })?;
    Ok(raws
        .into_iter()
        .map(|r| Feature {
            track_id: r.track_id,
            feature_id: r.feature_id,
            start_position: r.start_position,
            end_position: r.end_position,
            strand: r.strand,
            colour: r.colour,
            label: r.label,
        })
        .collect())
}

fn parse_links(v: Option<&Duckvalue>) -> Result<Vec<Link>, Duckerror> {
    let Some(json) = expect_list_json(v, "links")? else {
        return Ok(Vec::new());
    };
    let raws: Vec<RawLink> = serde_json::from_str(json).map_err(|e| {
        Duckerror::Invalidargument(format!(
            "synteny-renderer: cannot parse 'links' JSON: {e}"
        ))
    })?;
    Ok(raws
        .into_iter()
        .map(|r| Link {
            query_track: r.query_track,
            query_feature: r.query_feature,
            subject_track: r.subject_track,
            subject_feature: r.subject_feature,
            identity: r.identity,
            colour: r.colour,
        })
        .collect())
}
