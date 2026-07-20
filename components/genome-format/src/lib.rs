//! bio-genome-format — GenBank flat-file parsing as a wasm component
//! exporting `tegmentum:bio/genome-format@0.1.0`.
//!
//! One entry point: `parse(list<u8>) -> result<parsed-genbank, parse-error>`.
//! The four `list<...>` fields on `parsed-genbank` map 1:1 to the four
//! DuckDB relations a future `genbank_scan()` table function will expose:
//!
//! - records     — one row per LOCUS block
//! - features    — one row per FEATURES table entry, with a stable
//!                 `feature-index` within its enclosing record
//! - qualifiers  — one row per /name=value pair (multi-value qualifiers
//!                 arrive as multiple rows — flattening is a SQL concern)
//! - sequences   — one row per record, holding the ORIGIN bytes
//!
//! Now dual-exported: alongside `tegmentum:bio/genome-format` this
//! component also implements DuckLink's dispatch surface
//! (`duckdb:extension/{guest, callback-dispatch, table-stream-dispatch}`)
//! and registers a `genbank_scan(paths VARCHAR)` table function during
//! `load()`. See `src/ducklink.rs` for the bridge and
//! `wit/genome-format-world.wit` for the compound world declaration.

#[allow(warnings)]
mod bindings;

mod ducklink;
mod location;
mod model;
mod parser;

use bindings::exports::tegmentum::bio::genome_format::{
    GenbankFeature, GenbankQualifier, GenbankRecord, GenbankSequence,
    Guest as GenomeFormatGuest, ParseError as WitParseError, ParsedGenbank,
};

pub(crate) struct Component;

impl GenomeFormatGuest for Component {
    fn parse(data: Vec<u8>) -> Result<ParsedGenbank, WitParseError> {
        parser::parse_genbank(&data)
            .map(into_wit_parsed)
            .map_err(into_wit_error)
    }
}

// ---- boundary marshalling -------------------------------------------------
//
// The generated types are structurally identical to `model::*` but live in a
// different module and can't be blanket-`Into`'d. Explicit converters keep
// the boundary honest: any field-shape drift shows up as a compile error
// here rather than as silently wrong data on the far side.

fn into_wit_parsed(p: model::Parsed) -> ParsedGenbank {
    ParsedGenbank {
        records: p.records.into_iter().map(into_wit_record).collect(),
        features: p.features.into_iter().map(into_wit_feature).collect(),
        qualifiers: p.qualifiers.into_iter().map(into_wit_qualifier).collect(),
        sequences: p.sequences.into_iter().map(into_wit_sequence).collect(),
    }
}

fn into_wit_record(r: model::Record) -> GenbankRecord {
    GenbankRecord {
        record_id: r.record_id,
        accession: r.accession,
        version: r.version,
        organism: r.organism,
        length: r.length,
    }
}

fn into_wit_feature(f: model::Feature) -> GenbankFeature {
    GenbankFeature {
        record_id: f.record_id,
        feature_index: f.feature_index,
        feature_type: f.feature_type,
        start_position: f.start_position,
        end_position: f.end_position,
        strand: f.strand,
    }
}

fn into_wit_qualifier(q: model::Qualifier) -> GenbankQualifier {
    GenbankQualifier {
        record_id: q.record_id,
        feature_index: q.feature_index,
        name: q.name,
        value: q.value,
    }
}

fn into_wit_sequence(s: model::Sequence) -> GenbankSequence {
    GenbankSequence {
        record_id: s.record_id,
        data: s.data,
    }
}

fn into_wit_error(e: model::ParseError) -> WitParseError {
    match e {
        model::ParseError::Malformed(s) => WitParseError::Malformed(s),
        model::ParseError::UnsupportedVersion(s) => WitParseError::UnsupportedVersion(s),
    }
}

bindings::export!(Component with_types_in bindings);
