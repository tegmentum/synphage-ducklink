//! GenBank parser for DuckDB.
//!
//! THIN, GENERATED ducklink shim: a `wit_bindgen::generate!` block plus
//! one `datalink_extcore::duckdb_shim!` invocation. All logic + the
//! capability surface (including the `.gb`/`.gbk` replacement-scan
//! attribute on `genbank_read_path`) live ONCE in `genome-format-core`
//! (datalink); the registration ABI, handle table, dispatch arms, and
//! `Duckvalue` marshalling — plus the `files::register_replacement_scan`
//! call — are derived from the core's declaration.

wit_bindgen::generate!({
    path: "../../wit",
    world: "tegmentum:bio/genome-format-ducklink",
    generate_all,
});

datalink_extcore::duckdb_shim! {
    core = genome_format_core::Core;
    types = duckdb::extension::types;
    column_types = duckdb::extension::column_types;
    runtime = duckdb::extension::runtime;
    callback_dispatch = exports::duckdb::extension::callback_dispatch;
    guest = exports::duckdb::extension::guest;
    export = export;
    files = duckdb::extension::files;
}
