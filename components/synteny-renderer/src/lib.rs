//! Synteny-plot SVG renderer for DuckDB.
//!
//! THIN, GENERATED ducklink shim: a `wit_bindgen::generate!` block plus
//! one `datalink_extcore::duckdb_shim!` invocation. All logic + the
//! capability surface live ONCE in `synteny-renderer-core` (datalink);
//! the registration ABI, handle table, dispatch arm, and `Duckvalue`
//! marshalling are derived from the core's declaration.

wit_bindgen::generate!({
    path: "../../wit",
    world: "tegmentum:bio/synteny-renderer-ducklink",
    generate_all,
});

datalink_extcore::duckdb_shim! {
    core = synteny_renderer_core::Core;
    types = duckdb::extension::types;
    column_types = duckdb::extension::column_types;
    runtime = duckdb::extension::runtime;
    callback_dispatch = exports::duckdb::extension::callback_dispatch;
    guest = exports::duckdb::extension::guest;
    export = export;
}
