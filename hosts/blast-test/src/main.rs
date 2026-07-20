//! Standalone test host for `blast.wasm`.
//!
//! Loads the built component with wasmtime, mocks every DuckLink import with
//! trapping stubs (except `table-stream::register-filterable-table`, which
//! records the (name, callback-handle) pairs the guest registers), invokes
//! `guest::load()` so the two table functions register, then drives one
//! blastn scan end-to-end through the streaming table-dispatch cursor.
//!
//! This exists so we can verify the DuckLink dispatch surface works BEFORE
//! ducklink_{cli,core,loader}.wasm is built in the sibling `../ducklink/`
//! checkout — the acceptance-test scaffolding for the compound `blast`
//! world without the full runtime dependency.

use std::fs;
use std::io::Write;
use std::path::PathBuf;

use anyhow::{anyhow, Result};
use serde::Deserialize;

// wasmtime 46's Error is a NOT-a-`std::error::Error` custom type, so
// `anyhow::Context` can't be applied to `Result<T, wasmtime::Error>`
// directly. This macro lifts a wasmtime call into an `anyhow::Error` with
// a wrap message.
macro_rules! wasm_ctx {
    ($what:expr) => {
        |e: wasmtime::Error| ::anyhow::anyhow!("{}: {}", $what, e)
    };
}
use wasmtime::component::{Component, Linker, Resource, ResourceTable};
use wasmtime::{Config, Engine, Store};

mod bindings {
    // The compound `blast` world is defined in the tegmentum:bio package
    // (wit/world.wit). We point bindgen at `acceptance/wit-view/` — a
    // curated copy of the WIT files this binary needs (world.wit +
    // sequence-search.wit + genome-format.wit + synteny-renderer.wit +
    // deps/duckdb-extension/**) — to isolate ourselves from parallel-track
    // wit/ churn (new sibling *-world.wit files can name-collide with
    // existing interfaces before they land in the compound `blast` world).
    // Same package, same shape; just fewer files.
    wasmtime::component::bindgen!({
        path: "../../acceptance/wit-view",
        world: "blast",
    });
}

use bindings::duckdb::extension as ext;
use bindings::Blast;
use ext::types::{Complexvalue, Duckvalue};

// ------------------------------------------------------------------------
// Host state
// ------------------------------------------------------------------------

/// One entry per `table_stream::register_filterable_table` call. Captured
/// during `guest::load()` so the driver can look the blastn / blastp
/// callback-handle up by name afterwards.
#[derive(Debug, Clone)]
struct RegisteredTable {
    name: String,
    handle: u32,
}

struct State {
    table: ResourceTable,
    wasi: wasmtime_wasi::WasiCtx,
    /// Filled in by `register_filterable_table`.
    registrations: Vec<RegisteredTable>,
    /// Monotonic id we hand back from `register_filterable_table`; the
    /// component doesn't use it (see ducklink.rs), it's just a valid u32.
    next_reg_id: u32,
}

impl Default for State {
    fn default() -> Self {
        Self {
            table: ResourceTable::new(),
            wasi: wasmtime_wasi::WasiCtxBuilder::new()
                .inherit_stdio()
                .inherit_env()
                .build(),
            registrations: Vec::new(),
            next_reg_id: 1,
        }
    }
}

impl wasmtime_wasi::WasiView for State {
    fn ctx(&mut self) -> wasmtime_wasi::WasiCtxView<'_> {
        wasmtime_wasi::WasiCtxView {
            ctx: &mut self.wasi,
            table: &mut self.table,
        }
    }
}

// ------------------------------------------------------------------------
// Import trait impls
//
// The component's actual required imports (per `wasm-tools component wit
// blast.wasm`) are only `duckdb:extension/{types, table-stream, column-types}`
// plus wasi:{io,cli,random}. We still have to implement every import trait
// the compound `blast` world declares — wasmtime add_to_linker requires it —
// but the ones the guest never calls just trap.
// ------------------------------------------------------------------------

macro_rules! trap_import {
    ($what:literal) => {
        panic!(
            "unexpected host import call: {} — the test host expected the \
             component to only call `table_stream::register_filterable_table` \
             during load and the `table-stream-dispatch` cursor calls afterwards",
            $what
        )
    };
}

// types is a types-only interface (no free functions) — the Host trait is
// still generated as a marker.
impl ext::types::Host for State {}

// column-types is also a types-only interface but wasmtime bindgen still
// requires a marker impl for it (callback-dispatch references its `colvec`).
impl ext::column_types::Host for State {}

// logging
impl ext::logging::Host for State {
    fn log(&mut self, _level: ext::types::Loglevel, _message: String, _target: Option<String>) {
        trap_import!("logging::log");
    }
    fn log_fields(
        &mut self,
        _level: ext::types::Loglevel,
        _message: String,
        _fields: Vec<ext::types::Logfield>,
    ) {
        trap_import!("logging::log_fields");
    }
}

// config
impl ext::config::Host for State {
    fn provider_version(&mut self) -> String {
        trap_import!("config::provider_version");
    }
    fn list_keys(&mut self, _prefix: Option<String>) -> Vec<String> {
        trap_import!("config::list_keys");
    }
    fn get_string(&mut self, _path: String) -> Result<Option<String>, ext::types::Configerror> {
        trap_import!("config::get_string");
    }
    fn get_bool(&mut self, _path: String) -> Result<Option<bool>, ext::types::Configerror> {
        trap_import!("config::get_bool");
    }
    fn get_i64(&mut self, _path: String) -> Result<Option<i64>, ext::types::Configerror> {
        trap_import!("config::get_i64");
    }
    fn get_u64(&mut self, _path: String) -> Result<Option<u64>, ext::types::Configerror> {
        trap_import!("config::get_u64");
    }
    fn get_f64(&mut self, _path: String) -> Result<Option<f64>, ext::types::Configerror> {
        trap_import!("config::get_f64");
    }
    fn get_bytes(&mut self, _path: String) -> Result<Option<Vec<u8>>, ext::types::Configerror> {
        trap_import!("config::get_bytes");
    }
    fn get_string_list(
        &mut self,
        _path: String,
    ) -> Result<Option<Vec<String>>, ext::types::Configerror> {
        trap_import!("config::get_string_list");
    }
}

// files
impl ext::files::Host for State {
    fn register_replacement_scan(
        &mut self,
        _scan: ext::files::ReplacementScan,
    ) -> Result<ext::files::ReplacementScanId, String> {
        trap_import!("files::register_replacement_scan");
    }
    fn register_copy_handler(
        &mut self,
        _handler: ext::files::CopyHandler,
    ) -> Result<ext::files::CopyHandlerId, String> {
        trap_import!("files::register_copy_handler");
    }
}

// table-stream — THE ONLY IMPORT THE COMPONENT ACTUALLY CALLS.
impl ext::table_stream::Host for State {
    fn register_filterable_table(
        &mut self,
        name: String,
        _arguments: Vec<ext::types::Funcarg>,
        columns: Vec<ext::types::Columndef>,
        callback_handle: u32,
    ) -> Result<u32, ext::types::Duckerror> {
        let id = self.next_reg_id;
        self.next_reg_id += 1;
        eprintln!(
            "[host] register_filterable_table name={name:?} handle={callback_handle} \
             cols={} -> reg_id={id}",
            columns.len()
        );
        // We keep only (name, handle) in state — the driver looks the
        // handle up by name.
        self.registrations.push(RegisteredTable {
            name,
            handle: callback_handle,
        });
        Ok(id)
    }
}

// catalog
impl ext::catalog::Host for State {
    fn register_logical_type(
        &mut self,
        _ty: ext::catalog::LogicalType,
    ) -> Result<ext::catalog::LogicalTypeHandle, String> {
        trap_import!("catalog::register_logical_type");
    }
    fn register_cast(
        &mut self,
        _spec: ext::catalog::CastSpec,
        _callback: Resource<ext::runtime::CastCallback>,
    ) -> Result<(), String> {
        trap_import!("catalog::register_cast");
    }
    fn register_macro(&mut self, _def: ext::catalog::MacroDef) -> Result<(), String> {
        trap_import!("catalog::register_macro");
    }
}

// runtime: free fns + a resource per callback / registry family.
impl ext::runtime::Host for State {
    fn get_capability(
        &mut self,
        _kind: ext::types::Capabilitykind,
    ) -> Option<ext::runtime::Capability> {
        trap_import!("runtime::get_capability");
    }
    fn list_capabilities(&mut self) -> Vec<ext::types::Capabilitykind> {
        trap_import!("runtime::list_capabilities");
    }
}

impl ext::runtime::HostScalarCallback for State {
    fn new(&mut self, _handle: u32) -> Resource<ext::runtime::ScalarCallback> {
        trap_import!("runtime::scalar_callback::new");
    }
    fn call(
        &mut self,
        _self_: Resource<ext::runtime::ScalarCallback>,
        _args: Vec<Duckvalue>,
        _ctx: ext::types::Invokeinfo,
    ) -> Result<Duckvalue, ext::types::Duckerror> {
        trap_import!("runtime::scalar_callback::call");
    }
    fn drop(&mut self, _rep: Resource<ext::runtime::ScalarCallback>) -> wasmtime::Result<()> {
        trap_import!("runtime::scalar_callback::drop");
    }
}

impl ext::runtime::HostTableCallback for State {
    fn new(&mut self, _handle: u32) -> Resource<ext::runtime::TableCallback> {
        trap_import!("runtime::table_callback::new");
    }
    fn call(
        &mut self,
        _self_: Resource<ext::runtime::TableCallback>,
        _args: Vec<Duckvalue>,
    ) -> Result<ext::types::Resultset, ext::types::Duckerror> {
        trap_import!("runtime::table_callback::call");
    }
    fn drop(&mut self, _rep: Resource<ext::runtime::TableCallback>) -> wasmtime::Result<()> {
        trap_import!("runtime::table_callback::drop");
    }
}

impl ext::runtime::HostAggregateCallback for State {
    fn new(&mut self, _handle: u32) -> Resource<ext::runtime::AggregateCallback> {
        trap_import!("runtime::aggregate_callback::new");
    }
    fn call(
        &mut self,
        _self_: Resource<ext::runtime::AggregateCallback>,
        _rows: ext::types::Rowbatch,
    ) -> Result<Duckvalue, ext::types::Duckerror> {
        trap_import!("runtime::aggregate_callback::call");
    }
    fn drop(&mut self, _rep: Resource<ext::runtime::AggregateCallback>) -> wasmtime::Result<()> {
        trap_import!("runtime::aggregate_callback::drop");
    }
}

impl ext::runtime::HostPragmaCallback for State {
    fn new(&mut self, _handle: u32) -> Resource<ext::runtime::PragmaCallback> {
        trap_import!("runtime::pragma_callback::new");
    }
    fn call(
        &mut self,
        _self_: Resource<ext::runtime::PragmaCallback>,
        _args: Vec<Duckvalue>,
    ) -> Result<Option<Duckvalue>, ext::types::Duckerror> {
        trap_import!("runtime::pragma_callback::call");
    }
    fn drop(&mut self, _rep: Resource<ext::runtime::PragmaCallback>) -> wasmtime::Result<()> {
        trap_import!("runtime::pragma_callback::drop");
    }
}

impl ext::runtime::HostCastCallback for State {
    fn new(&mut self, _handle: u32) -> Resource<ext::runtime::CastCallback> {
        trap_import!("runtime::cast_callback::new");
    }
    fn call(
        &mut self,
        _self_: Resource<ext::runtime::CastCallback>,
        _value: Duckvalue,
    ) -> Result<Duckvalue, ext::types::Duckerror> {
        trap_import!("runtime::cast_callback::call");
    }
    fn drop(&mut self, _rep: Resource<ext::runtime::CastCallback>) -> wasmtime::Result<()> {
        trap_import!("runtime::cast_callback::drop");
    }
}

impl ext::runtime::HostScalarRegistry for State {
    fn register(
        &mut self,
        _self_: Resource<ext::runtime::ScalarRegistry>,
        _name: String,
        _arguments: Vec<ext::types::Funcarg>,
        _returns: ext::types::Logicaltype,
        _callback: Resource<ext::runtime::ScalarCallback>,
        _options: Option<ext::types::Funcopts>,
    ) -> Result<u32, ext::types::Duckerror> {
        trap_import!("runtime::scalar_registry::register");
    }
    fn drop(&mut self, _rep: Resource<ext::runtime::ScalarRegistry>) -> wasmtime::Result<()> {
        trap_import!("runtime::scalar_registry::drop");
    }
}

impl ext::runtime::HostTableRegistry for State {
    fn register(
        &mut self,
        _self_: Resource<ext::runtime::TableRegistry>,
        _name: String,
        _arguments: Vec<ext::types::Funcarg>,
        _columns: Vec<ext::types::Columndef>,
        _callback: Resource<ext::runtime::TableCallback>,
        _options: Option<ext::types::Extopts>,
    ) -> Result<u32, ext::types::Duckerror> {
        trap_import!("runtime::table_registry::register");
    }
    fn drop(&mut self, _rep: Resource<ext::runtime::TableRegistry>) -> wasmtime::Result<()> {
        trap_import!("runtime::table_registry::drop");
    }
}

impl ext::runtime::HostAggregateRegistry for State {
    fn register(
        &mut self,
        _self_: Resource<ext::runtime::AggregateRegistry>,
        _name: String,
        _arguments: Vec<ext::types::Funcarg>,
        _returns: ext::types::Logicaltype,
        _callback: Resource<ext::runtime::AggregateCallback>,
        _options: Option<ext::types::Funcopts>,
    ) -> Result<u32, ext::types::Duckerror> {
        trap_import!("runtime::aggregate_registry::register");
    }
    fn drop(&mut self, _rep: Resource<ext::runtime::AggregateRegistry>) -> wasmtime::Result<()> {
        trap_import!("runtime::aggregate_registry::drop");
    }
}

impl ext::runtime::HostPragmaRegistry for State {
    fn register_call(
        &mut self,
        _self_: Resource<ext::runtime::PragmaRegistry>,
        _name: String,
        _arguments: Vec<ext::types::Funcarg>,
        _returns: ext::types::Logicaltype,
        _callback: Resource<ext::runtime::PragmaCallback>,
        _options: Option<ext::types::Extopts>,
    ) -> Result<u32, ext::types::Duckerror> {
        trap_import!("runtime::pragma_registry::register_call");
    }
    fn drop(&mut self, _rep: Resource<ext::runtime::PragmaRegistry>) -> wasmtime::Result<()> {
        trap_import!("runtime::pragma_registry::drop");
    }
}

impl ext::runtime::HostMacroRegistry for State {
    fn register_scalar(
        &mut self,
        _self_: Resource<ext::runtime::MacroRegistry>,
        _name: String,
        _parameters: Vec<String>,
        _body_sql: String,
        _options: Option<ext::types::Extopts>,
    ) -> Result<bool, ext::types::Duckerror> {
        trap_import!("runtime::macro_registry::register_scalar");
    }
    fn drop(&mut self, _rep: Resource<ext::runtime::MacroRegistry>) -> wasmtime::Result<()> {
        trap_import!("runtime::macro_registry::drop");
    }
}

// ------------------------------------------------------------------------
// Sample sequences.
//
// The FASTA files at ../../examples/tiny-blastn/{queries,subjects}.fasta
// carry the same data — they exist so a reader can inspect the inputs
// without opening this source. The constants below are the runtime source
// of truth so the binary works whether or not the fasta files ship
// alongside it.
// ------------------------------------------------------------------------

const QUERIES: &[(&str, &str)] = &[
    (
        "gene_A",
        "ATGCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG",
    ),
    (
        "gene_B",
        "GGGCCCTTTAAAGGGCCCTTTAAAGGGCCCTTTAAAGGGCCCTTTAAAGGGCCCTTTAAA",
    ),
    (
        "gene_C",
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
    ),
];

const SUBJECTS: &[(&str, &str)] = &[
    (
        "genome_1",
        "TTATGCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGAA",
    ),
    (
        "genome_2",
        "CCAGGGCCCTTTAAAGGGCCCTTTAAAGGGCCCTTTAAAGGGCCCTTTAAAGGGCCCTTTAAATT",
    ),
    (
        "genome_3",
        "GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA",
    ),
];

// The type-expression strings the component's registration declares.
// Matched verbatim so a real DuckLink host and this test host encode the
// same `complex(complexvalue)` payload.
const TYPE_EXPR_SEQ_LIST: &str = "LIST(STRUCT(key VARCHAR, data VARCHAR))";
const TYPE_EXPR_OPTIONS: &str =
    "STRUCT(evalue_max DOUBLE, max_target_seqs INTEGER, min_identity DOUBLE)";

fn sequence_list_arg(seqs: &[(&str, &str)]) -> Duckvalue {
    let json = serde_json::Value::Array(
        seqs.iter()
            .map(|(k, d)| {
                serde_json::json!({
                    "key": k,
                    "data": d,
                })
            })
            .collect(),
    );
    Duckvalue::Complex(Complexvalue {
        type_expr: TYPE_EXPR_SEQ_LIST.to_string(),
        json: json.to_string(),
    })
}

fn options_arg() -> Duckvalue {
    // NULL is enough — parse_options() defaults to no filters. Kept as a
    // reference example of the STRUCT type-expr for when we want to
    // exercise the real path.
    let _ = TYPE_EXPR_OPTIONS; // silence dead-code lint if never used
    Duckvalue::Null
}

// ------------------------------------------------------------------------
// Component loading + driver
// ------------------------------------------------------------------------

fn wasm_path() -> PathBuf {
    let manifest = env!("CARGO_MANIFEST_DIR");
    PathBuf::from(manifest)
        .join("../../components/blast/target/wasm32-wasip2/release/blast.wasm")
}

fn load_component() -> Result<(Store<State>, Blast)> {
    let mut config = Config::new();
    config.wasm_component_model(true);
    let engine = Engine::new(&config).map_err(wasm_ctx!("Engine::new"))?;

    let mut linker: Linker<State> = Linker::new(&engine);
    Blast::add_to_linker::<State, wasmtime::component::HasSelf<State>>(&mut linker, |s| s)
        .map_err(wasm_ctx!("Blast::add_to_linker"))?;
    wasmtime_wasi::p2::add_to_linker_sync(&mut linker)
        .map_err(wasm_ctx!("wasi add_to_linker_sync"))?;

    let mut store = Store::new(&engine, State::default());

    let path = wasm_path();
    let bytes = std::fs::read(&path)
        .map_err(|e| anyhow!("reading blast.wasm at {}: {}", path.display(), e))?;
    // Component preamble: `\0asm` + version 0x0d 0x00 + layer 0x01 0x00.
    if bytes.len() < 8 || &bytes[0..4] != b"\0asm" || bytes[6..8] != [0x01, 0x00] {
        return Err(anyhow!(
            "{} does not look like a wasm component (bad magic / layer)",
            path.display()
        ));
    }
    let component = Component::new(&engine, &bytes).map_err(wasm_ctx!("Component::new"))?;
    let instance =
        Blast::instantiate(&mut store, &component, &linker).map_err(wasm_ctx!("instantiate"))?;
    Ok((store, instance))
}

fn duck_to_display(v: &Duckvalue) -> String {
    match v {
        Duckvalue::Null => "NULL".to_string(),
        Duckvalue::Boolean(b) => b.to_string(),
        Duckvalue::Int64(x) => x.to_string(),
        Duckvalue::Uint64(x) => x.to_string(),
        Duckvalue::Float64(x) => format!("{x:.6}"),
        Duckvalue::Text(s) => s.clone(),
        Duckvalue::Blob(b) => format!("<{} bytes>", b.len()),
        Duckvalue::Int32(x) => x.to_string(),
        Duckvalue::Timestamp(x) => x.to_string(),
        Duckvalue::Int8(x) => x.to_string(),
        Duckvalue::Int16(x) => x.to_string(),
        Duckvalue::Uint8(x) => x.to_string(),
        Duckvalue::Uint16(x) => x.to_string(),
        Duckvalue::Uint32(x) => x.to_string(),
        Duckvalue::Float32(x) => format!("{x:.6}"),
        Duckvalue::Date(x) => x.to_string(),
        Duckvalue::Time(x) => x.to_string(),
        Duckvalue::Timestamptz(x) => x.to_string(),
        Duckvalue::Decimal(d) => format!("decimal(lo={},hi={},w={},s={})", d.lower, d.upper, d.width, d.scale),
        Duckvalue::Interval(i) => format!("interval({}m,{}d,{}us)", i.months, i.days, i.micros),
        Duckvalue::Uuid(u) => format!("uuid({:x}-{:x})", u.hi, u.lo),
        Duckvalue::Complex(c) => format!("complex({}: {})", c.type_expr, c.json),
    }
}

/// The 13 hit columns registered by the component, in the order they emit.
/// Used for the printout header — the component returns them in this order
/// when projection is empty.
const HIT_COLUMNS: &[&str] = &[
    "query_key",
    "subject_key",
    "query_start",
    "query_end",
    "subject_start",
    "subject_end",
    "strand",
    "identity_count",
    "alignment_length",
    "percent_identity",
    "bit_score",
    "raw_score",
    "evalue",
];

// ------------------------------------------------------------------------
// CLI: optional file-based inputs. When no --queries / --subjects flags are
// supplied the built-in QUERIES / SUBJECTS constants are used (smoke test);
// otherwise the payloads are loaded from JSON files with the shape
//   [ { "key": "...", "data": "..." }, ... ]
// and the driven scan is written out as a TSV to --hits-out (or stdout).
// ------------------------------------------------------------------------

#[derive(Deserialize)]
struct SeqRecord {
    key: String,
    data: String,
}

#[derive(Default)]
struct Cli {
    queries_path: Option<PathBuf>,
    subjects_path: Option<PathBuf>,
    hits_out: Option<PathBuf>,
    function: String, // "blastn" | "blastp"
    /// Verbatim JSON that becomes the `options` STRUCT argument. If None,
    /// the third arg is DuckLink NULL and the component falls back to
    /// unfiltered defaults.
    options_json: Option<String>,
}

fn parse_cli() -> Result<Cli> {
    let mut cli = Cli {
        function: "blastn".to_string(),
        ..Default::default()
    };
    let mut args = std::env::args().skip(1);
    while let Some(a) = args.next() {
        match a.as_str() {
            "--queries" => {
                cli.queries_path = Some(PathBuf::from(
                    args.next().ok_or_else(|| anyhow!("--queries needs a value"))?,
                ));
            }
            "--subjects" => {
                cli.subjects_path = Some(PathBuf::from(
                    args.next().ok_or_else(|| anyhow!("--subjects needs a value"))?,
                ));
            }
            "--hits-out" => {
                cli.hits_out = Some(PathBuf::from(
                    args.next().ok_or_else(|| anyhow!("--hits-out needs a value"))?,
                ));
            }
            "--function" => {
                cli.function = args
                    .next()
                    .ok_or_else(|| anyhow!("--function needs a value"))?;
                if cli.function != "blastn" && cli.function != "blastp" {
                    return Err(anyhow!(
                        "--function must be 'blastn' or 'blastp', got {:?}",
                        cli.function
                    ));
                }
            }
            "--options-json" => {
                cli.options_json = Some(
                    args.next()
                        .ok_or_else(|| anyhow!("--options-json needs a value"))?,
                );
            }
            "-h" | "--help" => {
                println!(
                    "blast-test [--queries PATH.json] [--subjects PATH.json] \
                     [--hits-out PATH.tsv] [--function blastn|blastp] \
                     [--options-json '{{\"evalue_max\":1e-5,...}}']\n\
                     \n\
                     With no flags: runs the built-in smoke test against tiny \
                     hardcoded sequences and prints hits to stdout.\n\
                     \n\
                     With --queries + --subjects: loads JSON arrays of \
                     {{key, data}} objects, drives one scan through the \
                     table-stream cursor, and writes hits TSV to --hits-out \
                     (or stdout if omitted).\n\
                     \n\
                     --options-json takes the exact JSON shape the component's \
                     `options` STRUCT argument accepts: any subset of \
                     evalue_max (double), max_target_seqs (int), \
                     min_identity (double).",
                );
                std::process::exit(0);
            }
            other => return Err(anyhow!("unknown flag: {other}")),
        }
    }
    Ok(cli)
}

fn read_seq_json(path: &PathBuf) -> Result<Vec<(String, String)>> {
    let bytes = fs::read(path).map_err(|e| anyhow!("reading {}: {}", path.display(), e))?;
    let raws: Vec<SeqRecord> = serde_json::from_slice(&bytes)
        .map_err(|e| anyhow!("parsing {} as JSON array of {{key,data}}: {}", path.display(), e))?;
    Ok(raws.into_iter().map(|r| (r.key, r.data)).collect())
}

fn sequence_list_arg_owned(seqs: &[(String, String)]) -> Duckvalue {
    let json = serde_json::Value::Array(
        seqs.iter()
            .map(|(k, d)| {
                serde_json::json!({
                    "key": k,
                    "data": d,
                })
            })
            .collect(),
    );
    Duckvalue::Complex(Complexvalue {
        type_expr: TYPE_EXPR_SEQ_LIST.to_string(),
        json: json.to_string(),
    })
}

fn duck_to_tsv_cell(v: &Duckvalue) -> String {
    // TSV-safe: tabs, newlines, and CR replaced. `duck_to_display`'s
    // format is already flat text but we defensively strip separators
    // so downstream `read_csv(..., sep='\t')` in DuckDB is unambiguous.
    let s = duck_to_display(v);
    s.replace('\t', " ").replace('\n', " ").replace('\r', " ")
}

fn main() -> Result<()> {
    let cli = parse_cli()?;
    let (mut store, instance) = load_component()?;

    // --- 1) load(): registers blastn + blastp ---
    let guest = instance.duckdb_extension_guest();
    let load_res = guest
        .call_load(&mut store)
        .map_err(wasm_ctx!("guest::load host call"))?
        .map_err(|e| anyhow!("guest::load returned Duckerror: {e:?}"))?;
    eprintln!(
        "[host] guest.load -> name={:?} version={:?} requires={:?}",
        load_res.name, load_res.version, load_res.requires
    );

    // --- 2) look up the requested handle we captured during load ---
    let want_name = cli.function.as_str();
    let handle = store
        .data()
        .registrations
        .iter()
        .find(|r| r.name == want_name)
        .map(|r| r.handle)
        .ok_or_else(|| anyhow!("component did not register a {:?} table function", want_name))?;
    let all_regs: Vec<String> = store
        .data()
        .registrations
        .iter()
        .map(|r| format!("{}#{}", r.name, r.handle))
        .collect();
    eprintln!(
        "[host] registered table functions: [{}]",
        all_regs.join(", ")
    );

    // --- 3) build args ---
    let td = instance.duckdb_extension_table_stream_dispatch();

    let options_dv = match cli.options_json.as_deref() {
        None => options_arg(),
        Some(raw) => {
            // Validate the shape early — a bad string will surface as a
            // Duckerror::Invalidargument out of the component, but catching
            // it here gives a clearer host-side message.
            let _: serde_json::Value = serde_json::from_str(raw).map_err(|e| {
                anyhow!("--options-json is not valid JSON: {e}")
            })?;
            Duckvalue::Complex(Complexvalue {
                type_expr: TYPE_EXPR_OPTIONS.to_string(),
                json: raw.to_string(),
            })
        }
    };

    let (queries_owned, subjects_owned);
    let args = match (&cli.queries_path, &cli.subjects_path) {
        (None, None) => {
            // built-in smoke-test path
            eprintln!("[host] using built-in QUERIES/SUBJECTS constants");
            vec![
                sequence_list_arg(QUERIES),
                sequence_list_arg(SUBJECTS),
                options_dv,
            ]
        }
        (Some(qp), Some(sp)) => {
            queries_owned = read_seq_json(qp)?;
            subjects_owned = read_seq_json(sp)?;
            eprintln!(
                "[host] loaded queries={} subjects={} from {} and {}",
                queries_owned.len(),
                subjects_owned.len(),
                qp.display(),
                sp.display()
            );
            vec![
                sequence_list_arg_owned(&queries_owned),
                sequence_list_arg_owned(&subjects_owned),
                options_dv,
            ]
        }
        _ => {
            return Err(anyhow!(
                "--queries and --subjects must be provided together"
            ));
        }
    };

    let open = td
        .call_call_table_open_filtered(&mut store, handle, &args, &[], &[])
        .map_err(wasm_ctx!("call_table_open_filtered host call"))?
        .map_err(|e| anyhow!("call_table_open_filtered Duckerror: {e:?}"))?;
    eprintln!(
        "[host] open cursor={} columns={:?}",
        open.cursor,
        open.columns.iter().map(|c| &c.name).collect::<Vec<_>>()
    );

    // --- 4) drain the cursor in max-rows=1024 batches ---
    let mut hits: Vec<Vec<Duckvalue>> = Vec::new();
    loop {
        let batch = td
            .call_call_table_next(&mut store, handle, open.cursor, 1024)
            .map_err(wasm_ctx!("call_table_next host call"))?
            .map_err(|e| anyhow!("call_table_next Duckerror: {e:?}"))?;
        if batch.is_empty() {
            break;
        }
        for row in batch {
            hits.push(row);
        }
    }

    // --- 5) close ---
    td.call_call_table_close(&mut store, handle, open.cursor)
        .map_err(wasm_ctx!("call_table_close host call"))?
        .map_err(|e| anyhow!("call_table_close Duckerror: {e:?}"))?;

    // --- 6) emit ---
    let function_upper = cli.function.to_uppercase();
    if let Some(out) = cli.hits_out.as_ref() {
        let mut f = fs::File::create(out)
            .map_err(|e| anyhow!("creating {}: {}", out.display(), e))?;
        writeln!(f, "{}", HIT_COLUMNS.join("\t"))?;
        for row in &hits {
            let cells: Vec<String> = row.iter().map(duck_to_tsv_cell).collect();
            writeln!(f, "{}", cells.join("\t"))?;
        }
        eprintln!(
            "[host] wrote {} hit row(s) to {} ({} scan complete)",
            hits.len(),
            out.display(),
            function_upper
        );
    } else {
        println!();
        println!("{} scan complete: {} hit(s)", function_upper, hits.len());
        println!();
        println!("{}", HIT_COLUMNS.join("\t"));
        println!("{}", "-".repeat(HIT_COLUMNS.len() * 8));
        for row in &hits {
            let cells: Vec<String> = row.iter().map(duck_to_display).collect();
            println!("{}", cells.join("\t"));
        }
        if hits.is_empty() && cli.queries_path.is_none() {
            eprintln!(
                "[host] WARNING: no hits produced. Sample data may need adjusting, \
                 or the aligner's default thresholds are stricter than expected."
            );
        }
    }

    Ok(())
}
