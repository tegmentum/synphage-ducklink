#!/usr/bin/env bash
#
# End-to-end acceptance test for the synphage-ducklink pipeline against
# real NCBI phage genomes (lambda / T7 / P22).
#
# STATUS: RETIRED as of the blast v5.0.0 port.
#
# This driver uses hosts/blast-test/, a wasmtime-based mock host that
# reimplements just enough of DuckLink's *streaming* table dispatch
# (`table_stream::register_filterable_table` -> `table-stream-dispatch`
# cursor calls) to run one BLASTN scan without needing DuckLink to be
# built. It was scaffolding built when neither real-DuckLink driver
# below could run yet.
#
# blast has since been refactored to register through `runtime::TableRegistry`
# and dispatch through row-major `callback_dispatch::call_table` — the API
# ducklink-extension v5.0.0 actually surfaces into DuckDB. The mock host's
# `runtime::get_capability` import handler still traps, so this driver
# no longer runs end-to-end. It is left in the tree for reference but the
# three drivers below cover its purpose:
#
#   acceptance/run-native-duckdb.sh    <- native duckdb + ducklink extension
#   acceptance/run-ducklink.sh         <- DuckLink wasm CLI, Python-driven
#   acceptance/run-sql-only.sh         <- DuckLink wasm CLI, pure-SQL
#
# Restoring this driver means implementing runtime::TableRegistry +
# callback_dispatch::call_table dispatch in hosts/blast-test/src/main.rs.
# That's ~100 lines of Rust, worth it only if the mock host has to stand
# on its own again (e.g. when both DuckLink paths are unavailable).
#
# --- historical stages -----------------------------------------------------
#   1. Build the blast wasm component.
#   2. Build the wasmtime-based test host binary.
#   3. Parse each acceptance/data/*.gb file, extract CDS sequences.
#   4. Drive one BLASTN scan through the wasmtime host.
#   5. Load hits into DuckDB, apply the three sql/*.sql views.
#   6. Assert the same three markers the other drivers check.
set -euo pipefail

echo "[run.sh] RETIRED: mock-host driver no longer runs end-to-end after the" >&2
echo "         blast v5.0.0 port. Use:" >&2
echo "           bash acceptance/run-native-duckdb.sh   # daily-driver shape" >&2
echo "           bash acceptance/run-ducklink.sh        # DuckLink wasm CLI + Python" >&2
echo "           bash acceptance/run-sql-only.sh        # DuckLink wasm CLI, pure SQL" >&2
echo "         See the header of this file for context." >&2
exit 2

# ---------------------------------------------------------------------------
# Path anchoring: run from the repo root regardless of CWD when invoked.
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

BUILD_DIR="acceptance/build"
DATA_DIR="acceptance/data"
COMPONENT_WASM="components/blast/target/wasm32-wasip2/release/blast.wasm"
HOST_BIN="hosts/blast-test/target/release/blast-test"
HITS_TSV="${BUILD_DIR}/hits.tsv"
QUERIES_JSON="${BUILD_DIR}/queries.json"
SUBJECTS_JSON="${BUILD_DIR}/subjects.json"
FEATURES_TSV="${BUILD_DIR}/genome_features.tsv"
SQL_LOG="${BUILD_DIR}/pipeline.out"

EVALUE_MAX="${EVALUE_MAX:-1e-5}"
PER_GENOME="${PER_GENOME:-8}"
MAX_CDS_LEN="${MAX_CDS_LEN:-900}"

log() { printf '[run.sh] %s\n' "$*" >&2; }

mkdir -p "${BUILD_DIR}"

# ---------------------------------------------------------------------------
# 1) Build the blast wasm component if missing. (Skipped when the artifact
#    already exists so re-runs are fast; delete blast.wasm to force a
#    rebuild.)
# ---------------------------------------------------------------------------
if [ ! -f "${COMPONENT_WASM}" ]; then
    if ! command -v cargo-component >/dev/null 2>&1; then
        log "ERROR: ${COMPONENT_WASM} missing and cargo-component not on PATH"
        exit 1
    fi
    log "building blast wasm component"
    ( cd components/blast && cargo component build --target wasm32-wasip2 --release )
fi
log "component: ${COMPONENT_WASM} ($(wc -c < "${COMPONENT_WASM}") bytes)"

# ---------------------------------------------------------------------------
# 2) Build the wasmtime test host if missing. Same skip-if-present logic.
# ---------------------------------------------------------------------------
if [ ! -x "${HOST_BIN}" ]; then
    log "building blast-test wasmtime host binary"
    ( cd hosts/blast-test && cargo build --release )
fi
log "host: ${HOST_BIN}"

# ---------------------------------------------------------------------------
# 3) Parse the GenBank inputs.
# ---------------------------------------------------------------------------
if ! compgen -G "${DATA_DIR}/*.gb" >/dev/null && \
   ! compgen -G "${DATA_DIR}/*.gbk" >/dev/null; then
    log "ERROR: no *.gb/*.gbk files under ${DATA_DIR}"
    log "       See acceptance/README.md for how to fetch them from NCBI."
    exit 1
fi

log "parsing GenBank inputs"
python3 acceptance/gb_parse.py \
    --data-dir "${DATA_DIR}" \
    --out-dir "${BUILD_DIR}" \
    --per-genome "${PER_GENOME}" \
    --max-len "${MAX_CDS_LEN}"

# ---------------------------------------------------------------------------
# 4) Run BLASTN through the wasmtime host with the pushed-down evalue
#    threshold. NULL options would let every trivial match through, so the
#    filter is essential for a meaningful signal.
# ---------------------------------------------------------------------------
log "running BLASTN (evalue_max=${EVALUE_MAX})"
"${HOST_BIN}" \
    --queries   "${QUERIES_JSON}" \
    --subjects  "${SUBJECTS_JSON}" \
    --hits-out  "${HITS_TSV}" \
    --options-json "{\"evalue_max\":${EVALUE_MAX}}"
log "hits: ${HITS_TSV} ($(( $(wc -l < "${HITS_TSV}") - 1 )) rows)"

# ---------------------------------------------------------------------------
# 5) DuckDB pipeline. Runs sql/best_hits.sql -> sql/gene_conservation.sql
#    -> sql/summary.sql on top of the loaded hits and prints the results.
# ---------------------------------------------------------------------------
log "loading into DuckDB and applying conservation views"
if ! command -v duckdb >/dev/null 2>&1; then
    log "ERROR: duckdb not on PATH. Install via 'brew install duckdb' or"
    log "       download the release binary from https://duckdb.org/."
    exit 1
fi
duckdb -cmd "SET VARIABLE hits_path='${HITS_TSV}';" \
       -cmd "SET VARIABLE features_path='${FEATURES_TSV}';" \
       < acceptance/pipeline.sql | tee "${SQL_LOG}"

# ---------------------------------------------------------------------------
# 6) Verify the three assertion rows the SQL script emits. Each is a single
#    row of the form  "<name>  OK|FAIL  <count>".
# ---------------------------------------------------------------------------
log "checking assertion markers"
fail=0
for check in ASSERT_HITS_NONZERO ASSERT_ANY_CONSERVED ASSERT_PCT_IN_RANGE; do
    line=$(grep -E "^${check}\b" "${SQL_LOG}" | tail -1 || true)
    if [ -z "${line}" ]; then
        log "  ${check}: MISSING from pipeline output"
        fail=1
        continue
    fi
    status=$(printf '%s\n' "${line}" | awk '{print $2}')
    observed=$(printf '%s\n' "${line}" | awk '{print $3}')
    if [ "${status}" = "OK" ]; then
        log "  ${check}: OK (observed=${observed})"
    else
        log "  ${check}: FAIL (observed=${observed})"
        fail=1
    fi
done

if [ "${fail}" -ne 0 ]; then
    log "ACCEPTANCE FAILED"
    exit 1
fi
log "ACCEPTANCE PASSED"
