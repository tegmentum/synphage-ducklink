#!/usr/bin/env bash
#
# Native-DuckDB variant of the acceptance test: uses `duckdb` (the native
# binary from https://duckdb.org/) plus the `ducklink.duckdb_extension`
# native community extension, NOT the DuckLink wasm CLI. This is the shape
# an ordinary DuckDB user will run once ducklink lands in
# duckdb/community-extensions (`LOAD ducklink FROM community; ...`).
#
# Same pipeline as run-sql-only.sh (parse → BLAST → conservation → assert)
# and same 40-hits/8-conserved output, but zero WebAssembly on the CLI
# side — only the loaded extensions run wasm (blast.wasm, genome_format.wasm).
#
# Prereqs:
#   duckdb                                                 (native CLI on $PATH)
#   ../ducklink-extension/build/release/ducklink.duckdb_extension
#   components/{blast,genome-format}/target/wasm32-wasip2/release/*.wasm
#
# Override via env vars:
#   DUCKDB=/path/to/duckdb
#   DUCKLINK_EXT=/path/to/ducklink.duckdb_extension
#
# Usage:
#   bash acceptance/run-native-duckdb.sh
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

BUILD_DIR="acceptance/build"
DATA_DIR="${REPO_ROOT}/acceptance/data"
COMPONENT_BLAST="${REPO_ROOT}/components/blast/target/wasm32-wasip2/release/blast.wasm"
COMPONENT_GENOME="${REPO_ROOT}/components/genome-format/target/wasm32-wasip2/release/genome_format.wasm"
DUCKDB_BIN="${DUCKDB:-duckdb}"
DUCKLINK_EXT="${DUCKLINK_EXT:-/Users/zacharywhitley/git/ducklink-extension/build/release/ducklink.duckdb_extension}"
PIPELINE_SQL="${REPO_ROOT}/acceptance/pipeline-native.sql"

log() { printf '[run-native-duckdb.sh] %s\n' "$*" >&2; }

mkdir -p "${BUILD_DIR}"

# ---------------------------------------------------------------------------
# 1) Prereq checks — fail loudly with actionable messages.
# ---------------------------------------------------------------------------
if ! command -v "${DUCKDB_BIN}" >/dev/null 2>&1; then
    log "ERROR: '${DUCKDB_BIN}' not on PATH. Install DuckDB from https://duckdb.org/"
    log "       or set DUCKDB=<path>."
    exit 1
fi

if [ ! -f "${DUCKLINK_EXT}" ]; then
    log "ERROR: ducklink extension not found at ${DUCKLINK_EXT}"
    log "       Build with: (cd ../ducklink-extension && cargo build --release)"
    log "       Or set DUCKLINK_EXT=<path/to/ducklink.duckdb_extension>"
    exit 1
fi

for pair in "blast:${COMPONENT_BLAST}" "genome-format:${COMPONENT_GENOME}"; do
    name="${pair%%:*}"
    wasm="${pair##*:}"
    if [ ! -f "${wasm}" ]; then
        log "building ${name} wasm component"
        ( cd "components/${name}" && cargo component build --target wasm32-wasip2 --release )
    fi
    log "component: ${wasm} ($(wc -c < "${wasm}") bytes)"
done

for gb in NC_001416.gb NC_001604.gb NC_002371.gb; do
    if [ ! -f "${DATA_DIR}/${gb}" ]; then
        log "ERROR: missing ${DATA_DIR}/${gb} -- see acceptance/README.md for the fetch."
        exit 1
    fi
done

# ---------------------------------------------------------------------------
# 2) Run the pipeline.
#
# DUCKLINK_COMPONENTS names each wasm component up front (this is DuckDB's
# extension-load lifecycle: catalog registration happens at LOAD time).
#
# We pass data_dir as a SQL variable via `duckdb -cmd` so pipeline-native.sql
# stays free of absolute paths. -unsigned is required because both the
# ducklink extension and our components are locally-built (not signed by
# duckdb/community-extensions CI).
# ---------------------------------------------------------------------------
SQL_LOG="${BUILD_DIR}/pipeline-native.out"
log "running native DuckDB + ducklink extension"
DUCKLINK_COMPONENTS="blast=${COMPONENT_BLAST}:genome_format=${COMPONENT_GENOME}" \
"${DUCKDB_BIN}" -unsigned \
    -cmd "SET VARIABLE data_dir = '${DATA_DIR}';" \
    -cmd "LOAD '${DUCKLINK_EXT}';" \
    < "${PIPELINE_SQL}" \
    2>&1 | tee "${SQL_LOG}"

# ---------------------------------------------------------------------------
# 3) Verify the assertion rows the pipeline emits.
# ---------------------------------------------------------------------------
log "checking assertion markers"
fail=0
for check in ASSERT_HITS_NONZERO ASSERT_ANY_CONSERVED ASSERT_PCT_IN_RANGE; do
    line=$(grep -E "^\│ ${check} " "${SQL_LOG}" | tail -1 || true)
    if [ -z "${line}" ]; then
        # Also accept the ASCII table shape in case duckdb was invoked without
        # box drawing (e.g. .mode csv).
        line=$(grep -E "^\| ${check} " "${SQL_LOG}" | tail -1 || true)
    fi
    if [ -z "${line}" ]; then
        log "  ${check}: MISSING from pipeline output"
        fail=1
        continue
    fi
    status=$(printf '%s\n' "${line}" | awk -F'│|\\|' '{ gsub(/ /, "", $3); print $3 }')
    observed=$(printf '%s\n' "${line}" | awk -F'│|\\|' '{ gsub(/ /, "", $4); print $4 }')
    if [ "${status}" = "OK" ]; then
        log "  ${check}: OK (observed=${observed})"
    else
        log "  ${check}: FAIL (observed=${observed})"
        fail=1
    fi
done

if [ "${fail}" -ne 0 ]; then
    log "ACCEPTANCE-NATIVE-DUCKDB FAILED"
    exit 1
fi
log "ACCEPTANCE-NATIVE-DUCKDB PASSED"
