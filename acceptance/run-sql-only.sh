#!/usr/bin/env bash
#
# Pure-SQL variant of the acceptance test: NO Python driver, NO inlined-file
# SQL literals. The whole pipeline (read → parse → BLAST → conservation →
# summary → assertions) runs in one DuckLink invocation via
# acceptance/pipeline-sql-only.sql.
#
# The trick that makes this possible is `SET VARIABLE x = (SELECT ...);` +
# `getvariable('x')`. DuckDB accepts getvariable() as a TVF-argument
# "literal", so we can feed the output of `read_text(...)` into
# `genbank_scan(...)` and the output of an in-SQL `string_agg + json_quote`
# JSON builder into `blastn(...)` — the two TVF-arg subquery moves that
# would otherwise force us to build inputs out-of-band.
#
# Prereqs (see the repo README for build recipes):
#   ../ducklink/target/release/ducklink                            (native host)
#   ../ducklink/target/wasm32-wasip2/release/ducklink_{core,cli,loader}.wasm
#   ../ducklink/artifacts/extensions/jsonfns.wasm                  (for json_quote)
#   components/{blast, genome-format}/target/wasm32-wasip2/release/*.wasm
#
# Usage:
#   bash acceptance/run-sql-only.sh
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

BUILD_DIR="acceptance/build"
DATA_DIR="acceptance/data"
DUCKLINK_BIN="${DUCKLINK_BIN:-/Users/zacharywhitley/git/ducklink/target/release/ducklink}"
JSONFNS_WASM="${JSONFNS_WASM:-/Users/zacharywhitley/git/ducklink/artifacts/extensions/jsonfns.wasm}"
BLAST_WASM="components/blast/target/wasm32-wasip2/release/blast.wasm"
GENOME_WASM="components/genome-format/target/wasm32-wasip2/release/genome_format.wasm"
PIPELINE_SQL="acceptance/pipeline-sql-only.sql"

log() { printf '[run-sql-only.sh] %s\n' "$*" >&2; }

mkdir -p "${BUILD_DIR}"

# ---------------------------------------------------------------------------
# 1) Build both wasm components if missing. No Python parsing this time —
#    the SQL script does everything.
# ---------------------------------------------------------------------------
for pair in "blast:${BLAST_WASM}" "genome-format:${GENOME_WASM}"; do
    name="${pair%%:*}"
    wasm="${pair##*:}"
    if [ ! -f "${wasm}" ]; then
        log "building ${name} component"
        ( cd "components/${name}" && cargo component build --target wasm32-wasip2 --release )
    fi
    log "component: ${wasm} ($(wc -c < "${wasm}") bytes)"
done

if [ ! -f "${JSONFNS_WASM}" ]; then
    log "ERROR: ${JSONFNS_WASM} not found; sql/blast_macros.sql-style helpers need it."
    log "       Set JSONFNS_WASM=<path> if it lives elsewhere."
    exit 1
fi

if [ ! -x "${DUCKLINK_BIN}" ]; then
    log "ERROR: ${DUCKLINK_BIN} not found; build with:"
    log "         cd ../ducklink && cargo build --release -p ducklink-host --bin ducklink"
    log "       Set DUCKLINK_BIN=<path> to point elsewhere."
    exit 1
fi

# ---------------------------------------------------------------------------
# 2) Stage the extensions in a scratch dir the ducklink binary can scan.
# ---------------------------------------------------------------------------
EXT_DIR="${BUILD_DIR}/ext-sql-only"
mkdir -p "${EXT_DIR}"
cp "${BLAST_WASM}"  "${EXT_DIR}/blast.wasm"
cp "${GENOME_WASM}" "${EXT_DIR}/genome_format.wasm"
cp "${JSONFNS_WASM}" "${EXT_DIR}/jsonfns.wasm"
log "extensions staged in ${EXT_DIR}/"

# ---------------------------------------------------------------------------
# 3) Verify the GenBank inputs exist -- the SQL references them by name.
# ---------------------------------------------------------------------------
for gb in NC_001416.gb NC_001604.gb NC_002371.gb; do
    if [ ! -f "${DATA_DIR}/${gb}" ]; then
        log "ERROR: missing ${DATA_DIR}/${gb} -- see acceptance/README.md for the fetch."
        exit 1
    fi
done

# ---------------------------------------------------------------------------
# 4) Run the whole thing in one DuckLink invocation. --dir lets the CLI
#    read the .gb files via `read_text(...)`; the host preopens the repo
#    root under the same path so the SQL can use relative paths.
# ---------------------------------------------------------------------------
SQL_LOG="${BUILD_DIR}/pipeline-sql-only.out"
log "running DuckLink (SQL-only pipeline)"
"${DUCKLINK_BIN}" \
    --extensions-dir "${EXT_DIR}" \
    --dir "${REPO_ROOT}::${REPO_ROOT}" \
    -- "${REPO_ROOT}/${BUILD_DIR}/mem.duckdb" < "${PIPELINE_SQL}" \
    2>&1 \
    | grep -vE '^\[wasi-fs\]|^\[autoload\]|^\[resolver\]|^\[extension|^\[duckdb|^\[table-stream|^\[prefix\]' \
    | tee "${SQL_LOG}"

# clean the on-disk db so subsequent runs are fresh
rm -f "${REPO_ROOT}/${BUILD_DIR}/mem.duckdb"

# ---------------------------------------------------------------------------
# 5) Verify assertions (same shape as run-ducklink.sh).
# ---------------------------------------------------------------------------
log "checking assertion markers"
fail=0
for check in ASSERT_HITS_NONZERO ASSERT_ANY_CONSERVED ASSERT_PCT_IN_RANGE; do
    line=$(grep -E "\| ${check} " "${SQL_LOG}" | tail -1 || true)
    if [ -z "${line}" ]; then
        log "  ${check}: MISSING from pipeline output"
        fail=1
        continue
    fi
    status=$(printf '%s\n' "${line}" | awk -F'|' '{ gsub(/ /, "", $3); print $3 }')
    observed=$(printf '%s\n' "${line}" | awk -F'|' '{ gsub(/ /, "", $4); print $4 }')
    if [ "${status}" = "OK" ]; then
        log "  ${check}: OK (observed=${observed})"
    else
        log "  ${check}: FAIL (observed=${observed})"
        fail=1
    fi
done

if [ "${fail}" -ne 0 ]; then
    log "ACCEPTANCE-SQL-ONLY FAILED"
    exit 1
fi
log "ACCEPTANCE-SQL-ONLY PASSED"
