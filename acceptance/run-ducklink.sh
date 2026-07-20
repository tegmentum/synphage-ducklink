#!/usr/bin/env bash
#
# Real-DuckLink variant of acceptance/run.sh. Same pipeline (parse phage
# genbank → blastn → conservation SQL → summary + assertions), but the BLAST
# step goes through the actual `ducklink` binary and the `blast.wasm`
# extension instead of the wasmtime-mock host at hosts/blast-test.
#
# Locks in the "real DuckLink" end-to-end as a repeatable test.
#
# Prerequisites (see /Users/zacharywhitley/git/synphage-ducklink/README.md
# for the build recipe):
#   ../ducklink/target/release/ducklink              (native host binary)
#   ../ducklink/target/wasm32-wasip2/release/ducklink_{core,cli,loader}.wasm
#
# Usage:
#   bash acceptance/run-ducklink.sh
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

BUILD_DIR="acceptance/build"
DATA_DIR="acceptance/data"
COMPONENT_WASM="components/blast/target/wasm32-wasip2/release/blast.wasm"
DUCKLINK_BIN="${DUCKLINK_BIN:-/Users/zacharywhitley/git/ducklink/target/release/ducklink}"

EVALUE_MAX="${EVALUE_MAX:-1e-5}"
PER_GENOME="${PER_GENOME:-8}"
MAX_CDS_LEN="${MAX_CDS_LEN:-900}"

log() { printf '[run-ducklink.sh] %s\n' "$*" >&2; }

mkdir -p "${BUILD_DIR}"

# ---------------------------------------------------------------------------
# 1) Build the blast wasm component if missing.
# ---------------------------------------------------------------------------
if [ ! -f "${COMPONENT_WASM}" ]; then
    log "building blast wasm component"
    ( cd components/blast && cargo component build --target wasm32-wasip2 --release )
fi
log "component: ${COMPONENT_WASM} ($(wc -c < "${COMPONENT_WASM}") bytes)"

# ---------------------------------------------------------------------------
# 2) Stage the .wasm under a scratch extensions dir the ducklink binary
#    can scan (avoids polluting the user's default artifacts/extensions).
# ---------------------------------------------------------------------------
EXT_DIR="${BUILD_DIR}/ext"
mkdir -p "${EXT_DIR}"
cp "${COMPONENT_WASM}" "${EXT_DIR}/blast.wasm"
# Stage jsonfns too — DuckLink's core autoloads it and our sql/blast_macros.sql
# uses `to_json(...)` from it. Without this the CREATE MACRO calls that
# reference to_json fail (the aliases still work, but the sugar layer doesn't).
JSONFNS_WASM="${JSONFNS_WASM:-/Users/zacharywhitley/git/ducklink/artifacts/extensions/jsonfns.wasm}"
if [ -f "${JSONFNS_WASM}" ]; then
    cp "${JSONFNS_WASM}" "${EXT_DIR}/jsonfns.wasm"
    log "staged jsonfns from ${JSONFNS_WASM}"
else
    log "WARN: ${JSONFNS_WASM} not found; the to_json helpers in sql/blast_macros.sql will error at CREATE time (blastn_of alias still works)."
fi

# ---------------------------------------------------------------------------
# 3) Parse the GenBank inputs (same helper as the mock-mode driver).
# ---------------------------------------------------------------------------
if ! compgen -G "${DATA_DIR}/*.gb" >/dev/null && \
   ! compgen -G "${DATA_DIR}/*.gbk" >/dev/null; then
    log "ERROR: no *.gb/*.gbk files under ${DATA_DIR}"
    exit 1
fi

log "parsing GenBank inputs"
python3 acceptance/gb_parse.py \
    --data-dir "${DATA_DIR}" \
    --out-dir  "${BUILD_DIR}" \
    --per-genome "${PER_GENOME}" \
    --max-len  "${MAX_CDS_LEN}"

# ---------------------------------------------------------------------------
# 4) Build one big self-contained SQL script that DuckLink runs in a single
#    invocation. Everything -- queries, subjects, genome_features, the three
#    view definitions from sql/, and the assertion queries -- gets inlined
#    because file reads inside the DuckLink CLI wasm cannot reach paths on
#    the host filesystem without preopen fiddling that has not proven
#    reliable in testing.
# ---------------------------------------------------------------------------
log "assembling self-contained pipeline SQL"
PIPELINE_SQL="${BUILD_DIR}/ducklink-pipeline.sql"

python3 <<PY > "${PIPELINE_SQL}"
import json
from pathlib import Path

repo = Path("${REPO_ROOT}")
build = repo / "${BUILD_DIR}"
sql_dir = repo / "sql"

def sql_str(s: str) -> str:
    """Escape a Python str for use as a DuckDB single-quoted literal."""
    return "'" + s.replace("'", "''") + "'"

queries_json = (build / "queries.json").read_text()
subjects_json = (build / "subjects.json").read_text()
options_json = json.dumps({"evalue_max": ${EVALUE_MAX}})

print("LOAD blast;")
print("LOAD jsonfns;   -- gives us to_json() for the sql/blast_macros.sql sugar")
print()

# genome_features from the TSV, inlined as VALUES.
features_tsv = (build / "genome_features.tsv").read_text().strip().splitlines()
header = features_tsv[0].split("\t")
rows = [line.split("\t") for line in features_tsv[1:]]
print("CREATE TABLE genome_features (")
print("  genome_id TEXT, feature_key TEXT, feature_type TEXT,")
print("  start_position INTEGER, end_position INTEGER, strand INTEGER,")
print("  gene TEXT, product TEXT")
print(");")
print("INSERT INTO genome_features VALUES")
def row_sql(r):
    return "(" + ", ".join(sql_str(v) if not (v.lstrip('-').isdigit()) else v for v in r) + ")"
print(",\n".join(row_sql(r) for r in rows) + ";")
print()

# Load the blast macros (source-inline; the .read path can't reach ./sql).
print("-- blast_macros.sql")
print((sql_dir / "blast_macros.sql").read_text())
print()

# The actual BLASTN call, wrapped via the freshly-defined macros.
print("CREATE TABLE blast_hits AS")
print(f"SELECT * FROM blastn_of({sql_str(queries_json)}, {sql_str(subjects_json)}, {sql_str(options_json)});")
print()

# The three conservation views.
for name in ("best_hits", "gene_conservation", "summary"):
    print(f"-- {name}.sql")
    print((sql_dir / f"{name}.sql").read_text())
    print()

# Assertion rows.
print(r"""
SELECT 'ASSERT_HITS_NONZERO' AS assertion,
       CASE WHEN count(*) > 0 THEN 'OK' ELSE 'FAIL' END AS status,
       count(*) AS observed
FROM blast_hits;

SELECT 'ASSERT_ANY_CONSERVED' AS assertion,
       CASE WHEN sum(CASE WHEN n_conserved_in_genomes > 0 THEN 1 ELSE 0 END) > 0 THEN 'OK' ELSE 'FAIL' END AS status,
       sum(CASE WHEN n_conserved_in_genomes > 0 THEN 1 ELSE 0 END) AS observed
FROM gene_conservation_summary;

SELECT 'ASSERT_PCT_IN_RANGE' AS assertion,
       CASE WHEN count(*) FILTER (conservation_pct NOT BETWEEN 0 AND 100) = 0 THEN 'OK' ELSE 'FAIL' END AS status,
       count(*) FILTER (conservation_pct NOT BETWEEN 0 AND 100) AS observed
FROM gene_conservation_summary;
""")
PY

log "pipeline SQL: ${PIPELINE_SQL} ($(wc -c < "${PIPELINE_SQL}") bytes)"

# ---------------------------------------------------------------------------
# 5) Run the whole pipeline in one DuckLink invocation.
# ---------------------------------------------------------------------------
SQL_LOG="${BUILD_DIR}/ducklink-pipeline.out"
log "running DuckLink"
"${DUCKLINK_BIN}" --extensions-dir "${EXT_DIR}" -- :memory: < "${PIPELINE_SQL}" \
    2>&1 \
    | grep -vE '^\[wasi-fs\]|^\[autoload\]|^\[resolver\]|^\[extension|^\[duckdb|^\[table-stream' \
    | tee "${SQL_LOG}"

# ---------------------------------------------------------------------------
# 6) Verify the assertions.
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
    # Column layout: | assertion | status | observed |
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
    log "ACCEPTANCE-DUCKLINK FAILED"
    exit 1
fi
log "ACCEPTANCE-DUCKLINK PASSED"
