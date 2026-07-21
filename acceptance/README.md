# acceptance/

End-to-end acceptance test for the synphage-ducklink pipeline on real NCBI
phage genomes.

The tiny embedded samples in `hosts/blast-test/src/main.rs` prove the
component ABI works. This directory proves the shape and semantics of the
outputs make biological sense on inputs that resemble the ones Synphage
actually runs on.

## What it does

```
acceptance/data/*.gb            <- three real NCBI phage GenBank files
     │
     ▼  (Python)
acceptance/gb_parse.py          <- CDS + sequence extraction
     │
     ▼
acceptance/build/{queries.json, subjects.json, genome_features.tsv}
     │
     ▼  (wasmtime host)
hosts/blast-test/target/release/blast-test  \
    --queries ... --subjects ... --hits-out ... \
    --options-json '{"evalue_max":1e-5}'
     │
     ▼
acceptance/build/hits.tsv       <- BLASTN output emitted by the wasm component
     │
     ▼  (DuckDB)
acceptance/pipeline.sql
    - loads hits.tsv into `blast_hits`
    - loads genome_features.tsv into `genome_features`
    - .read sql/best_hits.sql
    - .read sql/gene_conservation.sql
    - .read sql/summary.sql
    - prints per-genome overview, per-gene conservation summary
    - emits three ASSERT_* rows
     │
     ▼
acceptance/build/pipeline.out   <- captured DuckDB output
     │
     ▼  (bash grep)
ACCEPTANCE PASSED / FAILED
```

## What we exercise

Same shape as the DESIGN.md acceptance query:

- `SELECT * FROM blastn(queries => ..., subjects => ..., options => {'evalue_max': 1e-5})`
  is driven through the component's DuckLink dispatch surface
  (`table-stream-dispatch::call_table_open_filtered` -> `call_table_next`
  -> `call_table_close`) using the wasmtime-based test host.
- The three shipped SQL views (`sql/best_hits.sql`,
  `sql/gene_conservation.sql`, `sql/summary.sql`) are `.read` from the
  canonical location, **unmodified**.
- All three of the design's "gene conservation" outputs -- per-gene
  presence-in-other-genomes, best-hit identity, and the flat rollup --
  are produced.

## Inputs

Three complete NCBI reference phage genomes, chosen because they cover a
useful mix of "closely related" and "distantly related":

| Accession    | Organism                       | Size  |
|--------------|--------------------------------|-------|
| `NC_001416`  | Enterobacteria phage lambda    | 49 kb |
| `NC_001604`  | Enterobacteria phage T7        | 40 kb |
| `NC_002371`  | Salmonella phage P22           | 42 kb |

Lambda and P22 are lambdoid siphoviruses -- close cousins with a
recognisable shared gene cassette (Nin proteins, cIII / c3, arc/repressor).
T7 is a podovirus from a different family -- we expect almost no
cross-hits to Lambda or P22 at BLASTN's nucleotide sensitivity.

The `.gb` files were fetched from NCBI eutils:

```sh
for acc in NC_001416 NC_001604 NC_002371; do
    curl -sSfL -o "acceptance/data/${acc}.gb" \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acc}&rettype=gb&retmode=text"
done
```

If they're missing, re-run that loop.

## What gets BLASTed

For tractability the driver samples up to 8 CDS per genome, preferring
CDS between 150 and 900 bp. That's 3 genomes x 8 CDS = 24 sequences,
alignedgeneviaanall-vs-all pairwise Smith-Waterman inside the component
(24 x 24 = 576 pairs x 2 strands = up to 1152 raw hits before
evalue filtering).

Full-genome-vs-full-genome pairwise SW would allocate an ~1.6-billion
cell DP matrix per pair; not feasible for a wasm-runtime test. Full
CDS-vs-CDS on all annotated CDS (~200 per genome) would allocate 40 000
DP matrices at up to ~4 000 x ~4 000 cells; several minutes in the
current wasmtime host. The 8-per-genome sample runs in about 0.5 seconds
and still surfaces the biologically expected Nin / cIII shared cassette
between lambda and P22.

Increase the sample size via env vars:

```sh
PER_GENOME=20 MAX_CDS_LEN=1500 ./acceptance/run.sh
```

Tightening the e-value threshold:

```sh
EVALUE_MAX=1e-10 ./acceptance/run.sh
```

## Assertions

`acceptance/pipeline.sql` emits three tagged rows the shell script
verifies:

- `ASSERT_HITS_NONZERO` -- at least one hit survived the e-value cutoff.
  Any real phage-vs-phage comparison should produce many.
- `ASSERT_ANY_CONSERVED` -- at least one query gene has a best-hit in one
  other genome. With three lambdoid-adjacent genomes and default lambda
  the count is ~8; if this drops to zero something is wrong with either
  the alignment or the join.
- `ASSERT_PCT_IN_RANGE` -- every `conservation_pct` value produced by
  `gene_conservation_summary` is in `[0, 100]` (or NULL). This catches
  division-by-zero regressions.

## Expected output

The full DuckDB output is captured to `acceptance/build/pipeline.out`.
A representative run at the defaults (24 CDS, evalue 1e-5) prints:

```
ASSERT_HITS_NONZERO  OK  40
ASSERT_ANY_CONSERVED OK  8
ASSERT_PCT_IN_RANGE  OK  0
```

and shows lambda's `NinD` / `NinE` / `NinF` / `cIII` genes hitting P22's
`ninD` / `ninE` / `ninF` / `c3` at 76-95% identity, with T7 genes finding
no homologs -- consistent with T7's phylogenetic distance from the
lambdoid group.

## Invocation

Two drivers ship. Both prep inputs the same way and produce the same
assertion output; they differ in how the BLAST call is routed.

### Mock host (fast, no DuckLink build required)

```sh
./acceptance/run.sh
```

Steps:
1. builds the blast wasm component if missing (needs `cargo-component`);
2. builds the wasmtime test host at `hosts/blast-test/` if missing;
3. parses the GenBank inputs (needs `python3` >= 3.11 -- no third-party
   dependencies);
4. runs the wasm BLASTN scan through the wasmtime mock host, which
   reimplements just enough of DuckLink's `table-stream` protocol to
   exercise the component ABI;
5. runs the DuckDB pipeline via a native `duckdb` binary
   (needs `duckdb` on `PATH` -- tested against 1.5.4);
6. verifies the assertion markers and prints `ACCEPTANCE PASSED / FAILED`.

### Real DuckLink (the locked-in end-to-end)

```sh
./acceptance/run-ducklink.sh
```

Same pipeline, but the BLAST call goes through the *actual* DuckLink
CLI (`ducklink_{core,cli,loader}.wasm` wired together via `wac plug`,
loading `blast.wasm` as an extension). The blast + conservation +
assertion queries all run in one DuckLink invocation.

### Native DuckDB (the daily-driver target)

```sh
./acceptance/run-native-duckdb.sh
```

Uses the native `duckdb` CLI plus the `ducklink.duckdb_extension` native
community extension (currently built locally at
`../ducklink-extension/build/release/`; once ducklink lands in
`duckdb/community-extensions` this becomes `LOAD ducklink FROM
community;`). This is the shape an ordinary DuckDB user will run.

Uses DuckDB's built-in `json` extension (autoloaded) — no `jsonfns`
needed here. Builds the queries JSON as
`to_json(list(struct_pack(...)))`.

Everything else is the same as `run-sql-only.sh`: `read_text` +
`SET VARIABLE` + `getvariable` + `genbank_scan` + `blastn` + conservation
views + assertions. 40 hits, 8 conserved. `ACCEPTANCE-NATIVE-DUCKDB
PASSED`.

### Pure SQL under the DuckLink wasm CLI

```sh
./acceptance/run-sql-only.sh
```

The whole pipeline (`read → parse → BLAST → conservation → summary →
assertions`) expressed as one DuckLink SQL script,
`acceptance/pipeline-sql-only.sql`. Nothing is inlined out-of-band:
`read_text('acceptance/data/NC_00….gb')` feeds a `SET VARIABLE`, and the
BLAST call takes a `getvariable('q_json')` populated by an in-SQL
`string_agg + json_quote` builder. This is what the design doc meant by
"the whole pipeline is one SQL query" — the two DuckDB-side rules
against subqueries in TVF argument positions are the reason we can't
compose the whole thing as literally one statement, but a single script
with variables and CTEs is as close as DuckDB permits today, and
substantially the same shape once those rules lift.

Also loads DuckLink's `jsonfns` extension for `json_quote` (JSON string
escaping) and `json_extract_string` (pulling `gene`/`product` out of the
qualifier map). No native DuckDB `json` extension needed — the compile-
time flag disables external native extensions in this build.

Prerequisites (see the repo `README.md` for the build recipe):
- `/Users/zacharywhitley/git/ducklink/target/release/ducklink` (native host)
- `/Users/zacharywhitley/git/ducklink/target/wasm32-wasip2/release/ducklink_{core,cli,loader}.wasm`
- `/Users/zacharywhitley/git/ducklink/artifacts/extensions/jsonfns.wasm`
  (auto-staged; provides `to_json` used by `sql/blast_macros.sql`)

Override paths via env vars if your DuckLink lives elsewhere:

```sh
DUCKLINK_BIN=/path/to/ducklink \
JSONFNS_WASM=/path/to/jsonfns.wasm \
    ./acceptance/run-ducklink.sh
```

Exit code 0 means all assertions passed; non-zero means at least one
assertion failed or a build step errored -- same shape as
`run.sh`.

## Layout

```
acceptance/
├── README.md               <- this file
├── run.sh                  <- mock-host driver (wasmtime + hosts/blast-test)
├── run-ducklink.sh         <- real-DuckLink-wasm driver (Python parses inputs)
├── run-sql-only.sh         <- pure-SQL driver on the DuckLink wasm CLI
├── run-native-duckdb.sh    <- pure-SQL driver on NATIVE DuckDB + ducklink extension  ★
├── gb_parse.py             <- GenBank -> CDS extractor (used by run.sh / run-ducklink.sh)
├── pipeline.sql            <- DuckDB pipeline for run.sh (native duckdb)
├── pipeline-sql-only.sql   <- DuckLink-CLI script for run-sql-only.sh
├── pipeline-native.sql     <- Native-DuckDB script for run-native-duckdb.sh  ★
├── data/                   <- input GenBank files (NC_001416, NC_001604, NC_002371)
├── wit-view/               <- curated copy of wit/ used by the test-host bindgen
│                              (isolated from parallel-track wit/ churn)
└── build/                  <- generated artifacts
```

★ run-native-duckdb.sh is the target daily-driver — the shape ordinary
DuckDB users will invoke once ducklink is in duckdb/community-extensions.

## Caveats

- `gb_parse.py` is a deliberately-narrow GenBank subset parser -- enough
  for well-formed NCBI reference genomes. Real fields with unusual
  location grammar (nested `join(complement(join(...)))`, remote
  references) may parse imperfectly. The future
  `components/genome-format/` DuckLink bridge will replace this driver
  script; this file exists so the acceptance test doesn't block on that.
- Only BLASTN is exercised (BLASTP would need protein sequences from
  `/translation="..."` qualifiers, which the parser also extracts on
  request but we don't invoke here).
- `acceptance/wit-view/` mirrors the subset of `wit/` needed to bindgen
  the compound `blast` world. It's copied at test-authoring time; if the
  canonical `wit/` files change materially, refresh it via `cp` from
  `wit/`.
- The 8-CDS-per-genome sample deliberately biases towards short-medium
  CDS to keep the pairwise SW cost bounded. This does bias the biology:
  the very-long capsid / tail-fibre / polymerase genes -- often the most
  conserved -- are not included in the default run. Raise `PER_GENOME` /
  `MAX_CDS_LEN` if you want to probe them.
