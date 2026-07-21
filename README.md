# synphage-ducklink

Rebuilding [Synphage](https://github.com/vestalisvirginis/synphage)'s comparative-genomics workflow as WebAssembly components that DuckDB — via DuckLink — calls as ordinary SQL table functions.

Not a fork. Not a wrapper. A reimplementation of the *shape*: the same scientific method, exposed as relations instead of a Python + Dagster + BLAST+ pipeline.

## The contrast

| Synphage today                                    | synphage-ducklink                          |
|---------------------------------------------------|--------------------------------------------|
| `pip install` + `apt install ncbi-blast+`         | (once merged) `LOAD ducklink FROM community; SELECT ducklink_load('blast');` — see below for today's local-build path |
| Dagster server, four sequential jobs              | One SQL query                              |
| BLAST JSON on disk, parsed via nested paths       | Typed hit rows straight from the aligner   |
| Biopython + Polars + Pandas glue                  | DuckDB relations end to end                |
| Linux/macOS/Windows setup differences             | Same `.wasm` everywhere (incl. browser)    |

The scientific outputs — best-hit selection, gene conservation, synteny diagrams — stay recognisably Synphage. What disappears is the machinery around them.

## Status

- **WIT drafted** — `wit/` defines three biological capabilities (`sequence-search`, `genome-format`, `synteny-renderer`) plus a compound world `blast` that also exports DuckLink's dispatch surface. DuckLink's own WIT is vendored under `wit/deps/duckdb-extension/`.
- **Three components built,** all targeting `wasm32-wasip2`:
  - **`blast`** (~319 KB) — exports both `tegmentum:bio/sequence-search` AND `duckdb:extension/{guest,callback-dispatch,table-stream-dispatch}`. Registers `blastn` and `blastp` DuckLink table functions.
  - **`genome-format`** (~106 KB) — handwritten GenBank parser, exports `tegmentum:bio/genome-format` only (no DuckLink surface yet).
  - **`synteny-renderer`** (~105 KB) — handwritten SVG generator, exports `tegmentum:bio/synteny-renderer` only.
- **Filter pushdown live in blast** — `evalue < X`, `percent_identity >= Y`, `query_key IN (…)`, `subject_key IN (…)`, `strand = 'plus'/'minus'` all recognised and folded into `SearchOptions` / early batch pruning / short-circuit. Contradictory clauses (empty key intersection, plus AND minus) trigger a zero-row cursor.
- **Standalone smoke test via `hosts/blast-test`** — a wasmtime-based Rust binary loads `blast.wasm`, mocks DuckLink's `table_stream` host import, and drives the full protocol (`guest.load` → `call_table_open_filtered` → `call_table_next` → `call_table_close`). Useful for iterating without the DuckLink base build.
- **Real end-to-end under DuckLink verified.** With `ducklink_{core,cli,loader}.wasm` built and the native `ducklink` host binary, `LOAD blast; SELECT … FROM blastn('[…]','[…]', NULL) WHERE evalue < 1e-6;` returns hits from a real DuckDB engine via wasm. Projection and filter pushdown both observed in the runtime logs (`projection=[6,12,0,1,9,10] nfilt=2`).
- **Real phage genomes, real biology.** Running BLASTN on 24 CDS from three NCBI phage references (lambda `NC_001416`, T7 `NC_001604`, P22 `NC_002371`) through real DuckLink produces 40 hits with the expected phylogeny: lambda↔P22 shares 6 conserved genes each way (Nin/cIII cassette), T7 finds nothing in either (correct — T7 is a distant podovirus). Full pipeline runs in ~1 s.
- **Three DuckLink-loadable components in the tree, all fully working end-to-end.**
  - `blast.blastn/blastp` — 40 hits across 24 CDS from three real NCBI phage genomes.
  - `synteny-renderer.render_synteny_svg` — produces 610 B of SVG on a tiny test model (`bytes_len=610`, `svg_len=610`); empty inputs correctly return zero rows.
  - `genome-format.genbank_scan(contents VARCHAR)` — parses all three phage genomes: NC_001416 → 285 features / 73 CDS, NC_001604 → 162 / 60, NC_002371 → 191 / 72.
- **Repeatable real-DuckLink acceptance test.** `acceptance/run-ducklink.sh` drives the full pipeline (24 CDS from lambda / T7 / P22 → `blastn` under real DuckLink → conservation views → assertions) in a single DuckLink invocation and prints `ACCEPTANCE-DUCKLINK PASSED` on success. The mock-host driver `acceptance/run.sh` remains for fast iteration when the DuckLink base isn't built.
- **Pure-SQL acceptance test (no Python).** `acceptance/run-sql-only.sh` runs the whole pipeline — read → parse → BLAST → conservation → summary → assertions — as **one DuckLink SQL script** (`acceptance/pipeline-sql-only.sql`). Uses `SET VARIABLE x = (SELECT content FROM read_text(...))` + `getvariable('x')` to pass computed values through DuckDB's TVF-argument binder (which otherwise rejects subqueries). Produces the same 40-hits/8-conserved result. This is the closest DuckDB currently permits to "the whole scientific pipeline is one query"; the remaining plural-statement shape is a DuckDB binder rule, not intrinsic to the design.
- **Native DuckDB variant.** `acceptance/run-native-duckdb.sh` runs the same pure-SQL pipeline through the native `duckdb` CLI + the `ducklink.duckdb_extension` native community extension (no DuckLink wasm CLI). This is the daily-driver shape once ducklink lands in `duckdb/community-extensions`. Uses DuckDB's built-in `to_json(list(struct_pack(...)))` for the queries payload — no `jsonfns` extension needed. Same 40 hits, same conservation, `ACCEPTANCE-NATIVE-DUCKDB PASSED`.
- **Thin SQL sugar** at `sql/blast_macros.sql` (`blastn_of`, `blastp_of`). Table macros over the three JSON string args, with default `opts_json := NULL`. Fuller sugar (`blastn(TABLE genes)`) is blocked by DuckDB's binder rejecting subqueries in TVF-arg positions — the moment that lifts, the natural signature drops in without extension-side changes.
- **Two DuckLink-shape workarounds now understood and documented,** both applied uniformly across the three components:
  - **`Logicaltype::Complex("LIST(STRUCT(…))")` flattens to `VARCHAR[]`** at the DuckDB binder in DuckLink 4.0.0 — the STRUCT payload is lost. All three components now take VARCHAR JSON args instead; the `Duckvalue::Complex` code path is still there for when DuckLink learns to preserve richer type-expressions.
  - **`std::fs::read` inside a component can't reach files** — WASI preopens (`ducklink --dir HOST::GUEST`) don't thread through to extension instances, and DuckLink's `files` interface is for replacement-scan registration only, not general file reads. Idiomatic fix: extension takes file *content*, not path. `genome-format.genbank_scan(contents VARCHAR)` accepts raw GenBank text; users would compose with DuckDB's native `read_text()` if TVF-arg subqueries were allowed (currently a separate DuckDB binder limitation — literal inlining works today, subquery-composition becomes daily-driver the moment that binder rule lifts).
- **Synphage conservation SQL ported** — `sql/best_hits.sql`, `sql/gene_conservation.sql`, `sql/summary.sql`. Run against DuckDB 1.5.4 fixture data cleanly.

### DuckLink caveat: complex arg types

DuckLink 4.0.0's `Logicaltype::Complex("LIST(STRUCT(...))")` currently flattens to `VARCHAR[]` at the DuckDB binder — the STRUCT payload is lost. Our first cut used the natural LIST(STRUCT) shape and hit a type-signature mismatch. Working around it: `blastn` now takes three VARCHAR args carrying JSON payloads. Consumers write:

```sql
LOAD blast;
SELECT * FROM blastn(
    '[{"key":"q1","data":"ACGT..."}]',
    '[{"key":"s1","data":"ACGT..."}]',
    NULL
);
```

A DuckDB SQL macro that hides the `to_json()` wrapping is the natural next sugar, and the parser still accepts `Duckvalue::Complex` so this switch is reversible once DuckLink learns to preserve arbitrary type-expressions.

See [`DESIGN.md`](DESIGN.md) for the rationale behind the WIT shapes and the target demo query.

## Layout

```
synphage-ducklink/
├── README.md, DESIGN.md
├── wit/
│   ├── world.wit             # 4 biological-only worlds + compound `blast` world
│   ├── sequence-search.wit   # BLASTN + BLASTP (implemented)
│   ├── genome-format.wit     # GenBank -> four relations (implemented, biology-only)
│   ├── synteny-renderer.wit  # tracks/features/links -> SVG (implemented, biology-only)
│   └── deps/duckdb-extension/  # vendored DuckLink WIT (52 files)
├── components/
│   ├── blast/                # sequence-search + DuckLink dispatch
│   │   └── src/{lib,align,scoring,strand,filter,ducklink,pushdown,bindings}.rs
│   ├── genome-format/        # handwritten GenBank parser (biology-only)
│   │   └── src/{lib,parser,location,model,bindings}.rs
│   └── synteny-renderer/     # handwritten SVG renderer (biology-only)
│       └── src/{lib,render,bindings}.rs
├── hosts/
│   └── blast-test/           # standalone wasmtime host for end-to-end validation
├── examples/
│   └── tiny-blastn/          # small sample FASTA files for the test-host
└── sql/
    ├── best_hits.sql         # per-(query, subject-genome) ranking view
    ├── gene_conservation.sql # (query gene × subject genome) grid + LEFT JOIN back to features
    ├── summary.sql           # per-gene rollup: n_other_genomes, conservation_pct, avg identity
    └── README.md             # what each file computes + assumed schemas
```

## Roadmap

1. ✅ **BLASTN + BLASTP wasm component** targeting `sequence-search-only`.
2. ✅ **DuckLink table-function registration** — `blastn(...)`, `blastp(...)` at the SQL surface.
3. ✅ **Filter pushdown** — evalue ceiling, identity floor, key restrictions, strand, with short-circuit.
4. ✅ **Standalone end-to-end smoke test** via `hosts/blast-test` (validates dispatch surface without DuckLink base build).
5. ✅ **Conservation SQL** ported into `sql/` — best-hit ranking, gene-conservation grid, per-gene summary.
6. ✅ **`genome-format` wasm component** — handwritten GenBank parser exporting biology-only.
7. ✅ **`synteny-renderer` wasm component** — handwritten SVG generator exporting biology-only.
8. **DuckLink base build + real E2E** — build `ducklink_{cli,core,loader}.wasm` in `../ducklink/`, then run `SELECT * FROM blastn(…)` under real DuckLink. Requires the sibling repo's substantial build.
9. **DuckLink bridges for `genome-format` and `synteny-renderer`** — add the compound world + dispatch impls the same way `blast` got them, so `SELECT * FROM genbank_scan('*.gb')` and `SELECT render_svg(…)` work.
10. **Acceptance test** — run the demo query in `DESIGN.md` on real Synphage example genomes; compare outputs to the reference Synphage pipeline.

## Related

- [`synphage`](https://github.com/vestalisvirginis/synphage) — the workflow this reimplements.
- [`scry-webfunctions-demo`](../scry-webfunctions-demo) — the rust-bio BLASTP wasm component this project extends.
- [`webfunction-wit`](../webfunction-wit), [`stardog-webfunction-wit`](../stardog-webfunction-wit) — the SPARQL-shaped WIT the BLASTP component currently targets. The `tegmentum:bio` package in this repo is the DuckLink-shaped sibling.
