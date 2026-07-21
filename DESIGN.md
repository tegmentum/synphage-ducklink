# synphage-ducklink design

## Intent

Reproduce [Synphage](https://github.com/vestalisvirginis/synphage)'s comparative-genomics workflow as a set of WebAssembly components that DuckDB, through DuckLink, calls as ordinary SQL table functions.

The point is not to repackage Synphage. It's to show that its four conceptual stages — acquire, parse, align, render — can each be exposed as a **relation**, so the whole pipeline becomes queryable rather than orchestrated.

The target user experience:

```sql
WITH genes AS (
    SELECT * FROM genbank_features('data/*.gb')
),
matches AS (
    SELECT *
    FROM blastn(
        queries  => (SELECT feature_key AS key, nucleotide_sequence AS data FROM genes),
        subjects => (SELECT feature_key AS key, nucleotide_sequence AS data FROM genes),
        scoring  => 'blastn-default',
        options  => {'evalue_max': 1e-5, 'max_target_seqs': 1}
    )
)
SELECT q.*, s.*, m.percent_identity
FROM genes q
LEFT JOIN matches m ON q.feature_key = m.query_key
LEFT JOIN genes s ON s.feature_key = m.subject_key;
```

No BLAST+ install. No Dagster server. No JSON on disk. Everything after `genbank_features()` is SQL.

## What the WIT declares

Three capabilities under one package (`tegmentum:bio@0.1.0`), each in its own file:

| Interface           | Owns                                                                    | Status         |
|---------------------|-------------------------------------------------------------------------|----------------|
| `sequence-search`   | BLASTN + BLASTP as one function with a scoring variant. Emits typed hit rows straight from the aligner — no BLAST JSON in sight. | First target. |
| `genome-format`     | GenBank parsing into four relations (records / features / qualifiers / sequences), replacing Biopython. | Draft. Later slice. |
| `synteny-renderer`  | Tracks + features + links in, SVG bytes out. The current Python renderer stays available in parallel; this is the pure-Wasm path. | Draft. Last slice. |

Four worlds are declared: one composite (`bio-components`) plus three single-capability worlds so components can ship independently. The initial `blast` component targets `sequence-search-only`.

## Why these shapes

### One `search` function, not `blastn` + `blastp`

BLASTN and BLASTP differ in exactly two places: scoring parameters and whether both strands are searched. Both are configuration, not interface. Encoding them as a `scoring` variant means:

- one Rust component satisfies both (shared aligner, shared Karlin-Altschul path, distinct scoring closures);
- adding new scoring schemes (BLASTX, TBLASTN, PAM matrices) is a WIT-additive change with no new function names;
- DuckLink can expose them as `blastn(...)`, `blastp(...)` table function aliases at the SQL surface — the naming lives in DuckLink's registration, not in the Wasm ABI.

### Unit variants for presets, data-carrying variants for overrides

`blastn-default` and `blastp-default` carry no data — they're a shorthand for "NCBI defaults, don't make me spell them out". `blastn(blastn-scoring)` and `blastp(blastp-scoring)` let power users override individual knobs without inventing per-parameter options-record fields.

The alternative — one big options struct with every knob optional — makes the "just give me the default" case verbose and the "which fields matter for blastn vs blastp" question un-answerable from the type alone.

### `list<T>` in, `list<T>` out (for now)

WIT's `stream<T>` is not in the stable component-model subset yet. `list<T>` is fine for the sizes Synphage-shape workflows produce (dozens of genomes, thousands of features each; hits are bounded by `max-target-seqs`). Migration path when streams stabilise:

- swap `list<sequence>` for `stream<sequence>`;
- swap `result<list<hit>, _>` for a `stream<hit>` plus a terminal error resource;
- DuckLink's Arrow-batch bridge maps cleanly onto both.

Interim, if throughput hurts: introduce a `sequence-source` / `hit-sink` resource pair with pull-based reads. That's a WIT-additive change too.

### Rank in SQL, not in the component

The component returns every hit surviving `options`. Best-hit selection (row_number over `query_key`, `subject_genome` ordered by bit_score / evalue) is a two-line CTE in DuckDB. Doing it in the component would:

- duplicate SQL that already reads well;
- fix a policy (best-per-subject-genome) that Synphage happens to use but is a domain choice;
- require the component to know what "subject genome" means, which it doesn't — it only knows subject keys.

Keeping ranking in SQL means the same component works for Synphage's best-hit workflow, ChromSyn's collinear-block workflow, and any future workflow that wants full hit rows.

### Positions are 1-indexed and inclusive at both ends

Matches NCBI BLAST's output convention. Half-open zero-indexed would match Rust idiom but forces every downstream consumer to remember to add 1 before comparing to a `.gff` file. The one place the component uses 0-indexed math (rust-bio internals) is behind the WIT boundary.

### `hit`, not `match`

`match` is a Rust keyword; the generated bindings would need renaming anyway. `hit` is BLAST-idiomatic ("high-scoring pair", "hit"), one syllable, no conflict.

## What a DuckLink host is expected to do

For each `SELECT * FROM blastn(queries => TABLE q, subjects => TABLE s, ...)` call:

1. Materialise the `queries` and `subjects` relations into `list<sequence>` batches — the two required columns are `key: VARCHAR` and `data: VARCHAR`.
2. Translate the SQL-level scoring argument (a struct literal or a preset string) into the WIT `scoring` variant.
3. Translate the SQL `options` MAP into `search-options`.
4. Invoke `sequence-search.search`.
5. Unpack the returned `list<hit>` into a DuckDB result set with columns matching the `hit` record's field names.

Errors from the `search-error` variant become SQL-visible errors, tagged with the variant name so users can `CASE` on them.

Everything else — the `WITH`s, joins, filters, best-hit CTEs, `COPY … TO 'x.parquet'` — is stock DuckDB.

## DuckLink integration

The `blast` compound world exports DuckLink's three dispatch interfaces alongside `sequence-search`, and the component registers two table functions in its `load()`:

```
blastn(queries LIST(STRUCT(key VARCHAR, data VARCHAR)),
       subjects LIST(STRUCT(key VARCHAR, data VARCHAR)),
       options STRUCT(evalue_max DOUBLE, max_target_seqs INTEGER, min_identity DOUBLE))
blastp(...)  -- same signature
```

Both dispatch to the same `crate::run_search` with different hardcoded scoring presets (`Scoring::BlastnDefault`, `Scoring::BlastpDefault`). Tunable scoring is future work — the first SQL surface stays clean.

Design choices:

- **Args are `list<struct>` via `complex(json)`.** DuckDB TVFs can't take subquery arguments and DuckLink's `logicaltype` enum is frozen (STRUCT/LIST ride the `complex(complexvalue { type-expr, json })` escape hatch). Callers wrap subqueries with `list(struct_pack(...))` at the SQL layer.
- **13 output columns match the `hit` WIT record.** Same names, same order, `s32` positions become `Uint32` on the DuckLink side (positions are always ≥ 1).
- **Streaming cursor, eager materialisation.** `call_table_open_filtered` runs the full search up-front, stores hits in a `HashMap<cursor_id, Vec<Hit>>`, and paginates through `call_table_next` in `max-rows` batches. Fine for Synphage-scale (dozens of genomes × thousands of features); if we ever need genome-scale streaming, the cursor can pull one query at a time from a lazy iterator.
- **Filter pushdown ignored on first pass.** Per DuckLink's freeze policy, ignoring pushed filters is still correct — the core re-checks above the scan. First slice punts on the optimisation.
- **`callback-dispatch` returns `Unsupported`** for every scalar/aggregate/cast entry point. We're a table-only extension.

## Deferred

- **Streaming.** Waiting on WIT `stream<T>` stabilisation. `list<T>` is enough for the acceptance test.
- **Additional matrices / scoring schemes.** BLOSUM62 first; others slot in as string names in `blastp-scoring.matrix`.
- **Tunable scoring at the DuckLink surface.** SQL callers today pick `blastn` or `blastp` (each a preset); custom (reward, penalty, gap-open, gap-extend) would surface as another STRUCT arg.
- **Filter pushdown.** Recognise `WHERE evalue < X`, `WHERE percent_identity >= Y`, `WHERE query_key = ...` and pre-filter in the component before rank/emit.
- **Palette / theme handling.** Currently baked into the renderer's `feature.colour` / `link.colour`. May extract to a `palette` capability once we have more than one renderer.
- **Multi-record GenBank edge cases.** The `parsed-genbank` shape is right for the common case; nested locations and joined features need real data to iron out.
- **End-to-end runtime test in DuckLink.** Requires `ducklink_{cli,core,loader}.wasm` built in the sibling `../ducklink/` checkout. Component ABI already verified against the WIT via `wasm-tools`.

## Acceptance test

The demo the project should be judged on:

DuckLink is a native DuckDB community extension (crate: [`ducklink-extension`](https://github.com/tegmentum/ducklink-extension), currently v5.0.0, WIT contract `duckdb:extension@4.0.0`). It embeds `wasmtime` inside a running DuckDB and bridges each loaded wasm component's registered functions into DuckDB's catalog.

**Today (pre-community-merge):** the extension isn't yet in `duckdb/community-extensions`, so users build the `.duckdb_extension` locally and `LOAD` it from a filesystem path. Components are named up front via the `DUCKLINK_COMPONENTS` env var (a `:`-separated list of `name=path` entries) — DuckDB's extension-load lifecycle registers catalog entries once at LOAD time, so the components have to be known before the LOAD runs:

```sh
DUCKLINK_COMPONENTS=blast=/path/to/blast.wasm \
    duckdb -unsigned -c "LOAD 'ducklink.duckdb_extension'; SELECT * FROM blastn(…);"
```

**In-SQL loading — `ducklink_load(<name>)` scalar:** ducklink-extension also registers a scalar function `ducklink_load('name')` that resolves the component by *catalog name* (not a filesystem path), fetches / verifies / caches the blob, and registers its functions:

```sql
LOAD 'ducklink.duckdb_extension';
SELECT ducklink_load('blast');   -- resolves 'blast' via the catalog
SELECT * FROM blastn(…);
```

DuckDB's own `LOAD <name>;` statement cannot be hijacked for a third-party name — `PhysicalLoad::GetDataInternal` calls `ExtensionHelper::LoadExternalExtension` unconditionally — which is why the entry point is an explicit function rather than another `LOAD` verb.

**Once ducklink lands in `duckdb/community-extensions`,** the bootstrap becomes:

```sql
LOAD ducklink FROM community;
SELECT ducklink_load('blast');
```

Two deployment models cover the same wasm layer:

- **DuckLink-as-DuckDB-extension** (above) — the ergonomic path for users already running native DuckDB. Native DuckDB loads `ducklink.duckdb_extension`, which embeds wasmtime and runs the wasm component inside the DuckDB process. A single portable component extends DuckDB on every platform without per-platform native extension builds.
- **Standalone `ducklink` binary** — a self-contained runner that wraps DuckDB-in-wasm + the wasm-component host in one CLI. Handy for containers, browsers, and CI. The CLI's `parse_load_names` preprocessor recognises `LOAD <name>;` inputs and rewrites them through the same wasm loading path, so under the standalone CLI `LOAD blast;` works too. `acceptance/run-ducklink.sh` uses this form because it needs zero user-side install for the smoke test — it stages `blast.wasm` under `acceptance/build/ext/` and runs `ducklink --extensions-dir acceptance/build/ext -- :memory: -c "LOAD blast; …"`.

The acceptance query (aspirational — see the caveats below the code):

```sql
LOAD ducklink FROM community;   -- once ducklink is merged into duckdb/community-extensions
SELECT ducklink_load('blast');

WITH genes AS (
    SELECT * FROM genbank_features('examples/*.gb')
),
matches AS (
    SELECT * FROM blastn(
        queries  => (SELECT feature_key AS key, sequence AS data FROM genes),
        subjects => (SELECT feature_key AS key, sequence AS data FROM genes),
        scoring  => 'blastn-default',
        options  => {'evalue_max': 1e-5, 'max_target_seqs': 1}
    )
),
ranked AS (
    SELECT *,
           row_number() OVER (
               PARTITION BY query_key, subject_genome
               ORDER BY bit_score DESC, evalue ASC
           ) AS rk
    FROM matches
),
conservation AS (
    SELECT q.*, s.*, r.percent_identity
    FROM genes q
    LEFT JOIN (SELECT * FROM ranked WHERE rk = 1) r
      ON q.feature_key = r.query_key
    LEFT JOIN genes s ON s.feature_key = r.subject_key
)
SELECT * FROM conservation;
```

**Caveats vs. what actually runs today** (see [`acceptance/run-ducklink.sh`](acceptance/run-ducklink.sh) for the working form):

- `SELECT * FROM genbank_features('examples/*.gb')` — replace with `genbank_scan('<inlined-genbank-text>')` today. DuckDB's TVF binder rejects both scalar subqueries and LATERAL in table-function argument positions, so the sugar `genbank_scan((SELECT content FROM read_text('*.gb')))` doesn't bind. The extension side is already positioned to accept it the moment that DuckDB binder rule lifts. Multi-file support is a DuckDB-side concatenation until then.
- `blastn(queries => (SELECT ...), subjects => (SELECT ...), scoring => 'blastn-default', options => {'evalue_max': 1e-5})` — replace with `blastn('<queries-json>', '<subjects-json>', '<opts-json>')` today. Same DuckDB binder rule, plus `Logicaltype::Complex("LIST(STRUCT(…))")` currently flattens to `VARCHAR[]` at DuckLink's DuckDB binder — the STRUCT payload is lost. Every DuckLink-loadable component in this repo takes JSON-string args instead. The `Duckvalue::Complex` code paths are still there for the day DuckLink learns to preserve richer type-expressions.

Compared against Synphage on the same inputs the pipeline should produce:

- identical best-hit selection;
- equivalent gene-presence percentages;
- reproducibly across Linux / macOS / Windows / browser (via the same `.wasm`);
- with no BLAST+ install and no Dagster.
