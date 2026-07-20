# sql/

DuckDB SQL that consumes the `blast_hits` rows emitted by the
`bio:blast` wasm component and the `genome_features` annotation table,
and computes Synphage's conservation-analysis outputs.

Pure SQL. Load the files in order (each defines a view the next uses):

```sql
.read sql/best_hits.sql
.read sql/gene_conservation.sql
.read sql/summary.sql
```

## Files

| File                          | Defines                       | Depends on                                    |
|-------------------------------|-------------------------------|-----------------------------------------------|
| `best_hits.sql`               | view `best_hits`              | table `blast_hits`, table `genome_features`   |
| `gene_conservation.sql`       | view `gene_conservation`      | view `best_hits`, table `genome_features`     |
| `summary.sql`                 | view `gene_conservation_summary` | view `gene_conservation`                   |

## Assumed schemas

`blast_hits` — as emitted by the `bio:blast` wasm component (13 columns):

```
query_key         TEXT
subject_key       TEXT
query_start       UINTEGER   -- 1-indexed inclusive
query_end         UINTEGER
subject_start     UINTEGER
subject_end       UINTEGER
strand            TEXT       -- 'plus' | 'minus'
identity_count    UINTEGER
alignment_length  UINTEGER
percent_identity  DOUBLE
bit_score         DOUBLE
raw_score         DOUBLE
evalue            DOUBLE
```

`genome_features` — annotation table (populated by the future
`bio:genome-format` component's `genbank_features()` TVF; hand-loaded from
GFF/GenBank meanwhile):

```
genome_id       TEXT
feature_key     TEXT   -- unique across all genomes; same key BLAST used
feature_type    TEXT   -- 'CDS' | 'gene' | ...
start_position  INTEGER
end_position    INTEGER
strand          INTEGER
gene            TEXT
product         TEXT
```

## What was ported from Synphage

Synphage's SQL, all of it, lives in
[`synphage/sql/`](https://github.com/vestalisvirginis/synphage/tree/main/synphage/sql):

- **`parse_blastn.sql`** — pulls flat rows out of BLAST's JSON via
  `read_json_auto` + `->>` path extraction, always taking `hits[0]` and
  `hsps[0]`. **Not ported.** Our component emits typed hit rows directly,
  so the JSON stage disappears. Its implicit "first hit wins" behaviour
  survives as the explicit ranking in `best_hits.sql`.
- **`gene_presence.sql`** — one `LEFT JOIN` of the concatenated BLAST
  results to the locus table, keyed on `key`. **Ported** into
  `gene_conservation.sql`, expanded to enumerate
  `(query_gene, subject_genome)` pairs so genes missing from a genome show
  up as rows with `conserved = false` rather than as absent rows. The
  Synphage-shaped one-line join is still recoverable:
  `SELECT * FROM genome_features f LEFT JOIN best_hits b ON b.query_key = f.feature_key`.

The Synphage pipeline's remaining computation (Polars concat of per-file
parquets into `blastn_summary.parquet`) doesn't exist here — the
component emits one relation, so consolidation is a no-op.

## Semantic choices worth flagging

- **Best-hit ordering:** `bit_score DESC, evalue ASC, subject_key ASC`.
  The first two match the DESIGN.md acceptance test; the third is a
  deterministic tie-break we add so best-hit selection is reproducible.
- **Subject-genome lookup:** derived from `genome_features.genome_id`
  keyed on `subject_key`. Synphage carries it implicitly (one BLAST run
  per file); we derive it from the data.
- **Self-hits:** dropped in `gene_conservation.sql`
  (`query_genome <> subject_genome`). Every gene trivially blasts against
  itself; leaving those rows in would peg `conservation_pct` at 100 for
  every gene.
- **Feature-type filter:** `WHERE feature_type = 'CDS'` in
  `gene_conservation.sql`, matching Synphage's protein/CDS focus. Remove
  the clause for RNA/pseudogene workflows.

## Future work

- Ranking-strategy switches (best-per-genome vs. best-per-subject,
  best-forward-hit-only, reciprocal best hits) as parameterised macros.
- A ChromSyn-style collinear-block macro that walks best_hits ordered by
  `query_start` and emits contiguous runs.
- A per-pair pivot view (`gene_conservation_matrix`) with one column per
  subject genome — useful for the synteny renderer's presence table.
