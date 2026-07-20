-- pipeline.sql
--
-- Runs the full acceptance query: loads hits.tsv + genome_features.tsv,
-- materialises the three views from ../sql/, and prints the summary.
--
-- Invoked by acceptance/run.sh with the acceptance/build/ inputs already
-- on disk. The paths are read from environment-substituted variables the
-- driver script sets via DuckDB's `SET VARIABLE` mechanism.

.echo off
.headers on
.mode column
.nullvalue NULL

-- 1) Load the raw hits emitted by the blast wasm component. The TSV has
--    exactly the 13 columns of the WIT `hit` record; DuckDB infers types
--    from the header but we're explicit about a few (evalue can be tiny;
--    percent_identity is fractional).
CREATE OR REPLACE TABLE blast_hits AS
SELECT
    query_key         :: VARCHAR AS query_key,
    subject_key       :: VARCHAR AS subject_key,
    query_start       :: UINTEGER AS query_start,
    query_end         :: UINTEGER AS query_end,
    subject_start     :: UINTEGER AS subject_start,
    subject_end       :: UINTEGER AS subject_end,
    strand            :: VARCHAR  AS strand,
    identity_count    :: UINTEGER AS identity_count,
    alignment_length  :: UINTEGER AS alignment_length,
    percent_identity  :: DOUBLE   AS percent_identity,
    bit_score         :: DOUBLE   AS bit_score,
    raw_score         :: DOUBLE   AS raw_score,
    evalue            :: DOUBLE   AS evalue
FROM read_csv(
    getvariable('hits_path'),
    delim='\t',
    header=true,
    auto_detect=true
);

-- 2) Load the CDS annotation table. Column order matches sql/README.md.
CREATE OR REPLACE TABLE genome_features AS
SELECT
    genome_id        :: VARCHAR AS genome_id,
    feature_key      :: VARCHAR AS feature_key,
    feature_type     :: VARCHAR AS feature_type,
    start_position   :: INTEGER AS start_position,
    end_position     :: INTEGER AS end_position,
    strand           :: INTEGER AS strand,
    gene             :: VARCHAR AS gene,
    product          :: VARCHAR AS product
FROM read_csv(
    getvariable('features_path'),
    delim='\t',
    header=true,
    auto_detect=true,
    nullstr=''
);

-- Sanity check the inputs.
SELECT 'blast_hits' AS tbl, count(*) AS n FROM blast_hits
UNION ALL
SELECT 'genome_features', count(*) FROM genome_features;

-- 3) Apply the shipped conservation SQL views (unchanged; loaded from the
--    canonical sql/ directory).
.read sql/best_hits.sql
.read sql/gene_conservation.sql
.read sql/summary.sql

-- 4) A quick per-genome overview of hits kept.
SELECT
    g.genome_id                                             AS query_genome,
    count(DISTINCT h.query_key)                             AS n_query_genes,
    count(*)                                                AS n_raw_hits,
    round(min(h.evalue)::DOUBLE,   6)                       AS min_evalue,
    round(max(h.percent_identity)::DOUBLE, 2)               AS max_pct_id
FROM blast_hits h
LEFT JOIN genome_features g ON g.feature_key = h.query_key
GROUP BY g.genome_id
ORDER BY g.genome_id;

-- 5) The main deliverable: per-gene conservation summary.
SELECT * FROM gene_conservation_summary ORDER BY query_genome, query_feature_key;

-- 6) Assertion rows -- the driver script parses these fixed markers.
SELECT
    'ASSERT_HITS_NONZERO' AS check_name,
    CASE WHEN (SELECT count(*) FROM blast_hits) > 0 THEN 'OK' ELSE 'FAIL' END AS status,
    (SELECT count(*) FROM blast_hits) AS observed;

SELECT
    'ASSERT_ANY_CONSERVED' AS check_name,
    CASE WHEN (SELECT count(*) FROM gene_conservation WHERE conserved) > 0 THEN 'OK' ELSE 'FAIL' END AS status,
    (SELECT count(*) FROM gene_conservation WHERE conserved) AS observed;

SELECT
    'ASSERT_PCT_IN_RANGE' AS check_name,
    CASE WHEN (
        SELECT count(*) FROM gene_conservation_summary
        WHERE conservation_pct IS NOT NULL
          AND (conservation_pct < 0 OR conservation_pct > 100)
    ) = 0 THEN 'OK' ELSE 'FAIL' END AS status,
    (
        SELECT count(*) FROM gene_conservation_summary
        WHERE conservation_pct IS NOT NULL
          AND (conservation_pct < 0 OR conservation_pct > 100)
    ) AS observed;
