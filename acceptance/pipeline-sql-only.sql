-- Pure-SQL synphage-ducklink pipeline. No Python driver. Runs in one
-- DuckLink invocation. Prose lives in acceptance/run-sql-only.sh header
-- and in this file's git blame: the DuckLink CLI parser dislikes
-- multi-line -- comments with parens/backticks, so we keep annotations
-- terse and single-line here.

LOAD blast;
LOAD genome_format;
LOAD jsonfns;

-- Read three GenBank files into variables via read_text.
SET VARIABLE gb_lambda = (SELECT content FROM read_text('acceptance/data/NC_001416.gb'));
SET VARIABLE gb_t7     = (SELECT content FROM read_text('acceptance/data/NC_001604.gb'));
SET VARIABLE gb_p22    = (SELECT content FROM read_text('acceptance/data/NC_002371.gb'));

-- Parse each via genbank_scan(getvariable(...)) and union.
CREATE TABLE all_features AS
    SELECT * FROM genbank_scan(getvariable('gb_lambda'))
    UNION ALL
    SELECT * FROM genbank_scan(getvariable('gb_t7'))
    UNION ALL
    SELECT * FROM genbank_scan(getvariable('gb_p22'));

-- Sample the 8 shortest CDS per genome, length 150..900.
CREATE TABLE gene_sample AS
WITH cds AS (
    SELECT accession, feature_index, feature_type, start_position, end_position,
           strand, qualifiers_json, sequence
    FROM all_features
    WHERE feature_type = 'CDS'
      AND sequence <> ''
      AND length(sequence) BETWEEN 150 AND 900
),
ranked AS (
    SELECT *, row_number() OVER (
        PARTITION BY accession ORDER BY length(sequence), feature_index
    ) AS rk
    FROM cds
)
SELECT * FROM ranked WHERE rk <= 8;

-- Build the queries JSON in SQL via string_agg + json_quote.
SET VARIABLE q_json = (
    SELECT '[' || string_agg(
        '{"key":' || json_quote(concat(accession, ':', feature_index)) ||
        ',"data":' || json_quote(sequence) || '}',
        ','
    ) || ']'
    FROM gene_sample
);

-- BLAST all-vs-all.
CREATE TABLE blast_hits AS
SELECT * FROM blastn(
    getvariable('q_json'),
    getvariable('q_json'),
    '{"evalue_max": 1e-5}'
);

-- genome_features shaped the way sql/best_hits.sql expects.
CREATE TABLE genome_features AS
SELECT
    accession                                             AS genome_id,
    concat(accession, ':', feature_index)                 AS feature_key,
    feature_type,
    start_position::INTEGER                               AS start_position,
    end_position::INTEGER                                 AS end_position,
    strand,
    json_extract_string(qualifiers_json, '$.gene')        AS gene,
    json_extract_string(qualifiers_json, '$.product')     AS product
FROM gene_sample;

-- Conservation views inlined from sql/*.sql (kept in sync there).
CREATE OR REPLACE VIEW best_hits AS
WITH ranked AS (
    SELECT h.query_key, h.subject_key, f.genome_id AS subject_genome,
           h.query_start, h.query_end, h.subject_start, h.subject_end,
           h.strand, h.identity_count, h.alignment_length,
           h.percent_identity, h.bit_score, h.raw_score, h.evalue,
           row_number() OVER (
               PARTITION BY h.query_key, f.genome_id
               ORDER BY h.bit_score DESC, h.evalue ASC, h.subject_key ASC
           ) AS rk
    FROM blast_hits h
    LEFT JOIN genome_features f ON f.feature_key = h.subject_key
)
SELECT query_key, subject_key, subject_genome, query_start, query_end,
       subject_start, subject_end, strand, identity_count,
       alignment_length, percent_identity, bit_score, raw_score, evalue
FROM ranked WHERE rk = 1;

CREATE OR REPLACE VIEW gene_conservation AS
WITH subject_genomes AS (
    SELECT DISTINCT genome_id AS subject_genome FROM genome_features
),
query_x_subject AS (
    SELECT q.genome_id AS query_genome, q.feature_key AS query_feature_key,
           q.feature_type AS query_feature_type, q.gene AS query_gene,
           q.product AS query_product,
           q.start_position AS query_start_position,
           q.end_position AS query_end_position,
           q.strand AS query_gene_strand,
           sg.subject_genome
    FROM genome_features q CROSS JOIN subject_genomes sg
    WHERE q.feature_type = 'CDS' AND q.genome_id <> sg.subject_genome
)
SELECT x.query_genome, x.subject_genome,
       x.query_feature_key, x.query_feature_type, x.query_gene, x.query_product,
       x.query_start_position, x.query_end_position, x.query_gene_strand,
       b.subject_key AS subject_feature_key,
       s.gene AS subject_gene, s.product AS subject_product,
       s.start_position AS subject_start_position,
       s.end_position AS subject_end_position, s.strand AS subject_gene_strand,
       b.query_start AS hit_query_start, b.query_end AS hit_query_end,
       b.subject_start AS hit_subject_start, b.subject_end AS hit_subject_end,
       b.strand AS hit_strand, b.alignment_length, b.identity_count,
       b.percent_identity, b.bit_score, b.raw_score, b.evalue,
       (b.subject_key IS NOT NULL) AS conserved
FROM query_x_subject x
LEFT JOIN best_hits b
  ON b.query_key = x.query_feature_key
 AND b.subject_genome = x.subject_genome
LEFT JOIN genome_features s ON s.feature_key = b.subject_key;

CREATE OR REPLACE VIEW gene_conservation_summary AS
SELECT query_genome, query_feature_key, query_gene, query_product,
       count(*) AS n_other_genomes,
       count(*) FILTER (WHERE conserved) AS n_conserved_in_genomes,
       round(100.0 * count(*) FILTER (WHERE conserved)
             / nullif(count(*), 0), 2) AS conservation_pct,
       avg(percent_identity) FILTER (WHERE conserved) AS avg_percent_identity,
       min(percent_identity) FILTER (WHERE conserved) AS min_percent_identity,
       max(percent_identity) FILTER (WHERE conserved) AS max_percent_identity,
       avg(evalue) FILTER (WHERE conserved) AS avg_evalue
FROM gene_conservation
GROUP BY query_genome, query_feature_key, query_gene, query_product
ORDER BY query_genome, query_feature_key;

.print === per-genome overview ===
SELECT query_genome,
       count(*)                                    AS query_genes,
       count(DISTINCT subject_genome)              AS n_other_genomes,
       count(*) FILTER (WHERE conserved)           AS conservation_events,
       round(avg(percent_identity) FILTER (WHERE conserved), 2) AS avg_id
FROM gene_conservation
GROUP BY query_genome ORDER BY query_genome;

.print === conserved genes ===
SELECT query_genome, query_feature_key, query_gene, query_product,
       n_conserved_in_genomes, round(avg_percent_identity, 2) AS avg_id
FROM gene_conservation_summary
WHERE n_conserved_in_genomes > 0
ORDER BY query_genome, query_feature_key;

.print === assertions ===

-- Assertion rows (grepped by the shell driver).
SELECT 'ASSERT_HITS_NONZERO' AS assertion,
       CASE WHEN count(*) > 0 THEN 'OK' ELSE 'FAIL' END AS status,
       count(*) AS observed
FROM blast_hits;

SELECT 'ASSERT_ANY_CONSERVED' AS assertion,
       CASE WHEN sum(CASE WHEN n_conserved_in_genomes > 0 THEN 1 ELSE 0 END) > 0
            THEN 'OK' ELSE 'FAIL' END AS status,
       sum(CASE WHEN n_conserved_in_genomes > 0 THEN 1 ELSE 0 END) AS observed
FROM gene_conservation_summary;

SELECT 'ASSERT_PCT_IN_RANGE' AS assertion,
       CASE WHEN count(*) FILTER (WHERE conservation_pct NOT BETWEEN 0 AND 100) = 0
            THEN 'OK' ELSE 'FAIL' END AS status,
       count(*) FILTER (WHERE conservation_pct NOT BETWEEN 0 AND 100) AS observed
FROM gene_conservation_summary;
