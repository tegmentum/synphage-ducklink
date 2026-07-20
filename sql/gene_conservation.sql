-- gene_conservation.sql
--
-- Ports Synphage's `gene_presence.sql`. Synphage's original was a one-line
-- LEFT JOIN of the concatenated BLAST results back to the locus table so
-- every locus shows up whether or not it had a hit:
--
--     SELECT A.*, B.gene
--     FROM locus B LEFT JOIN blastn A USING(key)
--
-- The reason it was that small: Synphage runs BLAST once per pair of
-- genomes and the (query_genome, subject_genome) context is baked into
-- which file each hit came from. Downstream Polars code un-pivots the
-- concatenated frame.
--
-- We do the equivalent join here directly in SQL, and also expand it into
-- the shape the design-doc acceptance test wants: for every (query gene,
-- subject genome) pair the query gene *could* have hit, one row -- with
-- the subject-side annotation and identity if the best hit exists, NULL
-- otherwise. That's "did this gene find a home in every other genome?".
--
-- Inputs:
--   genome_features  -- CDS/gene catalogue for all genomes
--                       (genome_id, feature_key, feature_type,
--                        start_position, end_position, strand,
--                        gene, product).
--   best_hits        -- view from best_hits.sql (which depends on
--                       blast_hits + genome_features).
--
-- Semantic notes:
--   * The cross join enumerates (query_gene, subject_genome) for every
--     subject genome != the query's own genome. If the caller only ran
--     BLAST for some pairs, the missing pairs simply come back with
--     conserved=false -- correct, meaning "we didn't find a hit for this
--     gene in that genome (because we didn't compare, or because BLAST
--     found nothing above the evalue cutoff)".
--   * Self-hits are excluded by `q.genome_id <> sg.subject_genome`. Every
--     gene trivially hits itself; keeping those rows would inflate the
--     conservation rate. Drop the WHERE clause if you want them.
--   * Feature-type filter is CDS by default (matches Synphage). Adjust or
--     remove for RNA / pseudogene / all-feature workflows.

CREATE OR REPLACE VIEW gene_conservation AS
WITH subject_genomes AS (
    SELECT DISTINCT genome_id AS subject_genome
    FROM genome_features
),
query_x_subject AS (
    SELECT
        q.genome_id       AS query_genome,
        q.feature_key     AS query_feature_key,
        q.feature_type    AS query_feature_type,
        q.gene            AS query_gene,
        q.product         AS query_product,
        q.start_position  AS query_start_position,
        q.end_position    AS query_end_position,
        q.strand          AS query_gene_strand,
        sg.subject_genome
    FROM genome_features q
    CROSS JOIN subject_genomes sg
    WHERE q.feature_type = 'CDS'
      AND q.genome_id <> sg.subject_genome
)
SELECT
    x.query_genome,
    x.subject_genome,

    x.query_feature_key,
    x.query_feature_type,
    x.query_gene,
    x.query_product,
    x.query_start_position,
    x.query_end_position,
    x.query_gene_strand,

    b.subject_key           AS subject_feature_key,
    s.gene                  AS subject_gene,
    s.product               AS subject_product,
    s.start_position        AS subject_start_position,
    s.end_position          AS subject_end_position,
    s.strand                AS subject_gene_strand,

    b.query_start           AS hit_query_start,
    b.query_end             AS hit_query_end,
    b.subject_start         AS hit_subject_start,
    b.subject_end           AS hit_subject_end,
    b.strand                AS hit_strand,
    b.alignment_length,
    b.identity_count,
    b.percent_identity,
    b.bit_score,
    b.raw_score,
    b.evalue,

    (b.subject_key IS NOT NULL) AS conserved
FROM query_x_subject x
LEFT JOIN best_hits b
  ON b.query_key       = x.query_feature_key
 AND b.subject_genome  = x.subject_genome
LEFT JOIN genome_features s
  ON s.feature_key = b.subject_key;
