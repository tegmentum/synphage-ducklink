-- summary.sql
--
-- Per-gene conservation rollup: for every (query_genome, query_gene),
-- how many of the other genomes carry a best-hit for it, and how strong
-- is the average match?
--
-- Not a direct port of anything in Synphage -- Synphage produces the flat
-- `gene_uniqueness.parquet` and reports at the Polars/Pandas layer. This
-- pushes the rollup into SQL so a caller can `SELECT * FROM
-- gene_conservation_summary WHERE conservation_pct = 100` without leaving
-- DuckDB.
--
-- Inputs:
--   gene_conservation  -- view from gene_conservation.sql
--
-- Output: one row per (query_genome, query_feature_key), with:
--   n_other_genomes         -- how many subject genomes it was checked against
--   n_conserved_in_genomes  -- how many of those had a best-hit
--   conservation_pct        -- 100 * conserved / other
--   avg/min/max_percent_identity -- identity distribution across the hits
--   avg_evalue              -- geometric-ish average of evalues (arith. mean)

CREATE OR REPLACE VIEW gene_conservation_summary AS
SELECT
    query_genome,
    query_feature_key,
    query_gene,
    query_product,
    count(*)                                             AS n_other_genomes,
    count(*) FILTER (WHERE conserved)                    AS n_conserved_in_genomes,
    round(
        100.0 * count(*) FILTER (WHERE conserved)
              / nullif(count(*), 0),
        2
    )                                                    AS conservation_pct,
    avg(percent_identity) FILTER (WHERE conserved)       AS avg_percent_identity,
    min(percent_identity) FILTER (WHERE conserved)       AS min_percent_identity,
    max(percent_identity) FILTER (WHERE conserved)       AS max_percent_identity,
    avg(evalue)           FILTER (WHERE conserved)       AS avg_evalue
FROM gene_conservation
GROUP BY query_genome, query_feature_key, query_gene, query_product
ORDER BY query_genome, query_feature_key;
