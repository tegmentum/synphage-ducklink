-- best_hits.sql
--
-- Ports Synphage's implicit "best hit" convention into an explicit ranking.
--
-- Synphage's `parse_blastn.sql` picks `hits[0]` and `hsps[0]` from BLAST's
-- JSON output for each query, relying on NCBI BLAST's own ordering. Our
-- `bio:blast` component emits every hit surviving the caller's options,
-- so we do the ranking ourselves in SQL. The design doc's acceptance-test
-- CTE names the same ordering; this file is that CTE promoted to a view.
--
-- Best hit per (query_key, subject_genome), ordered by:
--   bit_score DESC, evalue ASC, subject_key ASC (deterministic tie-break).
--
-- Inputs:
--   blast_hits       -- as emitted by the bio:blast wasm component
--                       (query_key, subject_key, query_start, query_end,
--                        subject_start, subject_end, strand,
--                        identity_count, alignment_length,
--                        percent_identity, bit_score, raw_score, evalue).
--   genome_features  -- annotation table; used only to look up which genome
--                       a subject_key belongs to. Any table exposing
--                       (feature_key, genome_id) works.
--
-- Output columns: every column from blast_hits plus subject_genome.

CREATE OR REPLACE VIEW best_hits AS
WITH ranked AS (
    SELECT
        h.query_key,
        h.subject_key,
        f.genome_id AS subject_genome,
        h.query_start,
        h.query_end,
        h.subject_start,
        h.subject_end,
        h.strand,
        h.identity_count,
        h.alignment_length,
        h.percent_identity,
        h.bit_score,
        h.raw_score,
        h.evalue,
        row_number() OVER (
            PARTITION BY h.query_key, f.genome_id
            ORDER BY h.bit_score DESC,
                     h.evalue ASC,
                     h.subject_key ASC
        ) AS rk
    FROM blast_hits h
    LEFT JOIN genome_features f
      ON f.feature_key = h.subject_key
)
SELECT
    query_key,
    subject_key,
    subject_genome,
    query_start,
    query_end,
    subject_start,
    subject_end,
    strand,
    identity_count,
    alignment_length,
    percent_identity,
    bit_score,
    raw_score,
    evalue
FROM ranked
WHERE rk = 1;
