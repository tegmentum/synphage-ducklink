LOAD blast;
LOAD jsonfns;   -- gives us to_json() for the sql/blast_macros.sql sugar

CREATE TABLE genome_features (
  genome_id TEXT, feature_key TEXT, feature_type TEXT,
  start_position INTEGER, end_position INTEGER, strand INTEGER,
  gene TEXT, product TEXT
);
INSERT INTO genome_features VALUES
('NC_001416', 'lambdap86', 'CDS', 33299, 33463, -1, 'cIII', 'CIII anti-termination'),
('NC_001416', 'lambdap35', 'CDS', 29118, 29285, -1, '', 'hypothetical protein'),
('NC_001416', 'lambdap67', 'CDS', 42269, 42439, 1, 'NinF', 'protein ninF'),
('NC_001416', 'lambdap65', 'CDS', 41950, 42123, 1, 'NinD', 'NinD protein'),
('NC_001416', 'lambdap39', 'CDS', 31169, 31351, -1, 'orf60a', 'DUF1317 domain-containing protein'),
('NC_001416', 'lambdap66', 'CDS', 42090, 42272, 1, 'NinE', 'NinE family protein'),
('NC_001416', 'lambdap91', 'CDS', 46186, 46368, 1, 'Rz1', 'Rz-like spanin'),
('NC_001416', 'lambdap37', 'CDS', 30839, 31024, -1, 'orf61', 'hypothetical protein'),
('NC_001604', 'T7p60', 'CDS', 39389, 39538, 1, '', 'hypothetical protein'),
('NC_001604', 'T7p02', 'CDS', 1278, 1433, 1, '', 'hypothetical protein'),
('NC_001604', 'T7p11', 'CDS', 7608, 7763, 1, '', 'hypothetical protein'),
('NC_001604', 'T7p34', 'CDS', 17359, 17517, 1, '', 'Gp5.9-like inihibitor of recBCD nuclease'),
('NC_001604', 'T7p06', 'CDS', 1636, 1797, 1, '', 'hypothetical protein'),
('NC_001604', 'T7p59', 'CDS', 38553, 38726, 1, '', 'hypothetical protein'),
('NC_001604', 'T7p16', 'CDS', 8898, 9092, 1, '', 'RNA polymerase inhibitor'),
('NC_001604', 'T7p53', 'CDS', 36344, 36547, 1, '', 'holin'),
('NC_002371', 'P22gp39', 'CDS', 28170, 28328, -1, 'c3', 'CIII anti-termination'),
('NC_002371', 'P22gp18', 'CDS', 14861, 15022, 1, 'arc', 'Arc-like repressor'),
('NC_002371', 'P22gp33', 'CDS', 26449, 26619, -1, 'orf56', 'hypothetical protein'),
('NC_002371', 'P22gp55', 'CDS', 35184, 35357, 1, 'ninD', 'NinD'),
('NC_002371', 'P22gp56', 'CDS', 35324, 35500, 1, 'ninE', 'NinE'),
('NC_002371', 'P22gp58', 'CDS', 35828, 36004, 1, 'ninF', 'NinF'),
('NC_002371', 'P22gp16', 'CDS', 14240, 14419, -1, 'orf59a', 'hypothetical protein'),
('NC_002371', 'P22gp26', 'CDS', 23482, 23661, -1, 'eaG', 'EaG');

-- blast_macros.sql
-- sql/blast_macros.sql
--
-- Thin sugar around the raw `blastn(VARCHAR, VARCHAR, VARCHAR)` and
-- `blastp(VARCHAR, VARCHAR, VARCHAR)` extension calls. The extension takes
-- three JSON string arguments (queries, subjects, options); these macros
-- give the calls names and default-argument slots without changing shape.
--
-- Sample usage:
--
--   LOAD blast;
--   .read sql/blast_macros.sql;
--
--   SELECT * FROM blastn_of(
--       '[{"key":"q1","data":"ACGTACGTACGTACGTACGT"}]',
--       '[{"key":"s1","data":"GGGACGTACGTACGTACGTACGTGGG"}]',
--       '{"evalue_max": 1e-5}'
--   );
--
-- Building the JSON payloads. If DuckDB's built-in `json` extension is
-- available (`LOAD json`), the tightest sugar is:
--
--   WITH q AS (SELECT feature_key AS key, sequence AS data FROM genes)
--   SELECT * FROM blastn_of(
--       (SELECT to_json(list(struct_pack(key := key, data := data))) FROM q),
--       (SELECT to_json(list(struct_pack(key := key, data := data))) FROM q),
--       to_json({'evalue_max': 1e-5})
--   );
--
-- (Note: `to_json` on arbitrary values is DuckDB core's `json` extension,
--  not DuckLink's `jsonfns` extension whose `to_json` is a VARCHAR->VARCHAR
--  validator only.) When the json extension is not loadable in a given
-- DuckLink build -- e.g. today's compile-time flag disables external
-- extension loading -- callers construct the JSON string themselves and
-- pass it in. Every extension consumer of `blast` sees the same three-
-- argument surface either way.
--
-- Why not a fuller `blastn(TABLE genes)` sugar? DuckDB's binder currently
-- rejects scalar subqueries and LATERAL column parameters inside table
-- function argument positions, so `blastn((SELECT ...))` fails at bind
-- time regardless of how the extension is registered. The moment that
-- lifts, the natural sugar drops in without changing the extension.

CREATE OR REPLACE MACRO blastn_of(queries_json, subjects_json, opts_json := NULL) AS TABLE
    SELECT * FROM blastn(queries_json, subjects_json, opts_json);

CREATE OR REPLACE MACRO blastp_of(queries_json, subjects_json, opts_json := NULL) AS TABLE
    SELECT * FROM blastp(queries_json, subjects_json, opts_json);


CREATE TABLE blast_hits AS
SELECT * FROM blastn_of('[{"key": "lambdap86", "data": "ATGCAATATGCCATTGCAGGGTGGCCTGTTGCTGGCTGCCCTTCCGAATCTTTACTTGAACGAATCACCCGTAAATTACGTGACGGATGGAAACGCCTTATCGACATACTTAATCAGCCAGGAGTCCCAAAGAATGGATCAAACACTTATGGCTATCCAGACTAA"}, {"key": "lambdap35", "data": "ATGCACTTCCGAGTCACAGGAGAATGGAATGGAGAGCCATTCAACAGAGTTATCGAAGCGGAGAACATCAACGACTGCTACGACCACTGGATGATATGGGCGCAGATAGCACATGCAGACGTAACCAATATTCGAATTGAAGAACTGAAAGAACACCAAGCCGCCTGA"}, {"key": "lambdap67", "data": "GTGATTGACCAAAATCGAAGTTACGAACAAGAAAGCGTCGAGCGAGCTTTAACGTGCGCTAACTGCGGTCAGAAGCTGCATGTGCTGGAAGTTCACGTGTGTGAGCACTGCTGCGCAGAACTGATGAGCGATCCGAATAGCTCGATGCACGAGGAAGAAGATGATGGCTAA"}, {"key": "lambdap65", "data": "ATGATGCGATGTTATCGGTGCGGTGAATGCAAAGAAGATAACCGCTTCCGACCAAATCAACCTTACTGGAATCGATGGTGTCTCCGGTGTGAAAGAACACCAACAGGGGTGTTACCACTACCGCAGGAAAAGGAGGACGTGTGGCGAGACAGCGACGAAGTATCACCGACATAA"}, {"key": "lambdap39", "data": "ATGACGCATCCTCACGATAATATCCGGGTAGGCGCAATCACTTTCGTCTACTCCGTTACAAAGCGAGGCTGGGTATTTCCCGGCCTTTCTGTTATCCGAAATCCACTGAAAGCACAGCGGCTGGCTGAGGAGATAAATAATAAACGAGGGGCTGTATGCACAAAGCATCTTCTGTTGAGTTAA"}, {"key": "lambdap66", "data": "GTGGCGAGACAGCGACGAAGTATCACCGACATAATCTGCGAAAACTGCAAATACCTTCCAACGAAACGCACCAGAAATAAACCCAAGCCAATCCCAAAAGAATCTGACGTAAAAACCTTCAACTACACGGCTCACCTGTGGGATATCCGGTGGCTAAGACGTCGTGCGAGGAAAACAAGGTGA"}, {"key": "lambdap91", "data": "ATGCTAAAGCTGAAAATGATGCTCTGCGTGATGATGTTGCCGCTGGTCGTCGTCGGTTGCACATCAAAGCAGTCTGTCAGTCAGTGCGTGAAGCCACCACCGCCTCCGGCGTGGATAATGCAGCCTCCCCCCGACTGGCAGACACCGCTGAACGGGATTATTTCACCCTCAGAGAGAGGCTGA"}, {"key": "lambdap37", "data": "ATGAGAGAAACCAGGTATGACAACCACGGAATGCATTTTTCTGGCAGCGGGCTTCATATTCTGTGTGCTTATGCTTGCCGACATGGGACTTGTTCAATGACACCTCAGCAGGAAAACGCCCTTCGCAGCATTGCCCGTCAGGCTAATTCTGAAATCAAAAAAAGCCAGACAGCAGTTTCCGGATAA"}, {"key": "T7p60", "data": "ATGTTCCGCTTATTGTTGAACCTACTGCGGCATAGAGTCACCTACCGATTTCTTGTGGTACTTTGTGCTGCCCTTGGGTACGCATCTCTTACTGGAGACCTCAGTTCACTGGAGTCTGTCGTTTGCTCTATACTCACTTGTAGCGATTAG"}, {"key": "T7p02", "data": "ATGTCTACTACCAACGTGCAATACGGTCTGACCGCTCAAACTGTACTTTTCTATAGCGACATGGTGCGCTGTGGCTTTAACTGGTCACTCGCAATGGCACAGCTCAAAGAACTGTACGAAAACAACAAGGCAATAGCTTTAGAATCTGCTGAGTGA"}, {"key": "T7p11", "data": "ATGTTTAAGAAGGTTGGTAAATTCCTTGCGGCTTTGGCAGCTATCCTGACGCTTGCGTATATTCTTGCGGTATACCCTCAAGTAGCACTAGTAGTAGTTGGCGCTTGTTACTTAGCGGCAGTGTGTGCTTGCGTGTGGAGTATAGTTAACTGGTAA"}, {"key": "T7p34", "data": "ATGTCTCGTGACCTTGTGACTATTCCACGCGATGTGTGGAACGATATACAGGGCTACATCGACTCTCTGGAACGTGAGAACGATAGCCTTAAGAATCAACTAATGGAAGCTGACGAATACGTAGCGGAACTAGAGGAGAAACTTAATGGCACTTCTTGA"}, {"key": "T7p06", "data": "ATGATGAAGCACTACGTTATGCCAATCCACACGTCCAACGGGGCAACCGTATGTACACCTGATGGGTTCGCAATGAAACAACGAATCGAACGCCTTAAGCGTGAACTCCGCATTAACCGCAAGATTAACAAGATAGGTTCCGGCTATGACAGAACGCACTGA"}, {"key": "T7p59", "data": "ATGGCTACTCCGATAAGACCCTTGAGTTACTCGCTAAGAAGGCAAAGCAATGGGGAGTCCAGACGGTTGTCTACGAGAGTAACTTCGGTGACGGTATGTTCGGTAAGGTATTCAGTCCTATCCTTCTTAAACACCACAACTGTGCGATGGAAGAGATTCGTGCCCGTGGTATGA"}, {"key": "T7p16", "data": "ATGTCAAACGTAAATACAGGTTCACTTAGTGTGGACAATAAGAAGTTTTGGGCTACCGTAGAGTCCTCGGAGCATTCCTTCGAGGTTCCAATCTACGCTGAGACCCTAGACGAAGCTCTGGAGTTAGCCGAATGGCAATACGTTCCGGCTGGCTTTGAGGTTACTCGTGTGCGTCCTTGTGTAGCACCGAAGTAA"}, {"key": "T7p53", "data": "GTGCTATCATTAGACTTTAACAACGAATTGATTAAGGCTGCTCCAATTGTTGGGACGGGTGTAGCAGATGTTAGTGCTCGACTGTTCTTTGGGTTAAGCCTTAACGAATGGTTCTACGTTGCTGCTATCGCCTACACAGTGGTTCAGATTGGTGCCAAGGTAGTCGATAAGATGATTGACTGGAAGAAAGCCAATAAGGAGTGA"}, {"key": "P22gp39", "data": "ATGATTATCGCAATCGCGGGAAGCGCTCGCATGGGCGTTTCCCAGTTACACGAATCACTTTTAGATCGCATCACCCGCAAATTACGCGCTGGCTGGAAACGGCTGGCAGACATCCTTAATCAGCCTGGAGTGCCGAGCCATGACTATTGTGCCTGTTAA"}, {"key": "P22gp18", "data": "ATGAAAGGAATGAGCAAAATGCCGCAGTTCAATTTGCGGTGGCCTAGAGAAGTATTGGATTTGGTACGCAAGGTAGCGGAAGAGAATGGTCGGTCTGTTAATTCTGAGATTTATCAGCGAGTAATGGAAAGCTTTAAGAAGGAAGGGCGCATTGGCGCGTAA"}, {"key": "P22gp33", "data": "ATGAGAGGACTTGCATACAATCCCGACATTCTACCAGCCGAACTTATTATTAGGCACAAAATTAAACCAATGCCTACACGCGAAGAATTATTGCAGCGCAATTCATTTCCTTCTATTAACGAGAATAAATATTTGAATGCGATACTGAGGAAAGATAAATGCAACAGGTAA"}, {"key": "P22gp55", "data": "ATGAAACACTGCTACCGCTGCGGAGAAAGCAAAGACGATTATCGATTCCGGCCAAATCAACCTTATTGGCACCAATGGTGTATCAGATGTGAGCGGTCGCCAGTAGGTAATTTCCCGCTGCCAGAGACGAAGGAGGACGTATGGCACGACAGCGACGAAGTATCACCGACATAA"}, {"key": "P22gp56", "data": "ATGGCACGACAGCGACGAAGTATCACCGACATAATCTGCGAAAACTGCAAATACCTTCCAACGAAACGTTCCAGAAATAAACGCAAGCCAATCCCAAAAGAGTCTGACGTAAAAACCTTCAATTACACGGCTCACCTGTGGGATATCCGGTGGCTAAGACATCGTGCGAGGAAATGA"}, {"key": "P22gp58", "data": "ATGCTTAGCCCATCACAATCCCTTCAATACCAGAAAGAAAGCGTCGAGCGGGCTTTAACGTGCGCTAACTGCGGTCAGAAGCTGCATGTGCTGGAAGTTCATGTATGTGAAGCGTGCTGCGCAGAACTGATGAGCGATCCGAATAGCTCAATGTACGAGGAAGAAGACGATGGCTAA"}, {"key": "P22gp16", "data": "ATGGCGAAAAAACCAGGTGAAAACACAGGAAAAAACGGCGGAATATACCAAGAAGTTGGCCCGCGAGGCGGTAAGAAAGACAATTTTGCCACCGTCAAGGACAACGAAAGGCTTCCACCAACAACAAAGCCAGGTAATGGCTGGGTATTAGATAAGCGAACTCCAGACAGCAAAAAGTAA"}, {"key": "P22gp26", "data": "ATGTCTTGTCCAAAATGCGGTTCTGGAAATATTGCAAAAGAAAAAACAATGCGTGGATGGTCTGATGATTATGTGTGCTGCGATTGCGGATACAACGACTCTAAAGACGCATTTGGAGAGCGTGGTAAAAACGAGTTTGTCAGAATTAATAAGGAACGCAAAGGCAACGAAAAAAGCTAA"}]', '[{"key": "lambdap86", "data": "ATGCAATATGCCATTGCAGGGTGGCCTGTTGCTGGCTGCCCTTCCGAATCTTTACTTGAACGAATCACCCGTAAATTACGTGACGGATGGAAACGCCTTATCGACATACTTAATCAGCCAGGAGTCCCAAAGAATGGATCAAACACTTATGGCTATCCAGACTAA"}, {"key": "lambdap35", "data": "ATGCACTTCCGAGTCACAGGAGAATGGAATGGAGAGCCATTCAACAGAGTTATCGAAGCGGAGAACATCAACGACTGCTACGACCACTGGATGATATGGGCGCAGATAGCACATGCAGACGTAACCAATATTCGAATTGAAGAACTGAAAGAACACCAAGCCGCCTGA"}, {"key": "lambdap67", "data": "GTGATTGACCAAAATCGAAGTTACGAACAAGAAAGCGTCGAGCGAGCTTTAACGTGCGCTAACTGCGGTCAGAAGCTGCATGTGCTGGAAGTTCACGTGTGTGAGCACTGCTGCGCAGAACTGATGAGCGATCCGAATAGCTCGATGCACGAGGAAGAAGATGATGGCTAA"}, {"key": "lambdap65", "data": "ATGATGCGATGTTATCGGTGCGGTGAATGCAAAGAAGATAACCGCTTCCGACCAAATCAACCTTACTGGAATCGATGGTGTCTCCGGTGTGAAAGAACACCAACAGGGGTGTTACCACTACCGCAGGAAAAGGAGGACGTGTGGCGAGACAGCGACGAAGTATCACCGACATAA"}, {"key": "lambdap39", "data": "ATGACGCATCCTCACGATAATATCCGGGTAGGCGCAATCACTTTCGTCTACTCCGTTACAAAGCGAGGCTGGGTATTTCCCGGCCTTTCTGTTATCCGAAATCCACTGAAAGCACAGCGGCTGGCTGAGGAGATAAATAATAAACGAGGGGCTGTATGCACAAAGCATCTTCTGTTGAGTTAA"}, {"key": "lambdap66", "data": "GTGGCGAGACAGCGACGAAGTATCACCGACATAATCTGCGAAAACTGCAAATACCTTCCAACGAAACGCACCAGAAATAAACCCAAGCCAATCCCAAAAGAATCTGACGTAAAAACCTTCAACTACACGGCTCACCTGTGGGATATCCGGTGGCTAAGACGTCGTGCGAGGAAAACAAGGTGA"}, {"key": "lambdap91", "data": "ATGCTAAAGCTGAAAATGATGCTCTGCGTGATGATGTTGCCGCTGGTCGTCGTCGGTTGCACATCAAAGCAGTCTGTCAGTCAGTGCGTGAAGCCACCACCGCCTCCGGCGTGGATAATGCAGCCTCCCCCCGACTGGCAGACACCGCTGAACGGGATTATTTCACCCTCAGAGAGAGGCTGA"}, {"key": "lambdap37", "data": "ATGAGAGAAACCAGGTATGACAACCACGGAATGCATTTTTCTGGCAGCGGGCTTCATATTCTGTGTGCTTATGCTTGCCGACATGGGACTTGTTCAATGACACCTCAGCAGGAAAACGCCCTTCGCAGCATTGCCCGTCAGGCTAATTCTGAAATCAAAAAAAGCCAGACAGCAGTTTCCGGATAA"}, {"key": "T7p60", "data": "ATGTTCCGCTTATTGTTGAACCTACTGCGGCATAGAGTCACCTACCGATTTCTTGTGGTACTTTGTGCTGCCCTTGGGTACGCATCTCTTACTGGAGACCTCAGTTCACTGGAGTCTGTCGTTTGCTCTATACTCACTTGTAGCGATTAG"}, {"key": "T7p02", "data": "ATGTCTACTACCAACGTGCAATACGGTCTGACCGCTCAAACTGTACTTTTCTATAGCGACATGGTGCGCTGTGGCTTTAACTGGTCACTCGCAATGGCACAGCTCAAAGAACTGTACGAAAACAACAAGGCAATAGCTTTAGAATCTGCTGAGTGA"}, {"key": "T7p11", "data": "ATGTTTAAGAAGGTTGGTAAATTCCTTGCGGCTTTGGCAGCTATCCTGACGCTTGCGTATATTCTTGCGGTATACCCTCAAGTAGCACTAGTAGTAGTTGGCGCTTGTTACTTAGCGGCAGTGTGTGCTTGCGTGTGGAGTATAGTTAACTGGTAA"}, {"key": "T7p34", "data": "ATGTCTCGTGACCTTGTGACTATTCCACGCGATGTGTGGAACGATATACAGGGCTACATCGACTCTCTGGAACGTGAGAACGATAGCCTTAAGAATCAACTAATGGAAGCTGACGAATACGTAGCGGAACTAGAGGAGAAACTTAATGGCACTTCTTGA"}, {"key": "T7p06", "data": "ATGATGAAGCACTACGTTATGCCAATCCACACGTCCAACGGGGCAACCGTATGTACACCTGATGGGTTCGCAATGAAACAACGAATCGAACGCCTTAAGCGTGAACTCCGCATTAACCGCAAGATTAACAAGATAGGTTCCGGCTATGACAGAACGCACTGA"}, {"key": "T7p59", "data": "ATGGCTACTCCGATAAGACCCTTGAGTTACTCGCTAAGAAGGCAAAGCAATGGGGAGTCCAGACGGTTGTCTACGAGAGTAACTTCGGTGACGGTATGTTCGGTAAGGTATTCAGTCCTATCCTTCTTAAACACCACAACTGTGCGATGGAAGAGATTCGTGCCCGTGGTATGA"}, {"key": "T7p16", "data": "ATGTCAAACGTAAATACAGGTTCACTTAGTGTGGACAATAAGAAGTTTTGGGCTACCGTAGAGTCCTCGGAGCATTCCTTCGAGGTTCCAATCTACGCTGAGACCCTAGACGAAGCTCTGGAGTTAGCCGAATGGCAATACGTTCCGGCTGGCTTTGAGGTTACTCGTGTGCGTCCTTGTGTAGCACCGAAGTAA"}, {"key": "T7p53", "data": "GTGCTATCATTAGACTTTAACAACGAATTGATTAAGGCTGCTCCAATTGTTGGGACGGGTGTAGCAGATGTTAGTGCTCGACTGTTCTTTGGGTTAAGCCTTAACGAATGGTTCTACGTTGCTGCTATCGCCTACACAGTGGTTCAGATTGGTGCCAAGGTAGTCGATAAGATGATTGACTGGAAGAAAGCCAATAAGGAGTGA"}, {"key": "P22gp39", "data": "ATGATTATCGCAATCGCGGGAAGCGCTCGCATGGGCGTTTCCCAGTTACACGAATCACTTTTAGATCGCATCACCCGCAAATTACGCGCTGGCTGGAAACGGCTGGCAGACATCCTTAATCAGCCTGGAGTGCCGAGCCATGACTATTGTGCCTGTTAA"}, {"key": "P22gp18", "data": "ATGAAAGGAATGAGCAAAATGCCGCAGTTCAATTTGCGGTGGCCTAGAGAAGTATTGGATTTGGTACGCAAGGTAGCGGAAGAGAATGGTCGGTCTGTTAATTCTGAGATTTATCAGCGAGTAATGGAAAGCTTTAAGAAGGAAGGGCGCATTGGCGCGTAA"}, {"key": "P22gp33", "data": "ATGAGAGGACTTGCATACAATCCCGACATTCTACCAGCCGAACTTATTATTAGGCACAAAATTAAACCAATGCCTACACGCGAAGAATTATTGCAGCGCAATTCATTTCCTTCTATTAACGAGAATAAATATTTGAATGCGATACTGAGGAAAGATAAATGCAACAGGTAA"}, {"key": "P22gp55", "data": "ATGAAACACTGCTACCGCTGCGGAGAAAGCAAAGACGATTATCGATTCCGGCCAAATCAACCTTATTGGCACCAATGGTGTATCAGATGTGAGCGGTCGCCAGTAGGTAATTTCCCGCTGCCAGAGACGAAGGAGGACGTATGGCACGACAGCGACGAAGTATCACCGACATAA"}, {"key": "P22gp56", "data": "ATGGCACGACAGCGACGAAGTATCACCGACATAATCTGCGAAAACTGCAAATACCTTCCAACGAAACGTTCCAGAAATAAACGCAAGCCAATCCCAAAAGAGTCTGACGTAAAAACCTTCAATTACACGGCTCACCTGTGGGATATCCGGTGGCTAAGACATCGTGCGAGGAAATGA"}, {"key": "P22gp58", "data": "ATGCTTAGCCCATCACAATCCCTTCAATACCAGAAAGAAAGCGTCGAGCGGGCTTTAACGTGCGCTAACTGCGGTCAGAAGCTGCATGTGCTGGAAGTTCATGTATGTGAAGCGTGCTGCGCAGAACTGATGAGCGATCCGAATAGCTCAATGTACGAGGAAGAAGACGATGGCTAA"}, {"key": "P22gp16", "data": "ATGGCGAAAAAACCAGGTGAAAACACAGGAAAAAACGGCGGAATATACCAAGAAGTTGGCCCGCGAGGCGGTAAGAAAGACAATTTTGCCACCGTCAAGGACAACGAAAGGCTTCCACCAACAACAAAGCCAGGTAATGGCTGGGTATTAGATAAGCGAACTCCAGACAGCAAAAAGTAA"}, {"key": "P22gp26", "data": "ATGTCTTGTCCAAAATGCGGTTCTGGAAATATTGCAAAAGAAAAAACAATGCGTGGATGGTCTGATGATTATGTGTGCTGCGATTGCGGATACAACGACTCTAAAGACGCATTTGGAGAGCGTGGTAAAAACGAGTTTGTCAGAATTAATAAGGAACGCAAAGGCAACGAAAAAAGCTAA"}]', '{"evalue_max": 1e-05}');

-- best_hits.sql
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


-- gene_conservation.sql
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


-- summary.sql
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

