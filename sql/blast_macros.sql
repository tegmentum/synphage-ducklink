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
