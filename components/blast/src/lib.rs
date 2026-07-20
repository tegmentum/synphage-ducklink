//! bio-blast — BLASTN + BLASTP as a DuckLink-loadable component AND a
//! reusable sequence-search capability.
//!
//! Exports two WIT surfaces from one .wasm:
//! - `tegmentum:bio/sequence-search` — the biological capability, callable
//!   by any component-model host that wants pure pairwise alignment.
//! - `duckdb:extension/{guest, callback-dispatch, table-stream-dispatch}` —
//!   DuckLink's dispatch surface. Two table functions, `blastn` and
//!   `blastp`, are registered during `load()` and delegate to the same
//!   underlying `run_search` implementation.

#[allow(warnings)]
mod bindings;

mod align;
mod ducklink;
mod filter;
mod pushdown;
mod scoring;
mod strand;

pub(crate) use bindings::exports::tegmentum::bio::sequence_search::{
    Guest as SequenceSearchGuest, Hit, Scoring, SearchError, SearchOptions, Sequence, Strand,
};

pub(crate) struct Component;

impl SequenceSearchGuest for Component {
    fn search(
        queries: Vec<Sequence>,
        subjects: Vec<Sequence>,
        scoring: Scoring,
        options: SearchOptions,
    ) -> Result<Vec<Hit>, SearchError> {
        run_search(&queries, &subjects, &scoring, &options)
    }
}

/// The one place both the sequence-search Guest impl and the DuckLink
/// bridge come together. Everything above is thin marshalling; everything
/// interesting lives in align.rs / scoring.rs / strand.rs / filter.rs.
pub(crate) fn run_search(
    queries: &[Sequence],
    subjects: &[Sequence],
    scoring: &Scoring,
    options: &SearchOptions,
) -> Result<Vec<Hit>, SearchError> {
    if queries.is_empty() {
        return Err(SearchError::EmptyQueries);
    }
    if subjects.is_empty() {
        return Err(SearchError::EmptySubjects);
    }

    let params = scoring::resolve(scoring).map_err(SearchError::InvalidScoring)?;
    let mut hits: Vec<Hit> = Vec::new();

    for q in queries {
        let query_bytes = q.data.as_bytes();
        if query_bytes.is_empty() {
            return Err(SearchError::InvalidSequence(format!(
                "query '{}' is empty",
                q.key
            )));
        }
        for s in subjects {
            let subject_bytes = s.data.as_bytes();
            if subject_bytes.is_empty() {
                return Err(SearchError::InvalidSequence(format!(
                    "subject '{}' is empty",
                    s.key
                )));
            }

            if let Some(hit) = align_one(
                &q.key,
                &s.key,
                query_bytes,
                subject_bytes,
                &params,
                Strand::Plus,
            ) {
                hits.push(hit);
            }

            if params.both_strands {
                let subject_rc = strand::revcomp(subject_bytes);
                if let Some(hit) = align_one(
                    &q.key,
                    &s.key,
                    query_bytes,
                    &subject_rc,
                    &params,
                    Strand::Minus,
                ) {
                    let subject_len = subject_bytes.len() as u32;
                    hits.push(strand::translate_minus_hit(hit, subject_len));
                }
            }
        }
    }

    Ok(filter::apply(hits, options))
}

/// One (query, subject, orientation) alignment. Returns None if the
/// aligner produces an empty (zero-scoring) local path — those cases
/// carry no useful bit score or coordinates.
fn align_one(
    query_key: &str,
    subject_key: &str,
    query: &[u8],
    subject: &[u8],
    params: &scoring::Params,
    strand: Strand,
) -> Option<Hit> {
    let raw = match &params.kind {
        scoring::Kind::Blastn {
            match_reward,
            mismatch_penalty,
        } => align::run_blastn(
            query,
            subject,
            *match_reward,
            *mismatch_penalty,
            params.gap_open,
            params.gap_extend,
            params.k,
            params.lambda,
        ),
        scoring::Kind::Blastp => align::run_blastp(
            query,
            subject,
            params.gap_open,
            params.gap_extend,
            params.k,
            params.lambda,
        ),
    };

    let raw = raw?;

    // 0-based half-open [start, end) → 1-based inclusive [start, end].
    // start_1 = start_0 + 1, end_1 = end_0. Matches NCBI BLAST output.
    Some(Hit {
        query_key: query_key.to_string(),
        subject_key: subject_key.to_string(),
        query_start: (raw.query_start_0 + 1) as u32,
        query_end: raw.query_end_0 as u32,
        subject_start: (raw.subject_start_0 + 1) as u32,
        subject_end: raw.subject_end_0 as u32,
        strand,
        identity_count: raw.identity_count,
        alignment_length: raw.alignment_length,
        percent_identity: if raw.alignment_length == 0 {
            0.0
        } else {
            (raw.identity_count as f64 / raw.alignment_length as f64) * 100.0
        },
        bit_score: raw.bit_score,
        raw_score: raw.raw_score as f64,
        evalue: raw.e_value,
    })
}

bindings::export!(Component with_types_in bindings);
