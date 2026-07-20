//! Filter-pushdown planner for the DuckLink `blast{n,p}` table functions.
//!
//! DuckDB pushes a conjunctive set of predicates below the scan as a
//! `list<table-filter>`; each filter is `(column, op, values)` where
//! `column` indexes our registered `hit_columns()`:
//!
//! ```text
//! 0: query_key (Text)         7:  identity_count (Uint32)
//! 1: subject_key (Text)       8:  alignment_length (Uint32)
//! 2: query_start (Uint32)     9:  percent_identity (Float64)
//! 3: query_end (Uint32)       10: bit_score (Float64)
//! 4: subject_start (Uint32)   11: raw_score (Float64)
//! 5: subject_end (Uint32)     12: evalue (Float64)
//! 6: strand (Text)
//! ```
//!
//! We recognise five profitable patterns and let the DuckLink core re-check
//! everything else above the scan (per the freeze policy, dropping a
//! pushed filter only forgoes optimisation — never correctness):
//!
//! 1. `evalue < X` / `evalue <= X`             → tighten `SearchOptions.evalue_max`.
//! 2. `percent_identity >= X` / `... > X`      → tighten `SearchOptions.min_identity`.
//! 3. `query_key = K` / `query_key IN (...)`   → prune the queries batch pre-search.
//! 4. `subject_key = K` / `subject_key IN ...` → prune the subjects batch pre-search.
//! 5. `strand = 'plus' | 'minus'`              → drop non-matching hits post-search.
//!
//! Anything else is ignored (the split above keeps ducklink.rs tidy and
//! the recognition logic unit-testable without any WASM host).
//!
//! Multiple clauses on the same column combine as a *conjunction*:
//! evalue ceilings meet by MIN, identity floors by MAX, key sets by
//! set-intersection, strand by equality (contradictory clauses trip the
//! `short_circuit` flag and the scan opens on an empty cursor).

use std::collections::HashSet;

use crate::bindings::duckdb::extension::types::Duckvalue;
use crate::bindings::exports::duckdb::extension::table_stream_dispatch::{FilterOp, TableFilter};
use crate::bindings::exports::tegmentum::bio::sequence_search::{SearchOptions, Sequence, Strand};

// Column indices — mirror `ducklink::hit_columns()`. Kept as a `const`
// block so the schema and the pushdown planner drift together or not at
// all.
const COL_QUERY_KEY: u32 = 0;
const COL_SUBJECT_KEY: u32 = 1;
const COL_STRAND: u32 = 6;
const COL_PERCENT_IDENTITY: u32 = 9;
const COL_EVALUE: u32 = 12;

/// The strand a `strand = 'plus'|'minus'` clause pinned the scan to. Kept
/// separate from the WIT `Strand` so tests don't have to spin up a
/// `Sequence` just to talk about strand filters.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StrandFilter {
    Plus,
    Minus,
}

impl StrandFilter {
    pub fn matches(self, s: Strand) -> bool {
        matches!(
            (self, s),
            (StrandFilter::Plus, Strand::Plus) | (StrandFilter::Minus, Strand::Minus)
        )
    }
}

/// The tightened / restricted view of the scan after every recognisable
/// pushed filter has been folded in. Callers apply it: rewrite
/// `SearchOptions`, `retain` the query/subject batches, run the search,
/// then post-filter hits by `strand_keep`.
///
/// If `short_circuit` is set the whole plan is infeasible — the scan
/// yields zero rows and `run_search` can be skipped entirely.
#[derive(Debug, Default, Clone)]
pub struct PushdownPlan {
    pub evalue_ceiling: Option<f64>,
    pub identity_floor: Option<f64>,
    /// `None` = every query key allowed. `Some(set)` = only keys in `set`.
    pub query_keys: Option<HashSet<String>>,
    /// `None` = every subject key allowed. `Some(set)` = only keys in `set`.
    pub subject_keys: Option<HashSet<String>>,
    pub strand_keep: Option<StrandFilter>,
    pub short_circuit: bool,
}

/// Classify one filter and merge it into `plan`. Unrecognised patterns
/// (columns we don't know how to prune on, unshippable value types, ops
/// we don't handle) are silently dropped — the DuckLink core re-checks
/// them above the scan.
fn absorb(plan: &mut PushdownPlan, filter: &TableFilter) {
    match filter.column {
        COL_EVALUE => absorb_evalue(plan, filter),
        COL_PERCENT_IDENTITY => absorb_percent_identity(plan, filter),
        COL_QUERY_KEY => absorb_text_key(&mut plan.query_keys, &mut plan.short_circuit, filter),
        COL_SUBJECT_KEY => absorb_text_key(&mut plan.subject_keys, &mut plan.short_circuit, filter),
        COL_STRAND => absorb_strand(plan, filter),
        _ => {} // ignore — core re-checks above the scan
    }
}

fn absorb_evalue(plan: &mut PushdownPlan, filter: &TableFilter) {
    if !matches!(filter.op, FilterOp::Lt | FilterOp::Le) {
        return;
    }
    let Some(x) = as_f64(filter.values.first()) else {
        return;
    };
    plan.evalue_ceiling = Some(match plan.evalue_ceiling {
        Some(cur) => cur.min(x),
        None => x,
    });
}

fn absorb_percent_identity(plan: &mut PushdownPlan, filter: &TableFilter) {
    if !matches!(filter.op, FilterOp::Ge | FilterOp::Gt) {
        return;
    }
    let Some(x) = as_f64(filter.values.first()) else {
        return;
    };
    plan.identity_floor = Some(match plan.identity_floor {
        Some(cur) => cur.max(x),
        None => x,
    });
}

fn absorb_text_key(
    slot: &mut Option<HashSet<String>>,
    short_circuit: &mut bool,
    filter: &TableFilter,
) {
    let incoming: HashSet<String> = match filter.op {
        FilterOp::Eq => match as_text(filter.values.first()) {
            Some(s) => std::iter::once(s).collect(),
            None => return, // unshippable constant — skip
        },
        FilterOp::IsIn => {
            let set: HashSet<String> = filter.values.iter().filter_map(as_text_ref).collect();
            if set.is_empty() {
                // IN () can't match anything.
                *short_circuit = true;
                return;
            }
            set
        }
        _ => return, // ne / lt / gt / null-checks on text keys: ignore
    };

    *slot = Some(match slot.take() {
        Some(existing) => {
            let merged: HashSet<String> = existing.intersection(&incoming).cloned().collect();
            if merged.is_empty() {
                *short_circuit = true;
            }
            merged
        }
        None => incoming,
    });
}

fn absorb_strand(plan: &mut PushdownPlan, filter: &TableFilter) {
    if !matches!(filter.op, FilterOp::Eq) {
        return;
    }
    let Some(s) = as_text(filter.values.first()) else {
        return;
    };
    let incoming = match s.as_str() {
        "plus" => StrandFilter::Plus,
        "minus" => StrandFilter::Minus,
        _ => return, // unknown label — leave for core to reject above the scan
    };
    match plan.strand_keep {
        None => plan.strand_keep = Some(incoming),
        Some(cur) if cur == incoming => {}
        Some(_) => plan.short_circuit = true, // plus AND minus is impossible
    }
}

/// Build a `PushdownPlan` from the conjunctive filter set the host pushed
/// into `call-table-open-filtered`.
pub fn plan(filters: &[TableFilter]) -> PushdownPlan {
    let mut p = PushdownPlan::default();
    for f in filters {
        absorb(&mut p, f);
    }
    p
}

/// Fold the plan's numeric bounds into `opts`. Tightening only:
/// existing user options are never loosened.
pub fn tighten_options(plan: &PushdownPlan, opts: &mut SearchOptions) {
    if let Some(e) = plan.evalue_ceiling {
        opts.evalue_max = Some(match opts.evalue_max {
            Some(cur) => cur.min(e),
            None => e,
        });
    }
    if let Some(i) = plan.identity_floor {
        opts.min_identity = Some(match opts.min_identity {
            Some(cur) => cur.max(i),
            None => i,
        });
    }
}

/// Drop queries/subjects the key-restriction clauses ruled out. Returns
/// `true` if either batch is now empty (nothing left to search).
pub fn prune_batches(
    plan: &PushdownPlan,
    queries: &mut Vec<Sequence>,
    subjects: &mut Vec<Sequence>,
) -> bool {
    if let Some(ks) = &plan.query_keys {
        queries.retain(|q| ks.contains(&q.key));
    }
    if let Some(ks) = &plan.subject_keys {
        subjects.retain(|s| ks.contains(&s.key));
    }
    queries.is_empty() || subjects.is_empty()
}

// ---- value coercion helpers ---------------------------------------------

/// Coerce a `Duckvalue` numeric literal to `f64`. Returns `None` for
/// anything the DuckDB planner shouldn't have pushed at us for a
/// float column (Text, Null, Blob, ...).
fn as_f64(v: Option<&Duckvalue>) -> Option<f64> {
    match v? {
        Duckvalue::Float64(x) => Some(*x),
        Duckvalue::Float32(x) => Some(*x as f64),
        Duckvalue::Int64(x) => Some(*x as f64),
        Duckvalue::Int32(x) => Some(*x as f64),
        Duckvalue::Int16(x) => Some(*x as f64),
        Duckvalue::Int8(x) => Some(*x as f64),
        Duckvalue::Uint64(x) => Some(*x as f64),
        Duckvalue::Uint32(x) => Some(*x as f64),
        Duckvalue::Uint16(x) => Some(*x as f64),
        Duckvalue::Uint8(x) => Some(*x as f64),
        _ => None,
    }
}

fn as_text(v: Option<&Duckvalue>) -> Option<String> {
    match v? {
        Duckvalue::Text(s) => Some(s.clone()),
        _ => None,
    }
}

fn as_text_ref(v: &Duckvalue) -> Option<String> {
    match v {
        Duckvalue::Text(s) => Some(s.clone()),
        _ => None,
    }
}

// ---- tests --------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn f(column: u32, op: FilterOp, values: Vec<Duckvalue>) -> TableFilter {
        TableFilter {
            column,
            op,
            values,
        }
    }

    fn text(s: &str) -> Duckvalue {
        Duckvalue::Text(s.to_string())
    }

    fn seq(k: &str) -> Sequence {
        Sequence {
            key: k.to_string(),
            data: "ACGT".to_string(),
        }
    }

    #[test]
    fn evalue_lt_and_le_tighten_ceiling_by_min() {
        // Two conjunctive evalue ceilings — the tighter one wins.
        let p = plan(&[
            f(COL_EVALUE, FilterOp::Lt, vec![Duckvalue::Float64(1e-3)]),
            f(COL_EVALUE, FilterOp::Le, vec![Duckvalue::Float64(1e-5)]),
        ]);
        assert_eq!(p.evalue_ceiling, Some(1e-5));
        assert!(!p.short_circuit);

        // And tightening never loosens an existing user option.
        let mut opts = SearchOptions {
            evalue_max: Some(1e-8),
            max_target_seqs: None,
            min_identity: None,
        };
        tighten_options(&p, &mut opts);
        assert_eq!(opts.evalue_max, Some(1e-8));

        // With no prior option, the plan's ceiling is what sticks.
        let mut opts2 = SearchOptions {
            evalue_max: None,
            max_target_seqs: None,
            min_identity: None,
        };
        tighten_options(&p, &mut opts2);
        assert_eq!(opts2.evalue_max, Some(1e-5));
    }

    #[test]
    fn percent_identity_ge_gt_tighten_floor_by_max() {
        let p = plan(&[
            f(COL_PERCENT_IDENTITY, FilterOp::Ge, vec![Duckvalue::Float64(70.0)]),
            f(COL_PERCENT_IDENTITY, FilterOp::Gt, vec![Duckvalue::Float64(90.0)]),
        ]);
        assert_eq!(p.identity_floor, Some(90.0));

        // Never loosens an existing higher floor.
        let mut opts = SearchOptions {
            evalue_max: None,
            max_target_seqs: None,
            min_identity: Some(95.0),
        };
        tighten_options(&p, &mut opts);
        assert_eq!(opts.min_identity, Some(95.0));
    }

    #[test]
    fn key_eq_and_is_in_prune_batches() {
        let p = plan(&[
            f(COL_QUERY_KEY, FilterOp::Eq, vec![text("qA")]),
            f(COL_SUBJECT_KEY, FilterOp::IsIn, vec![text("s1"), text("s3")]),
        ]);
        assert_eq!(
            p.query_keys.as_ref().map(|s| s.iter().cloned().collect::<Vec<_>>()),
            Some(vec!["qA".to_string()])
        );
        assert_eq!(p.subject_keys.as_ref().map(|s| s.len()), Some(2));

        let mut qs = vec![seq("qA"), seq("qB"), seq("qC")];
        let mut ss = vec![seq("s1"), seq("s2"), seq("s3")];
        let empty = prune_batches(&p, &mut qs, &mut ss);
        assert!(!empty);
        assert_eq!(qs.iter().map(|q| q.key.as_str()).collect::<Vec<_>>(), vec!["qA"]);
        assert_eq!(
            ss.iter().map(|s| s.key.as_str()).collect::<Vec<_>>(),
            vec!["s1", "s3"]
        );
    }

    #[test]
    fn contradictory_key_clauses_short_circuit() {
        // query_key = 'qA' AND query_key IN ('qB','qC') is unsatisfiable.
        let p = plan(&[
            f(COL_QUERY_KEY, FilterOp::Eq, vec![text("qA")]),
            f(COL_QUERY_KEY, FilterOp::IsIn, vec![text("qB"), text("qC")]),
        ]);
        assert!(p.short_circuit);
        assert_eq!(p.query_keys.as_ref().map(|s| s.len()), Some(0));
    }

    #[test]
    fn strand_eq_recognised_and_contradictions_short_circuit() {
        let p = plan(&[f(COL_STRAND, FilterOp::Eq, vec![text("minus")])]);
        assert_eq!(p.strand_keep, Some(StrandFilter::Minus));
        assert!(!p.short_circuit);
        assert!(p.strand_keep.unwrap().matches(Strand::Minus));
        assert!(!p.strand_keep.unwrap().matches(Strand::Plus));

        // plus AND minus can't both hold.
        let p = plan(&[
            f(COL_STRAND, FilterOp::Eq, vec![text("plus")]),
            f(COL_STRAND, FilterOp::Eq, vec![text("minus")]),
        ]);
        assert!(p.short_circuit);
    }

    #[test]
    fn unrecognised_filters_are_dropped_not_errored() {
        // query_start > 100 : column we don't prune on -> ignored.
        // strand IN (...)   : op we don't handle on strand -> ignored.
        // evalue = 1e-5     : Eq on a float column -> ignored (Lt/Le only).
        let p = plan(&[
            f(2, FilterOp::Gt, vec![Duckvalue::Uint32(100)]),
            f(COL_STRAND, FilterOp::IsIn, vec![text("plus"), text("minus")]),
            f(COL_EVALUE, FilterOp::Eq, vec![Duckvalue::Float64(1e-5)]),
        ]);
        assert_eq!(p.evalue_ceiling, None);
        assert_eq!(p.identity_floor, None);
        assert!(p.query_keys.is_none());
        assert!(p.subject_keys.is_none());
        assert!(p.strand_keep.is_none());
        assert!(!p.short_circuit);
    }
}
