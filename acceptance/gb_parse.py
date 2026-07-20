#!/usr/bin/env python3
"""
Minimal GenBank parser for the synphage-ducklink acceptance test.

Extracts, from a set of .gb files, one CDS-level DataFrame we can feed to
BLASTN and to the `genome_features` table the SQL views expect:

  * `genome_id`       -- LOCUS name (e.g. NC_001416)
  * `feature_key`     -- unique across all genomes (locus_tag preferred,
                          otherwise `{genome_id}:{gene}:{i}` fallback)
  * `feature_type`    -- always 'CDS' (only feature we emit)
  * `start_position`  -- 1-based inclusive
  * `end_position`    -- 1-based inclusive
  * `strand`          -- +1 for forward, -1 for complement
  * `gene`            -- /gene qualifier or NULL
  * `product`         -- /product qualifier or NULL
  * `nucleotide_sequence` -- the CDS DNA sequence (reverse-complemented on
                             complement strand)

This is NOT a general-purpose GenBank parser -- it deliberately handles only
the subset of the LOCUS / FEATURES / ORIGIN grammar needed for
well-formed NCBI reference phage genomes:

  * simple ranges         `191..736`
  * complement wrap       `complement(20147..20767)`
  * join wrap             `join(1..10, 20..30)`
  * complement(join(...)) recognised as join across complement

Fuzzy positions (`<1..736`, `191..>736`) are treated as exact by stripping
the < / >. Nested join(complement(...)) and non-contiguous joins are
supported by concatenating the extracted segments in order.

The equivalent code in the future `components/genome-format/` DuckLink
bridge will be the real parser; this driver script exists so the
acceptance test doesn't depend on that bridge being ready.
"""
from __future__ import annotations

import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

COMPLEMENT_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT_TABLE)[::-1]


@dataclass
class Cds:
    genome_id: str
    feature_key: str
    feature_type: str  # always 'CDS'
    start_position: int
    end_position: int
    strand: int  # +1 or -1
    gene: str | None
    product: str | None
    nucleotide_sequence: str


# ---------------------------------------------------------------------------
# Location parser
# ---------------------------------------------------------------------------

_LOCATION_RANGE = re.compile(r"<?(?P<start>\d+)\.\.>?(?P<end>\d+)")
_SINGLE_POS = re.compile(r"^\d+$")


def _extract_ranges(loc: str) -> list[tuple[int, int]]:
    """Return list of (start, end) 1-based inclusive spans.

    Ignores strand -- the caller handles complement separately.
    """
    inner = loc
    # peel one layer of complement()
    if inner.startswith("complement(") and inner.endswith(")"):
        inner = inner[len("complement(") : -1]
    if inner.startswith("join(") and inner.endswith(")"):
        inner = inner[len("join(") : -1]
    ranges: list[tuple[int, int]] = []
    for part in inner.split(","):
        part = part.strip()
        # Strip nested complement inside a join(complement(),...) — for our
        # purposes we handle strand at the join level, so just take the
        # inner range coordinates.
        if part.startswith("complement(") and part.endswith(")"):
            part = part[len("complement(") : -1]
        m = _LOCATION_RANGE.match(part)
        if m:
            ranges.append((int(m.group("start")), int(m.group("end"))))
            continue
        if _SINGLE_POS.match(part):
            p = int(part)
            ranges.append((p, p))
    return ranges


def _location_is_complement(loc: str) -> bool:
    return loc.startswith("complement(")


# ---------------------------------------------------------------------------
# Feature block scanner
# ---------------------------------------------------------------------------


def _iter_feature_blocks(lines: list[str]) -> Iterator[tuple[str, str, list[str]]]:
    """Yield (feature_type, location_string, qualifier_lines) from a
    FEATURES ... ORIGIN section.

    A feature record spans:
      * one header line: `     CDS             191..736`
      * zero+ location-continuation lines (unquoted text, no leading '/')
      * one+ qualifier lines: `                     /gene="A"` and their
        continuation lines (which do NOT start with '/', so we distinguish
        by whether any qualifier has begun yet in this record).
    """
    current: list[str] | None = None
    current_type: str | None = None
    current_location: str | None = None
    saw_qualifier = False  # once True, non-'/' continuations belong to the
                           # LAST qualifier, not to the location.

    for line in lines:
        if len(line) < 21:
            continue
        # feature header: cols 6..21 hold the feature key
        head = line[5:20].strip()
        rest = line[21:].rstrip()
        if head:  # new feature key
            # flush previous record
            if (
                current is not None
                and current_type is not None
                and current_location is not None
            ):
                yield (current_type, current_location, current)
            current = []
            current_type = head
            current_location = rest
            saw_qualifier = False
        else:
            # continuation line for the record in progress
            if current is None or current_type is None:
                continue
            body = line[21:].rstrip()
            if body.startswith("/"):
                current.append(body)
                saw_qualifier = True
            elif saw_qualifier:
                # continuation of the previous qualifier value
                if current:
                    current[-1] = current[-1] + body.strip()
            else:
                # location wrap -- append to location, no separator
                current_location = (current_location or "") + body.strip()

    if (
        current is not None
        and current_type is not None
        and current_location is not None
    ):
        yield (current_type, current_location, current)


# ---------------------------------------------------------------------------
# GenBank file parser
# ---------------------------------------------------------------------------

_QUAL_RE = re.compile(r'/(?P<key>\w+)=(?P<val>"[^"]*"|[^\s]+)')


def _extract_qualifier(qual_lines: list[str], key: str) -> str | None:
    """Return the value of /key=... concatenated across continuation lines
    (GenBank sometimes wraps long /product strings), stripped of surrounding
    quotes. Returns None if key not present.
    """
    prefix = f"/{key}="
    hit = None
    collected: list[str] | None = None
    for line in qual_lines:
        s = line.strip()
        if collected is not None:
            # still collecting a wrapped value
            if s.startswith("/"):
                # start of a new qualifier -- stop collecting
                break
            collected.append(s)
            if s.endswith('"'):
                break
            continue
        if s.startswith(prefix):
            hit = s[len(prefix) :]
            if hit.startswith('"') and hit.endswith('"') and len(hit) >= 2:
                return hit[1:-1]
            if hit.startswith('"'):
                # begin multi-line value
                collected = [hit]
                continue
            return hit
    if collected:
        joined = " ".join(collected)
        if joined.startswith('"') and joined.endswith('"'):
            return joined[1:-1]
        return joined
    return None


def parse_gb(path: Path) -> list[Cds]:
    """Return the list of CDS records extracted from a single GenBank file.

    The `feature_key` is a genome-globally unique string: locus_tag if
    present, else `{genome_id}_cds{ordinal}`.
    """
    text = path.read_text()
    # split into header / FEATURES / ORIGIN
    m_features = re.search(r"^FEATURES\s+.*$", text, re.MULTILINE)
    m_origin = re.search(r"^ORIGIN.*$", text, re.MULTILINE)
    if m_features is None or m_origin is None:
        raise ValueError(f"{path}: missing FEATURES or ORIGIN section")
    features_block = text[m_features.end() : m_origin.start()].splitlines()

    # sequence: everything after ORIGIN up to the next // line
    origin_body = text[m_origin.end() :]
    stop = origin_body.find("//")
    if stop >= 0:
        origin_body = origin_body[:stop]
    seq_chars = re.sub(r"[^A-Za-z]", "", origin_body).upper()

    # LOCUS is the primary genome id; fall back to file stem
    m_locus = re.search(r"^LOCUS\s+(\S+)", text, re.MULTILINE)
    genome_id = m_locus.group(1) if m_locus else path.stem

    cds_records: list[Cds] = []
    ordinal = 0
    for feature_type, location, quals in _iter_feature_blocks(features_block):
        if feature_type != "CDS":
            continue
        ordinal += 1
        ranges = _extract_ranges(location)
        if not ranges:
            continue
        complement = _location_is_complement(location)
        # extract concatenated DNA
        parts: list[str] = []
        for s, e in ranges:
            if s < 1 or e > len(seq_chars) or s > e:
                # skip malformed range
                continue
            parts.append(seq_chars[s - 1 : e])
        dna = "".join(parts)
        if not dna:
            continue
        if complement:
            dna = revcomp(dna)
        gene = _extract_qualifier(quals, "gene")
        product = _extract_qualifier(quals, "product")
        locus_tag = _extract_qualifier(quals, "locus_tag")
        feature_key = locus_tag or f"{genome_id}_cds{ordinal:03d}"
        cds_records.append(
            Cds(
                genome_id=genome_id,
                feature_key=feature_key,
                feature_type="CDS",
                start_position=min(r[0] for r in ranges),
                end_position=max(r[1] for r in ranges),
                strand=-1 if complement else 1,
                gene=gene,
                product=product,
                nucleotide_sequence=dna,
            )
        )
    return cds_records


# ---------------------------------------------------------------------------
# Selection heuristic
# ---------------------------------------------------------------------------


def pick_cds(
    all_cds: list[Cds], per_genome: int, max_len: int
) -> list[Cds]:
    """Sample `per_genome` CDS from each genome, preferring lengths in
    [~150, max_len].

    Very short CDS (< 150 bp) tend to be spurious short ORFs; very long ones
    make the pairwise Smith-Waterman cost explode (O(m*n) memory and time),
    so we cap length. Within the accepted range we take the shortest to
    keep the run fast — homologous phage genes still land in this bucket
    (many capsid / tail / lysis genes are 300-900 bp).
    """
    by_genome: dict[str, list[Cds]] = {}
    for c in all_cds:
        by_genome.setdefault(c.genome_id, []).append(c)
    out: list[Cds] = []
    for gid, items in by_genome.items():
        eligible = [
            c
            for c in items
            if 150 <= len(c.nucleotide_sequence) <= max_len
        ]
        eligible.sort(key=lambda c: len(c.nucleotide_sequence))
        out.extend(eligible[:per_genome])
    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _emit_features_tsv(cds_records: Iterable[Cds], path: Path) -> None:
    header = [
        "genome_id",
        "feature_key",
        "feature_type",
        "start_position",
        "end_position",
        "strand",
        "gene",
        "product",
    ]
    with path.open("w") as f:
        f.write("\t".join(header) + "\n")
        for c in cds_records:
            f.write(
                "\t".join(
                    [
                        c.genome_id,
                        c.feature_key,
                        c.feature_type,
                        str(c.start_position),
                        str(c.end_position),
                        str(c.strand),
                        c.gene or "",
                        (c.product or "").replace("\t", " ").replace("\n", " "),
                    ]
                )
                + "\n"
            )


def _emit_seq_json(cds_records: Iterable[Cds], path: Path) -> None:
    payload = [
        {"key": c.feature_key, "data": c.nucleotide_sequence}
        for c in cds_records
    ]
    path.write_text(json.dumps(payload))


def main(argv: list[str]) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument(
        "--data-dir",
        default="acceptance/data",
        help="directory containing .gb files",
    )
    ap.add_argument(
        "--out-dir",
        default="acceptance/build",
        help="directory to write queries.json / subjects.json / features.tsv",
    )
    ap.add_argument(
        "--per-genome",
        type=int,
        default=8,
        help="max number of CDS to sample from each genome",
    )
    ap.add_argument(
        "--max-len",
        type=int,
        default=900,
        help="skip CDS longer than this (bp)",
    )
    args = ap.parse_args(argv)

    data_dir = Path(args.data_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    gb_files = sorted(data_dir.glob("*.gb")) + sorted(data_dir.glob("*.gbk"))
    if not gb_files:
        print(
            f"ERROR: no *.gb or *.gbk files under {data_dir}",
            file=sys.stderr,
        )
        return 1

    print(f"[gb_parse] parsing {len(gb_files)} GenBank file(s) from {data_dir}")
    all_cds: list[Cds] = []
    for gb in gb_files:
        records = parse_gb(gb)
        print(
            f"[gb_parse]   {gb.name}: {len(records)} CDS "
            f"(seq lengths min={min((len(r.nucleotide_sequence) for r in records), default=0)}, "
            f"max={max((len(r.nucleotide_sequence) for r in records), default=0)})"
        )
        all_cds.extend(records)

    print(
        f"[gb_parse] total CDS parsed: {len(all_cds)} across "
        f"{len({c.genome_id for c in all_cds})} genome(s)"
    )
    selected = pick_cds(
        all_cds, per_genome=args.per_genome, max_len=args.max_len
    )
    print(
        f"[gb_parse] sampled {len(selected)} CDS "
        f"(per-genome cap={args.per_genome}, max-len={args.max_len})"
    )
    # Emit BOTH the full annotation table (all CDS) and the sampled
    # queries / subjects. genome_features.tsv is used by the SQL views
    # via `feature_key` lookup, so it must contain every subject we
    # BLAST against; but the BLAST payload is the sampled subset only.
    _emit_features_tsv(selected, out_dir / "genome_features.tsv")
    # queries == subjects for the all-vs-all conservation query
    _emit_seq_json(selected, out_dir / "queries.json")
    _emit_seq_json(selected, out_dir / "subjects.json")
    print(
        f"[gb_parse] wrote {out_dir}/{{queries,subjects}}.json "
        f"and {out_dir}/genome_features.tsv"
    )
    # summary per-genome
    per_g: dict[str, int] = {}
    for c in selected:
        per_g[c.genome_id] = per_g.get(c.genome_id, 0) + 1
    for gid, n in sorted(per_g.items()):
        print(f"[gb_parse]   {gid}: {n} CDS in sample")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
