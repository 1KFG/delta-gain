#!/usr/bin/env python3
"""Bin DIAMOND_NOVELTY_SCREEN hits into strong-hit / weak-hit / no-hit per
BFD query protein (design doc Step 3 framing fix: E-value + coverage, not
raw %identity -- a raw identity cutoff wrongly flags paralogs and
fast-evolving orthologs as "novel").

Takes one or more DIAMOND outfmt6+stitle TSVs, each labelled with which
target database it searched (nr_cluster, uniprot_fungi, ...), plus the full
query FASTA so proteins with NO hit at all in any target still get a row
(DIAMOND's own output only ever lists queries that hit something).

Per target, per query, keeps only the best hit (lowest E-value) and
classifies it:
  strong_hit  evalue <= --evalue-strong AND qcovhsp/scovhsp >= --cov-strong
  weak_hit    evalue <= --evalue-weak, but not strong
  no_hit      no hit at all, or hit exists but evalue > --evalue-weak

Overall per-query status is the best (strong > weak > no_hit) across all
targets -- a hit in ANY reference (nr_cluster or uniprot_fungi) counts as
"already known."

THRESHOLD DEFAULTS ARE PLACEHOLDERS (HANDOFF.md item 4: "better made
deliberately than defaulted", not yet backed by real score-distribution
data from this pipeline). evalue-strong reuses params.diamond_evalue's
existing 1e-5 search cutoff; 70% bidirectional coverage for "strong" is a
conventional profile-search rule of thumb, not a value derived from this
dataset. Revisit both once a real pilot run's hit-score distribution is in
hand (same spirit as HANDOFF.md item 2's plan to revisit the ANI thresholds
empirically after the Stage 1 pilot).
"""
import argparse
import csv
import sys
from collections import defaultdict

OUTFMT6_STITLE_COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "qcovhsp", "scovhsp", "stitle",
]


def iter_query_ids(fasta_path):
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                yield line[1:].split()[0]


def classify(evalue, qcov, scov, evalue_strong, evalue_weak, cov_strong):
    if evalue <= evalue_strong and qcov >= cov_strong and scov >= cov_strong:
        return "strong_hit"
    if evalue <= evalue_weak:
        return "weak_hit"
    return "no_hit"


RANK = {"strong_hit": 2, "weak_hit": 1, "no_hit": 0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--queries", required=True, help="deduplicated BFD query FASTA")
    ap.add_argument(
        "--hits", required=True, action="append",
        metavar="TARGET=PATH",
        help="one per target db, e.g. --hits nr_cluster=chunk1.nr_cluster.tsv "
             "--hits uniprot_fungi=chunk1.uniprot_fungi.tsv (repeat per chunk "
             "and target; same target label accumulates across chunks)",
    )
    ap.add_argument("--evalue-strong", type=float, default=1e-5)
    ap.add_argument("--evalue-weak", type=float, default=1e-3)
    ap.add_argument("--cov-strong", type=float, default=70.0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    targets = sorted({spec.split("=", 1)[0] for spec in args.hits})

    # best[target][qseqid] = (evalue, qcov, scov)
    best = {t: {} for t in targets}
    for spec in args.hits:
        target, path = spec.split("=", 1)
        with open(path, newline="") as fh:
            for line in fh:
                if not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                row = dict(zip(OUTFMT6_STITLE_COLS, fields))
                q = row["qseqid"]
                evalue = float(row["evalue"])
                qcov = float(row["qcovhsp"])
                scov = float(row["scovhsp"])
                cur = best[target].get(q)
                if cur is None or evalue < cur[0]:
                    best[target][q] = (evalue, qcov, scov)

    query_ids = list(iter_query_ids(args.queries))

    counts = defaultdict(int)
    with open(args.out, "w", newline="") as out_fh:
        w = csv.writer(out_fh, delimiter="\t")
        header = ["query_id"]
        for t in targets:
            header += [f"{t}_status", f"{t}_evalue", f"{t}_qcov", f"{t}_scov"]
        header.append("overall_status")
        w.writerow(header)

        for q in query_ids:
            row = [q]
            overall = "no_hit"
            for t in targets:
                hit = best[t].get(q)
                if hit is None:
                    status, ev, qc, sc = "no_hit", "", "", ""
                else:
                    evalue, qcov, scov = hit
                    status = classify(
                        evalue, qcov, scov,
                        args.evalue_strong, args.evalue_weak, args.cov_strong,
                    )
                    ev, qc, sc = evalue, qcov, scov
                row += [status, ev, qc, sc]
                if RANK[status] > RANK[overall]:
                    overall = status
            row.append(overall)
            counts[overall] += 1
            w.writerow(row)

    print(
        f"bin_novelty_hits: {len(query_ids)} queries -> "
        f"strong_hit={counts['strong_hit']} weak_hit={counts['weak_hit']} "
        f"no_hit={counts['no_hit']} (targets: {', '.join(targets)})",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
