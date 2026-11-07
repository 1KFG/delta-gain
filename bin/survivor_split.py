#!/usr/bin/env python3
"""Split tier-50 representatives for the survivor search (docs/DELTA_GAIN_PLAN.md, stage 5).

query  = representatives of survivor clusters (no baseline member at 50 %) seen in >= --min-genomes genomes
target = representatives of clusters that DO contain a baseline protein (all tier-50 reps minus survivor reps)
Streams the FASTA once; memory ~ one id set per survivor class.
"""
import argparse, csv, gzip, sys

ap = argparse.ArgumentParser()
ap.add_argument("--reps", required=True, help="tier_50_reps.faa")
ap.add_argument("--survivors", required=True, help="survivors_clusters.tsv (rep, n_new_proteins, n_genomes)")
ap.add_argument("--min-genomes", type=int, default=5)
ap.add_argument("--max-genomes", type=int, default=0, help="0 = no upper bound")
ap.add_argument("--sample", type=int, default=0, help="random sample of N eligible queries (seed 5)")
ap.add_argument("--max-query", type=int, default=0, help="0 = all; else keep the first N queries (test runs)")
ap.add_argument("--out-query", required=True)
ap.add_argument("--out-target", required=True)
a = ap.parse_args()
qset, allsurv = set(), set()
with open(a.survivors, newline="") as fh:
    rd = csv.reader(fh, delimiter="\t"); next(rd)
    for rep, n, g in rd:
        allsurv.add(rep)
        if int(g) >= a.min_genomes and (not a.max_genomes or int(g) <= a.max_genomes) and \
           (not a.max_query or len(qset) < a.max_query):
            qset.add(rep)
if a.sample and len(qset) > a.sample:
    import random
    random.seed(5)
    qset = set(random.sample(sorted(qset), a.sample))
print("survivor clusters", len(allsurv), "queries", len(qset), file=sys.stderr)
nq = nt = 0
out = None
with open(a.reps) as fh, open(a.out_query, "w") as fq, open(a.out_target, "w") as ft:
    for line in fh:
        if line.startswith(">"):
            rid = line[1:].split()[0]
            if rid in qset:
                out = fq; nq += 1
            elif rid in allsurv:
                out = None
            else:
                out = ft; nt += 1
        if out is not None:
            out.write(line)
print("wrote queries", nq, "targets", nt, file=sys.stderr)
