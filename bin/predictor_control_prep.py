#!/usr/bin/env python3
"""Predictor-mismatch control for the funnel (docs/DELTA_GAIN_PLAN.md §9; reviewer risk #1).

The funnel compares NEW genomes annotated by the BFD pipeline against BASELINE genomes annotated by
NCBI. If gene-model differences alone create "novel" proteins, the headline is inflated. Control:
take T baseline genomes that have BOTH an NCBI protein set and a BFD re-annotation, and a reference R of
other baseline genomes (NCBI proteins, the test genomes are NOT in R). Run the funnel twice against the same R:
  A: query = NCBI proteins of the T genomes     (baseline-style predictor)
  B: query = BFD proteins of the same T genomes (new-genome predictor)
Per-genome novelty(B) - novelty(A) is the predictor effect (same genomes, same reference).
Writes status_ref.tsv, manifest_A.tsv, manifest_B.tsv in --outdir.
"""
import argparse, csv, os, random

ap = argparse.ArgumentParser()
ap.add_argument("--status", required=True)
ap.add_argument("--classes", required=True)
ap.add_argument("--manifest", required=True, help="bfd_pep_manifest.tsv")
ap.add_argument("--asm-dir", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--test", type=int, default=300)
ap.add_argument("--ref", type=int, default=1000)
ap.add_argument("--seed", type=int, default=11)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

bfd = {r["ASMID"].strip(): r["pep_path"].strip() for r in csv.DictReader(open(a.manifest), delimiter="\t")}
cls = {r["ASMID"].strip(): (r["class"].strip(), r["pep_status"].strip())
       for r in csv.DictReader(open(a.classes), delimiter="\t")}
stat = {}
for r in csv.DictReader(open(a.status), delimiter="\t"):
    f = r["folder"].strip()
    if r["status"] == "baseline_ok" and r["protein_file_ok"] == "1" and \
       os.path.exists(os.path.join(a.asm_dir, f, f + "_protein.faa.gz")):
        stat[f] = r
cand = sorted(g for g, (c, p) in cls.items() if c == "baseline" and p == "matched" and g in stat
              and g in bfd and os.path.exists(bfd[g]))
print("baseline genomes with NCBI + BFD proteins:", len(cand))
assert len(cand) > 4000, "ASMID/folder mismatch?"
random.seed(a.seed)
random.shuffle(cand)
test, ref = cand[:a.test], cand[a.test:a.test + a.ref]
rows = [stat[g] for g in ref]
with open(os.path.join(a.outdir, "status_ref.tsv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)
for lab, path_of in (("A", lambda g: os.path.join(a.asm_dir, g, g + "_protein.faa.gz")), ("B", lambda g: bfd[g])):
    with open(os.path.join(a.outdir, f"manifest_{lab}.tsv"), "w") as fh:
        fh.write("ASMID\ttag\tpep_path\tstatus\n")
        for g in test:
            fh.write(f"{g}\t{g}\t{path_of(g)}\tmatched\n")
print("test", len(test), "ref", len(ref))
