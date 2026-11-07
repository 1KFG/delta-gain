#!/usr/bin/env python3
"""Print family_key<TAB>group for the F4 QC: enriched_pfam (top 200 by F0 species enrichment),
control_pfam (results/pfam_control_families.txt), cazy, merops (all extracted families of the wide run)."""
import csv, sys, os
R = sys.argv[1]
out = {}
rows = [r for r in csv.DictReader(open(f"{R}/family_census/pfam_census.tsv"), delimiter="\t")
        if r["q_species"] != "NA" and int(r["g_new"]) >= 100 and float(r["q_species"]) < 0.01]
rows.sort(key=lambda r: -float(r["log2_prev_ratio_species"]))
for r in rows[:200]:
    out["pfam_" + r["family"]] = "enriched_pfam"
for f in open(f"{R}/pfam_control_families.txt").read().split():
    out[f.replace(":", "_")] = "control_pfam"
for r in csv.DictReader(open(f"{R}/family_regions_wide/regions_summary.tsv"), delimiter="\t"):
    k = r["family"]
    if k.startswith("cazy_"):
        out[k] = "cazy"
    elif k.startswith("merops_"):
        out[k] = "merops"
for k, g in sorted(out.items()):
    print(f"{k}\t{g}")
