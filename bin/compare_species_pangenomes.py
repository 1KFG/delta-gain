#!/usr/bin/env python3
"""Cross-species comparison of pangenomes (follows bin/species_pangenome.py).

Inputs: the incidence matrix, samples.csv (taxonomy), and the table written by
species_pangenome.py. The species set is taken from that table and each
species' genome count is re-derived here and must match the table's n_genomes,
so a different --exclude-asmids between the two steps fails loudly.

Outputs (in --outdir):
  species_cluster_sharing.tsv  one row per species pair sharing >= 1 cluster:
      shared_clusters / jaccard over the species' FULL cluster sets (exact, via
      a sparse species x cluster product), and shared_clusters_narrow /
      jaccard_narrow over clusters present in at most --max-species-per-cluster
      species (broad clusters say little about relatedness), plus same_genus /
      same_family flags.
  taxon_pangenome_summary.tsv  per taxon at each --ranks level: species counts
      by openness (open / closed / unresolved), fraction open among resolved
      species (closed species stay in the denominator), medians of alpha,
      core_per_mean_genome, rarefied_pan_size and rarefied_core.
  lineage_conflicts.tsv        species whose genomes disagree on a taxonomic rank
      (the majority label is used).
"""
import argparse
import csv
import sys
from collections import Counter
from statistics import median

import numpy as np
import scipy.sparse as sp

from species_pangenome import load_pairs, load_species, load_exclude, group_species

RANKS = ["GENUS", "FAMILY", "ORDER", "CLASS", "PHYLUM"]


def load_lineage(samples_csv):
    per = {}
    with open(samples_csv) as fh:
        for row in csv.DictReader(fh):
            name = (row.get("SPECIES") or "").strip()
            if name:
                per.setdefault(name, []).append(tuple((row.get(r) or "").strip() for r in RANKS))
    lineage, conflicts = {}, []
    for name, tuples in per.items():
        lin = {}
        for i, r in enumerate(RANKS):
            c = Counter(t[i] for t in tuples)
            lin[r] = c.most_common(1)[0][0]
            if len([k for k in c if k]) > 1:
                conflicts.append({"species": name, "rank": r, "labels": ";".join(f"{k}:{v}" for k, v in c.items())})
        lineage[name] = lin
    return lineage, conflicts


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def taxon_summary(rows, lineage, ranks):
    out = []
    for rank in ranks:
        groups = {}
        for r in rows:
            label = lineage.get(r["species"], {}).get(rank, "")
            if label:
                groups.setdefault(label, []).append(r)
        for label, rs in sorted(groups.items()):
            cnt = Counter(r["openness"] for r in rs)
            resolved = cnt["open"] + cnt["closed"]
            fl = lambda k: [float(r[k]) for r in rs if r[k] != ""]
            med = lambda v: median(v) if v else ""
            out.append({
                "rank": rank, "taxon": label, "n_species": len(rs),
                "n_open": cnt["open"], "n_closed": cnt["closed"], "n_unresolved": cnt["unresolved"],
                "frac_open": cnt["open"] / resolved if resolved else "",
                "median_alpha": med(fl("alpha")),
                "median_core_per_mean_genome": med(fl("core_per_mean_genome")),
                "median_rarefied_pan_size": med(fl("rarefied_pan_size")),
                "median_rarefied_core": med(fl("rarefied_core")),
            })
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--incidence-matrix", required=True)
    ap.add_argument("--samples-csv", required=True)
    ap.add_argument("--species-table", required=True, help="output of species_pangenome.py")
    ap.add_argument("--exclude-asmids", default="")
    ap.add_argument("--max-species-per-cluster", type=int, default=100)
    ap.add_argument("--ranks", default=",".join(RANKS))
    ap.add_argument("--block", type=int, default=500, help="species per sparse-product block")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    with open(args.species_table) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    sp_names = sorted(r["species"] for r in rows)
    n_expected = {r["species"]: int(r["n_genomes"]) for r in rows}

    a_codes, c_codes, asm_names, n_clusters = load_pairs(args.incidence_matrix)
    groups, _, _ = group_species(asm_names, load_species(args.samples_csv),
                                 load_exclude(args.exclude_asmids), 1, True)
    bad = [(s, n_expected[s], len(groups.get(s, []))) for s in sp_names
           if len(groups.get(s, [])) != n_expected[s]]
    if bad:
        sys.exit("species genome counts differ from the species table (different "
                 f"--exclude-asmids or inputs?): {bad[:5]} (species, table, here)")

    sp_of_asm = np.full(len(asm_names), -1, dtype=np.int32)
    for i, s in enumerate(sp_names):
        sp_of_asm[groups[s]] = i
    s_pair = sp_of_asm[a_codes]
    keep = s_pair >= 0
    ones = np.ones(int(keep.sum()), dtype=np.int32)
    S = sp.coo_matrix((ones, (s_pair[keep], c_codes[keep])), shape=(len(sp_names), n_clusters)).tocsr()
    S.sum_duplicates()
    S.data[:] = 1                                   # presence/absence per species
    n_species_per_cluster = np.asarray(S.sum(axis=0)).ravel()
    narrow_cols = np.flatnonzero((n_species_per_cluster >= 1) & (n_species_per_cluster <= args.max_species_per_cluster))
    Sn = S[:, narrow_cols].tocsr()
    pan = np.asarray(S.sum(axis=1)).ravel()
    pan_n = np.asarray(Sn.sum(axis=1)).ravel()
    ST, SnT = S.T.tocsr(), Sn.T.tocsr()

    lineage, conflicts = load_lineage(args.samples_csv)
    write_tsv(f"{args.outdir}/lineage_conflicts.tsv", conflicts, ["species", "rank", "labels"])
    if conflicts:
        print(f"WARNING: {len(conflicts)} species/rank lineage conflicts (see lineage_conflicts.tsv)", file=sys.stderr)

    fields = ["species_a", "species_b", "shared_clusters", "jaccard", "shared_clusters_narrow",
              "jaccard_narrow", "same_genus", "same_family"]
    n_pairs = 0
    with open(f"{args.outdir}/species_cluster_sharing.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(fields)
        for i0 in range(0, len(sp_names), args.block):
            i1 = min(i0 + args.block, len(sp_names))
            B = (S[i0:i1] @ ST).tocoo()
            Bn = (Sn[i0:i1] @ SnT).tocsr()
            m = B.col > (B.row + i0)                 # upper triangle only
            r, c, v = B.row[m] + i0, B.col[m], B.data[m]
            vn = np.asarray(Bn[r - i0, c]).ravel()
            for ri, ci, sh, shn in zip(r, c, v, vn):
                a, b = sp_names[ri], sp_names[ci]
                la, lb = lineage[a], lineage[b]
                w.writerow([a, b, int(sh), f"{sh / (pan[ri] + pan[ci] - sh):.6g}", int(shn),
                            f"{(shn / (pan_n[ri] + pan_n[ci] - shn)) if (pan_n[ri] + pan_n[ci] - shn) else 0:.6g}",
                            int(bool(la['GENUS']) and la['GENUS'] == lb['GENUS']),
                            int(bool(la['FAMILY']) and la['FAMILY'] == lb['FAMILY'])])
                n_pairs += 1
    ts = taxon_summary(rows, lineage, [x for x in args.ranks.split(",") if x])
    write_tsv(f"{args.outdir}/taxon_pangenome_summary.tsv", ts,
              ["rank", "taxon", "n_species", "n_open", "n_closed", "n_unresolved", "frac_open",
               "median_alpha", "median_core_per_mean_genome", "median_rarefied_pan_size",
               "median_rarefied_core"])
    print(f"species: {len(sp_names)}; pairs sharing >=1 cluster: {n_pairs}; taxon rows: {len(ts)}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
