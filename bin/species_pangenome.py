#!/usr/bin/env python3
"""Per-species pangenome summaries from the genome x cluster incidence matrix.

Groups genomes by samples.csv SPECIES, then for each species with at least
--min-genomes genomes reports pangenome size, core/soft-core/shell/cloud
cluster counts, a rarefied pangenome size at a common genome count (so
species with different sampling can be compared), and the power-law fit of
the new-cluster rate (reuses lib_accumulation, same model as the global
curve: new clusters(N) ~= kappa * N^-alpha; alpha < 1 means an open
pangenome).

Clusters are whatever clustering produced the incidence matrix (currently
MMSEQS linclust 95% id / 90% cov). Pass --exclude-asmids to drop redundant
strains (e.g. Stage 1 redundant_strain calls) before counting.
"""
import argparse
import csv

import numpy as np
import pyarrow.parquet as pq

from lib_accumulation import (
    run_all_permutations,
    fit_power_law_per_permutation,
)

SOFT_CORE_FRAC = 0.95
CLOUD_FRAC = 0.15


def load_incidence(path):
    table = pq.read_table(path, columns=["asmid", "cluster_id"]).to_pydict()
    asmid_clusters = {}
    for a, c in zip(table["asmid"], table["cluster_id"]):
        asmid_clusters.setdefault(a, set()).add(c)
    return asmid_clusters


def load_species(samples_csv):
    species_of = {}
    with open(samples_csv) as fh:
        for row in csv.DictReader(fh):
            species_of[row["ASMID"]] = (row.get("SPECIES") or "").strip()
    return species_of


def load_exclude(path):
    if not path:
        return set()
    with open(path) as fh:
        return {line.strip().split("\t")[0] for line in fh if line.strip()}


def presence_classes(asmids, asmid_clusters):
    n = len(asmids)
    counts = {}
    for a in asmids:
        for c in asmid_clusters[a]:
            counts[c] = counts.get(c, 0) + 1
    vals = np.array(list(counts.values()))
    return {
        "pan_size": len(counts),
        "core": int((vals == n).sum()),
        "soft_core": int(((vals >= SOFT_CORE_FRAC * n) & (vals < n)).sum()),
        "cloud": int((vals < CLOUD_FRAC * n).sum()),
        "singleton_clusters": int((vals == 1).sum()),
    }


def summarize_species(species, asmids, asmid_clusters, n_perm, rarefy_n, seed):
    n = len(asmids)
    row = {"species": species, "n_genomes": n}
    row.update(presence_classes(asmids, asmid_clusters))
    row["shell"] = row["pan_size"] - row["core"] - row["soft_core"] - row["cloud"]
    row["core_fraction"] = row["core"] / row["pan_size"]
    row["mean_genome_clusters"] = float(np.mean([len(asmid_clusters[a]) for a in asmids]))
    walks = run_all_permutations(sorted(asmids), asmid_clusters, {}, n_perm, seed)
    if n >= rarefy_n:
        row["rarefied_pan_size"] = float(np.mean([w[rarefy_n - 1]["internal_cumulative"] for w in walks]))
    else:
        row["rarefied_pan_size"] = ""
    # A species whose genomes add no new clusters after the first (closed
    # pangenome) leaves fewer than 2 nonzero points, so no fit exists. Report
    # blanks rather than aborting the whole run.
    try:
        fit = fit_power_law_per_permutation(walks, "new_internal_count")
    except (IndexError, ValueError):
        fit = None
    for key in ("alpha_mean", "alpha_ci_low", "alpha_ci_high", "kappa_mean"):
        row[key] = fit[key] if fit else ""
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--incidence-matrix", required=True)
    ap.add_argument("--samples-csv", required=True)
    ap.add_argument("--exclude-asmids", default="", help="file of ASMIDs to drop (first column)")
    ap.add_argument("--min-genomes", type=int, default=5)
    ap.add_argument("--rarefy-n", type=int, default=5)
    ap.add_argument("--n-permutations", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    if args.rarefy_n > args.min_genomes:
        ap.error("--rarefy-n must be <= --min-genomes")

    asmid_clusters = load_incidence(args.incidence_matrix)
    species_of = load_species(args.samples_csv)
    exclude = load_exclude(args.exclude_asmids)

    by_species, unknown = {}, []
    for asmid in asmid_clusters:
        if asmid in exclude:
            continue
        sp = species_of.get(asmid, "")
        if not sp:
            unknown.append(asmid)
            continue
        by_species.setdefault(sp, []).append(asmid)

    rows = [summarize_species(sp, a, asmid_clusters, args.n_permutations,
                              args.rarefy_n, args.seed)
            for sp, a in sorted(by_species.items()) if len(a) >= args.min_genomes]
    if not rows:
        raise SystemExit("no species met --min-genomes; nothing to write")
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"species written: {len(rows)}; genomes with no species: {len(unknown)}; "
          f"excluded: {len(exclude & set(asmid_clusters))}")


if __name__ == "__main__":
    main()
