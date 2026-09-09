#!/usr/bin/env python3
"""CLI wrapping lib_accumulation.py: writes the accumulation curve summary,
pan-proteome power-law fit, per-clade Shapley contribution, and raw
per-permutation marginal records (design doc: Components, Pipeline
integration)."""
import argparse
import csv

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from lib_accumulation import (
    run_all_permutations,
    aggregate_clade_contribution,
    fit_power_law_per_permutation,
)


def _load_inputs(incidence_path, external_status_path, genome_metadata_path):
    incidence = pq.read_table(incidence_path).to_pylist()
    asmid_clusters = {}
    for row in incidence:
        asmid_clusters.setdefault(row["asmid"], set()).add(row["cluster_id"])

    ext_status_table = pq.read_table(external_status_path).to_pylist()
    external_status = {r["cluster_id"]: r["external_status"] for r in ext_status_table}

    genome_metadata = pq.read_table(genome_metadata_path).to_pylist()
    clade_of = {r["asmid"]: (r["clade_rank"], r["clade_label"]) for r in genome_metadata}

    return asmid_clusters, external_status, clade_of


def _write_summary_tsv(walks, path):
    n_genomes = len(walks[0])
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["N", "internal_mean", "internal_p2_5", "internal_p97_5",
                          "external_mean", "external_p2_5", "external_p97_5"])
        for position in range(1, n_genomes + 1):
            internal_vals = np.array([w[position - 1]["internal_cumulative"] for w in walks], dtype=float)
            external_vals = np.array([w[position - 1]["external_cumulative"] for w in walks], dtype=float)
            writer.writerow([
                position,
                internal_vals.mean(), np.percentile(internal_vals, 2.5), np.percentile(internal_vals, 97.5),
                external_vals.mean(), np.percentile(external_vals, 2.5), np.percentile(external_vals, 97.5),
            ])


def _write_fit_tsv(walks, path):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["series", "alpha_mean", "alpha_ci_low", "alpha_ci_high",
                          "kappa_mean", "r_squared_mean"])
        for series, count_key in [("internal", "new_internal_count"), ("external", "new_external_count")]:
            fit = fit_power_law_per_permutation(walks, count_key)
            writer.writerow([series, fit["alpha_mean"], fit["alpha_ci_low"], fit["alpha_ci_high"],
                              fit["kappa_mean"], fit["r_squared_mean"]])


def _write_clade_tsv(walks, clade_of, path):
    rows = aggregate_clade_contribution(walks, clade_of)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _write_raw_marginals_parquet(walks, permutation_ids, path):
    records = []
    for perm_id, walk in zip(permutation_ids, walks):
        for record in walk:
            records.append({"permutation_id": perm_id, **record})
    table = pa.table({key: [r[key] for r in records] for key in records[0]})
    pq.write_table(table, path, compression="zstd")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--incidence-matrix", required=True)
    ap.add_argument("--external-status", required=True)
    ap.add_argument("--genome-metadata", required=True)
    ap.add_argument("--n-permutations", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--version-tag", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    asmid_clusters, external_status, clade_of = _load_inputs(
        args.incidence_matrix, args.external_status, args.genome_metadata
    )
    asmids = sorted(asmid_clusters)
    walks = run_all_permutations(asmids, asmid_clusters, external_status,
                                  args.n_permutations, args.seed)
    permutation_ids = list(range(args.n_permutations))

    _write_summary_tsv(walks, f"{args.outdir}/accumulation_curve_summary.tsv")
    _write_fit_tsv(walks, f"{args.outdir}/pangenome_powerlaw_fit.tsv")
    _write_clade_tsv(walks, clade_of, f"{args.outdir}/clade_contribution.tsv")
    _write_raw_marginals_parquet(walks, permutation_ids, f"{args.outdir}/permutation_marginals.parquet")


if __name__ == "__main__":
    main()
