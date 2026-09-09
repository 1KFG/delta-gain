#!/usr/bin/env python3
"""Permutation-resampling core for the protein-content accumulation curve
(design doc: Computational strategy, Novelty definition, Components).

Every statistic in this analysis is computed by resampling over a
precomputed genome x cluster incidence matrix -- no bioinformatics tool is
re-run per permutation.
"""
import random
import numpy as np


def run_permutation(asmid_order, asmid_clusters, external_status):
    seen_clusters = set()
    seen_external_count = 0
    records = []
    for position, asmid in enumerate(asmid_order, start=1):
        genome_clusters = asmid_clusters[asmid]
        new_clusters = genome_clusters - seen_clusters
        new_internal_count = len(new_clusters)
        new_external_count = sum(
            1 for c in new_clusters if external_status.get(c) == "no_hit"
        )
        seen_clusters |= new_clusters
        seen_external_count += new_external_count
        records.append({
            "position": position,
            "asmid": asmid,
            "internal_cumulative": len(seen_clusters),
            "external_cumulative": seen_external_count,
            "new_internal_count": new_internal_count,
            "new_external_count": new_external_count,
        })
    return records


def run_all_permutations(asmids, asmid_clusters, external_status,
                          n_permutations, seed):
    rng = random.Random(seed)
    walks = []
    for _ in range(n_permutations):
        order = list(asmids)
        rng.shuffle(order)
        walks.append(run_permutation(order, asmid_clusters, external_status))
    return walks


def aggregate_clade_contribution(walks, clade_of):
    n_permutations = len(walks)
    genome_marginal_sum = {}  # asmid -> [sum_internal, sum_external]
    for walk in walks:
        for record in walk:
            asmid = record["asmid"]
            acc = genome_marginal_sum.setdefault(asmid, [0.0, 0.0])
            acc[0] += record["new_internal_count"]
            acc[1] += record["new_external_count"]

    # Genome-level Shapley value = mean marginal contribution across permutations.
    genome_shapley = {
        asmid: (total[0] / n_permutations, total[1] / n_permutations)
        for asmid, total in genome_marginal_sum.items()
    }

    clade_groups = {}  # (rank, label) -> list of asmid
    for asmid, key in clade_of.items():
        clade_groups.setdefault(key, []).append(asmid)

    rows = []
    for (rank, label), asmids in clade_groups.items():
        internal_vals = [genome_shapley[a][0] for a in asmids]
        external_vals = [genome_shapley[a][1] for a in asmids]
        rows.append({
            "clade_rank": rank,
            "clade_label": label,
            "n_genomes": len(asmids),
            "mean_marginal_internal": sum(internal_vals) / len(asmids),
            "total_marginal_internal": sum(internal_vals),
            "mean_marginal_external": sum(external_vals) / len(asmids),
            "total_marginal_external": sum(external_vals),
        })
    return rows


def _fit_one_walk(walk, count_key):
    """Log-log OLS fit of count_key(N) ~= kappa * N^-alpha for one
    permutation's marginal series. Points with count==0 are excluded (log(0)
    is undefined -- see design doc validation notes on this edge case)."""
    positions = np.array([r["position"] for r in walk], dtype=float)
    counts = np.array([r[count_key] for r in walk], dtype=float)
    mask = counts > 0
    if mask.sum() < 2:
        return None
    log_n = np.log(positions[mask])
    log_count = np.log(counts[mask])
    # log(count) = log(kappa) - alpha * log(N)  ->  linear regression.
    slope, intercept = np.polyfit(log_n, log_count, 1)
    alpha = -slope
    kappa = np.exp(intercept)
    predicted = intercept + slope * log_n
    ss_res = np.sum((log_count - predicted) ** 2)
    ss_tot = np.sum((log_count - log_count.mean()) ** 2)
    r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return alpha, kappa, r_squared


def fit_power_law_per_permutation(walks, count_key):
    fits = [f for f in (_fit_one_walk(w, count_key) for w in walks) if f is not None]
    alphas = np.array([f[0] for f in fits])
    kappas = np.array([f[1] for f in fits])
    r_squareds = np.array([f[2] for f in fits])
    return {
        "alpha_mean": float(alphas.mean()),
        "alpha_ci_low": float(np.percentile(alphas, 2.5)),
        "alpha_ci_high": float(np.percentile(alphas, 97.5)),
        "kappa_mean": float(kappas.mean()),
        "r_squared_mean": float(r_squareds.mean()),
        "alphas": alphas.tolist(),
    }
