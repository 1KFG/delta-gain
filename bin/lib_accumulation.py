#!/usr/bin/env python3
"""Permutation-resampling core for the protein-content accumulation curve
(design doc: Computational strategy, Novelty definition, Components).

Every statistic in this analysis is computed by resampling over a
precomputed genome x cluster incidence matrix -- no bioinformatics tool is
re-run per permutation.
"""
import random


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
