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
