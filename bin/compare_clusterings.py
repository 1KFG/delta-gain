#!/usr/bin/env python3
"""Compare two dereplications of the SAME BFD protein set -- MMseqs2 linclust
vs. DIAMOND cluster (2026-09-08 methodology question: does DIAMOND's
DeepClust-style cascaded clustering, Buchfink et al. 2026 Nat Methods,
capture materially different structure than the linclust dereplication
this pipeline already uses?).

Both inputs are `<representative>\\t<member>` tables, no header (MMseqs2
`createtsv` and DIAMOND `cluster --out` both emit this format). Reports,
per method: cluster count, singleton count, cluster-size mean/median/max;
and between methods: the Adjusted Rand Index (ARI) over the shared member
set -- the standard measure for "do two clusterings of the same items agree,
correcting for chance," computed via the pairwise contingency table rather
than any O(n^2) pairwise enumeration (tractable at pilot scale and at full
~176M-protein scale alike).
"""
import argparse
import statistics
import sys
from collections import Counter, defaultdict


def load_membership(path):
    """Return {member_id: rep_id}."""
    member_to_rep = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            rep, member = line.rstrip("\n").split("\t")[:2]
            member_to_rep[member] = rep
    return member_to_rep


def cluster_stats(member_to_rep):
    sizes = Counter(member_to_rep.values())
    size_list = list(sizes.values())
    return {
        "n_members": len(member_to_rep),
        "n_clusters": len(sizes),
        "n_singletons": sum(1 for s in size_list if s == 1),
        "mean_size": statistics.mean(size_list),
        "median_size": statistics.median(size_list),
        "max_size": max(size_list),
    }


def choose2(n):
    return n * (n - 1) // 2


def adjusted_rand_index(a_to_rep, b_to_rep, shared_members):
    contingency = defaultdict(int)
    a_sizes = defaultdict(int)
    b_sizes = defaultdict(int)
    for m in shared_members:
        ra, rb = a_to_rep[m], b_to_rep[m]
        contingency[(ra, rb)] += 1
        a_sizes[ra] += 1
        b_sizes[rb] += 1

    n = len(shared_members)
    sum_comb_c = sum(choose2(v) for v in contingency.values())
    sum_comb_a = sum(choose2(v) for v in a_sizes.values())
    sum_comb_b = sum(choose2(v) for v in b_sizes.values())
    comb_n = choose2(n)

    expected_index = (sum_comb_a * sum_comb_b) / comb_n if comb_n else 0
    max_index = 0.5 * (sum_comb_a + sum_comb_b)
    denom = max_index - expected_index
    if denom == 0:
        return 1.0 if sum_comb_c == expected_index else 0.0
    return (sum_comb_c - expected_index) / denom


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", required=True, metavar="LABEL=PATH")
    ap.add_argument("--b", required=True, metavar="LABEL=PATH")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    a_label, a_path = args.a.split("=", 1)
    b_label, b_path = args.b.split("=", 1)

    a_to_rep = load_membership(a_path)
    b_to_rep = load_membership(b_path)

    shared_members = set(a_to_rep) & set(b_to_rep)
    only_a = set(a_to_rep) - shared_members
    only_b = set(b_to_rep) - shared_members

    ari = adjusted_rand_index(a_to_rep, b_to_rep, shared_members)

    a_stats = cluster_stats(a_to_rep)
    b_stats = cluster_stats(b_to_rep)

    with open(args.out, "w") as fh:
        fh.write(f"# clustering comparison: {a_label} vs {b_label}\n")
        fh.write(f"shared_members\t{len(shared_members)}\n")
        fh.write(f"members_only_in_{a_label}\t{len(only_a)}\n")
        fh.write(f"members_only_in_{b_label}\t{len(only_b)}\n")
        fh.write(f"adjusted_rand_index\t{ari:.4f}\n")
        fh.write("\nmetric\t{}\t{}\n".format(a_label, b_label))
        for key in ("n_members", "n_clusters", "n_singletons", "mean_size", "median_size", "max_size"):
            fh.write(f"{key}\t{a_stats[key]}\t{b_stats[key]}\n")

    print(
        f"compare_clusterings: {a_label} n_clusters={a_stats['n_clusters']} "
        f"vs {b_label} n_clusters={b_stats['n_clusters']}, "
        f"ARI={ari:.4f} (1.0=identical clustering, 0.0=chance-level agreement)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
