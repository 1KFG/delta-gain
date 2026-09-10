#!/usr/bin/env python3
"""Validation #4: redundant_strain genomes must show near-zero marginal
contribution, CONDITIONED on that genome's own recorded best_ani_partner
already having been inserted earlier in that permutation (design doc:
"Positive-control cross-check with Stage 1", docs/superpowers/specs/
2026-09-09-protein-novelty-accumulation-curve-design.md). An unconditioned
average is wrong: in roughly half of permutations a redundant_strain genome
is inserted before its ANI partner and gets full first-discovery credit for
their shared content, so its unconditional average will not be near zero
even when the pipeline is correct.

Conditioning on best_ani_partner specifically, NOT on "any genome in the
same ANI component" (fixed 2026-09-09, found running this for real on the
n=100 pilot data): a component from find_components.py's union-find graph
can be a loose chain of 10+ genomes connected transitively via mash/skani
distance, not a clique -- most members can be only just above the
clustering floor from each other even when the specific best_ani_partner
pair is 99.99% identical. Conditioning on "any component member" accepts a
distantly-related chain-mate as satisfying the condition, which lets a
genuinely near-duplicate pair's real relationship go untested while a
spurious, unrelated marginal count gets recorded instead -- exactly the
false-positive violation pattern seen on real pilot data before this fix
(30-90%+ "conditioned" marginals for pairs at 99.99% ANI).

Every redundant_strain genome is tested independently. Two mutual
best-ANI-partners can both be labeled redundant_strain (see
`classify_genomes.py`), and testing both is correct, not double-counting:
"is A redundant given B present" and "is B redundant given A present" are
two distinct, legitimate questions, and a genuinely, symmetrically redundant
pair should show near-zero conditioned marginal in BOTH directions.
"""
import argparse
import sys

import pyarrow.parquet as pq


def check_positive_control(raw_marginals_path, genome_metadata_path, max_marginal_frac=0.05):
    raw = pq.read_table(raw_marginals_path).to_pylist()
    genome_meta = {r["asmid"]: r for r in pq.read_table(genome_metadata_path).to_pylist()}

    # Group raw records by permutation, so we know each genome's position
    # relative to its ANI-component-mates within that specific permutation.
    by_permutation = {}
    for record in raw:
        by_permutation.setdefault(record["permutation_id"], []).append(record)

    redundant_asmids = {
        asmid for asmid, meta in genome_meta.items()
        if meta["cluster_class"] == "redundant_strain"
    }
    partner_of = {asmid: meta.get("best_ani_partner") for asmid, meta in genome_meta.items()}

    conditioned_marginals = {asmid: [] for asmid in redundant_asmids}
    for perm_id, records in by_permutation.items():
        position_of = {r["asmid"]: r["position"] for r in records}
        for record in records:
            asmid = record["asmid"]
            if asmid not in redundant_asmids:
                continue
            partner = partner_of.get(asmid)
            if partner is None or partner not in position_of:
                continue  # partner not in this genome universe -- untestable, not a violation.
            if position_of[partner] < record["position"]:
                conditioned_marginals[asmid].append(record["new_internal_count"])

    violations = []
    for asmid, marginals in conditioned_marginals.items():
        if not marginals:
            continue  # no conditioned observations available -- not a violation, just untestable.
        median_marginal = sorted(marginals)[len(marginals) // 2]
        n_proteins = genome_meta[asmid]["n_proteins"]
        frac = median_marginal / n_proteins if n_proteins else 0
        if frac >= max_marginal_frac:
            violations.append({"asmid": asmid, "median_marginal": median_marginal,
                                "n_proteins": n_proteins, "frac": frac})

    return {"passed": len(violations) == 0, "violations": violations}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-marginals", required=True)
    ap.add_argument("--genome-metadata", required=True)
    ap.add_argument("--max-marginal-frac", type=float, default=0.05)
    args = ap.parse_args()

    result = check_positive_control(args.raw_marginals, args.genome_metadata, args.max_marginal_frac)
    if not result["passed"]:
        for v in result["violations"]:
            print(f"VIOLATION: {v['asmid']} conditioned median marginal "
                  f"{v['median_marginal']}/{v['n_proteins']} proteins "
                  f"({v['frac']:.1%}) >= threshold {args.max_marginal_frac:.1%}", file=sys.stderr)
        sys.exit(1)
    print("positive-control check passed")
