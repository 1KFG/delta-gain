#!/usr/bin/env python3
"""Validation #4: redundant_strain genomes must show near-zero marginal
contribution, CONDITIONED on at least one same-ANI-component genome already
having been inserted earlier in that permutation (design doc: "Positive-
control cross-check with Stage 1", docs/superpowers/specs/
2026-09-09-protein-novelty-accumulation-curve-design.md). An unconditioned
average is wrong: in roughly half of permutations a redundant_strain genome
is inserted before its ANI partner and gets full first-discovery credit for
their shared content, so its unconditional average will not be near zero
even when the pipeline is correct.

Component representative, not every redundant_strain member independently:
`classify_genomes.py` labels a genome redundant_strain based on its OWN best
within-component ANI hit, so two mutual best-ANI-partners in a size-2
component can BOTH end up labeled redundant_strain. Testing both of them
independently double-counts a single underlying ANI relationship: in
permutations where genome X is inserted second, X's conditioned marginal
reflects the pair's shared content; in permutations where the roles reverse,
the SAME shared-content relationship shows up again, now as the partner's
conditioned marginal. These are two views of one relationship, not two
independent pieces of evidence. So for each ANI component we designate one
deterministic representative (the alphabetically-first redundant_strain
asmid) as the tested genome; every other member of that component -- of any
cluster_class -- still counts as a valid "inserted earlier" anchor, it is
just never independently tested itself. For the common case of a component
with a single redundant_strain member, this is a no-op: that lone member is
its own representative.
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

    component_of = {asmid: meta["component_id"] for asmid, meta in genome_meta.items()}

    # One tested representative per component (see module docstring); other
    # same-component genomes are still used as position anchors below, they
    # are just excluded from `tested_asmids`.
    redundant_by_component = {}
    for asmid, meta in genome_meta.items():
        if meta["cluster_class"] == "redundant_strain":
            redundant_by_component.setdefault(meta["component_id"], []).append(asmid)
    tested_asmids = {sorted(asmids)[0] for asmids in redundant_by_component.values()}

    conditioned_marginals = {asmid: [] for asmid in tested_asmids}
    for perm_id, records in by_permutation.items():
        position_of = {r["asmid"]: r["position"] for r in records}
        for record in records:
            asmid = record["asmid"]
            if asmid not in tested_asmids:
                continue
            component = component_of[asmid]
            partner_positions = [
                position_of[other] for other, comp in component_of.items()
                if comp == component and other != asmid and other in position_of
            ]
            if partner_positions and min(partner_positions) < record["position"]:
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
