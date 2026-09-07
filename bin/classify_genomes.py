#!/usr/bin/env python3
"""Combine mash components + per-component skani triangle ANI + existing
labels into the Step 1 genome classification table (design doc axis 4 input).

This is BFD-internal classification only (see MASH_PREFILTER_GLOBAL's header
comment for why reference genomes are deliberately NOT mixed in here): it
answers "which BFD genomes are redundant strains of which other BFD genomes,"
and flags label disagreements as an ANI-vs-metadata QC signal. It does NOT by
itself decide whether a species is already known to UniProt/RefSeq -- that is
a separate taxonomy lookup (see docs/TAXONOMIC_BREADTH.md, not yet built),
independent of genome ANI.

cluster_class values:
  singleton_isolated      component of size 1 -- no other BFD genome within
                           the mash prefilter ANI floor. Candidate novel
                           lineage OR simply the only genome for a well-known
                           species; disambiguate via the taxonomy lookup.
  redundant_strain         genome ANI >= ani_cluster_threshold to another
                           genome of the SAME existing_species label in its
                           component.
  novel_strain_same_species genome is in a multi-member component whose
                           majority label matches its own label, but its best
                           within-component ANI falls between the outlier and
                           cluster thresholds (a more divergent strain).
  label_mismatch           genome's existing_species label disagrees with the
                           majority label of its own high-ANI component --
                           likely mislabeled, or a cryptic/novel species
                           sitting inside what looks like a known-species
                           cluster by name only.
  unresolved_below_ani_floor  best within-component ANI (skani) is below
                           params.ani_outlier_threshold despite mash having
                           grouped it into a multi-member component -- treat
                           as unresolved, defer to phylogenomic placement
                           (nf_phyling), not to a raw ANI number.
"""
import argparse
import csv
import glob
import os
from collections import Counter, defaultdict


def load_components(path):
    rows = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rows[row["asmid"]] = row
    return rows


def _asmid_from_path(path_field):
    # skani echoes the path it was given; genome files are <ASMID>.fa.gz.
    base = os.path.basename(path_field)
    return base.split(".fa")[0]


def load_skani_dir(skani_dir):
    """`skani triangle --sparse` output: sparse (query_path, ref_path, ANI%%)
    triplets, one file per component (no header) -- same format as the
    existing Fungi_BFD ANI pipeline's *.ani.tsv (see
    Fungi_BFD_runs/results/ANI/skani/SPECIES/*/*.ani.tsv for a real example).
    Returns {comp_id: {(asmid_a, asmid_b): ani_pct}} (both directions)."""
    best_ani = defaultdict(dict)
    for path in glob.glob(os.path.join(skani_dir, "comp_*.skani.tsv")):
        comp_id = os.path.basename(path).split(".")[0]
        if not os.path.getsize(path):
            continue
        with open(path) as fh:
            for line in fh:
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                a, b = _asmid_from_path(parts[0]), _asmid_from_path(parts[1])
                try:
                    ani = float(parts[2])
                except ValueError:
                    continue
                best_ani[comp_id][(a, b)] = ani
                best_ani[comp_id][(b, a)] = ani
    return best_ani


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--components", required=True)
    ap.add_argument("--skani-dir", required=True, help="dir containing comp_*.skani.tsv")
    ap.add_argument("--ani-cluster-threshold", type=float, default=95.0)
    ap.add_argument("--ani-outlier-threshold", type=float, default=90.0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    components = load_components(args.components)
    skani = load_skani_dir(args.skani_dir)

    # Majority label per component (for label_mismatch detection).
    comp_labels = defaultdict(Counter)
    for asmid, row in components.items():
        comp_labels[row["component_id"]][row["existing_species"]] += 1

    with open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow([
            "asmid", "component_id", "component_size", "existing_species",
            "majority_species_in_component", "best_within_component_ani",
            "best_ani_partner", "cluster_class",
        ])
        for asmid, row in sorted(components.items()):
            comp_id = row["component_id"]
            size = int(row["component_size"])
            species = row["existing_species"]
            majority_species = comp_labels[comp_id].most_common(1)[0][0]

            if size == 1:
                w.writerow([asmid, comp_id, size, species, majority_species,
                            "", "", "singleton_isolated"])
                continue

            partner_anis = {
                partner: ani
                for (a, partner), ani in skani.get(comp_id, {}).items()
                if a == asmid
            }
            if not partner_anis:
                w.writerow([asmid, comp_id, size, species, majority_species,
                            "", "", "unresolved_below_ani_floor"])
                continue

            best_partner = max(partner_anis, key=partner_anis.get)
            best_ani = partner_anis[best_partner]

            if best_ani < args.ani_outlier_threshold:
                cls = "unresolved_below_ani_floor"
            elif species != majority_species:
                cls = "label_mismatch"
            elif best_ani >= args.ani_cluster_threshold:
                cls = "redundant_strain"
            else:
                cls = "novel_strain_same_species"

            w.writerow([asmid, comp_id, size, species, majority_species,
                        f"{best_ani:.2f}", best_partner, cls])


if __name__ == "__main__":
    main()
