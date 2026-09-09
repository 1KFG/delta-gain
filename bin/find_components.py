#!/usr/bin/env python3
"""Union-find connected components from a global all-vs-all mash distance table.

This is the label-free clustering step: edges are drawn purely from mash
distance (approximate ANI), with NO reference to any pre-existing species/genus
name. A component that mixes genomes from more than one taxonomic label is
exactly the kind of disagreement DeltaGain wants surfaced (Step 1 of the
design doc) -- it is either a mislabeled genome or a genuinely novel lineage
sitting closer to an unexpected neighbor than to its nominal species.

Input:
  --mash-dist   `mash dist` output (query, ref, distance, p-value, shared-hashes),
                no header, all-vs-all (produced by MASH_PREFILTER_GLOBAL).
  --samples     Fungi_BFD_runs/samples.csv (columns: ASMID, SPECIES_IN, STRAIN,
                BIOPROJECT, NCBI_TAXONID, BUSCO_LINEAGE, PHYLUM..SPECIES, ...) --
                reused as-is from the annotation pipeline, not reinvented here.
                Used only to attach existing labels to the output for downstream
                QC, never to gate the clustering.
  --genome-dir  directory holding <ASMID>.fa.gz (Fungi_BFD_runs/input_clean_genomes
                convention) -- used to resolve genome paths for the per-component
                genome-list files fed to SKANI_TRIANGLE.
  --prefilter-ani  ANI%% floor for drawing an edge (default 80.0 -- below this,
                mash distance itself is not reliable; genomes with no edge above
                this threshold are left as singleton components and flagged for
                the phylogenomic fallback rather than forced into a cluster).

Output (--outdir):
  components.tsv               asmid, component_id, component_size, existing_species
  component_<id>.genomes.tsv    asmid, genome_path   (one per multi-member component;
                                 fed to SKANI_TRIANGLE for exact-ANI refinement)
"""
import argparse
import csv
import os
import sys
from collections import defaultdict


class UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self, x):
        self.parent.setdefault(x, x)
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


def asmid_from_mash_id(mash_path):
    """Recover the bare ASMID from a `mash dist` query/ref field.

    Real bug found 2026-09-08 running the first real multi-genome pilot
    (never surfaced at the earlier --n_test 2 validation, which produced
    only singleton components and so never actually exercised
    SKANI_TRIANGLE): `mash sketch`'s internal sketch ID is the INPUT FASTA
    FILENAME given on its command line (here, the staged `<ASMID>.fa.gz`
    Nextflow input), NOT the `-o` output-prefix argument -- so `mash dist`
    reports query/ref as `<ASMID>.fa.gz`, not `<ASMID>.msh`. Stripping only
    a literal `.msh` suffix left the `.fa.gz` in place, so the genome path
    this script later builds (`<genome_dir>/<asmid>.fa.gz`) came out as
    `.../<ASMID>.fa.gz.fa.gz` -- skani then skipped every genome in every
    multi-member component ("not a valid fasta/fastq file"). Stripping both
    possible suffixes (order-independent, since only one will ever be
    present) is robust to either naming behavior.
    """
    name = os.path.basename(mash_path)
    if name.endswith(".msh"):
        name = name[: -len(".msh")]
    if name.endswith(".fa.gz"):
        name = name[: -len(".fa.gz")]
    return name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mash-dist", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--genome-dir", required=True)
    ap.add_argument("--prefilter-ani", type=float, default=80.0)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    samples = {}
    with open(args.samples, newline="") as fh:
        for row in csv.DictReader(fh):
            asmid = row["ASMID"].strip()
            samples[asmid] = row

    uf = UnionFind()
    # Seeded ONLY from genomes that actually appear in the mash-dist input
    # (i.e. actually mash-sketched this run), NOT from every row of
    # --samples. Real bug found 2026-09-08 running the first n_test-limited
    # pilot (--n_test 100): seeding from the full samples.csv (23,683 rows)
    # made every genome NOT included in this run's subset show up as a
    # `singleton_isolated` component -- falsely implying "compared, found no
    # neighbor above the ANI floor" for ~23,600 genomes that were never
    # mash-sketched at all in an n_test-limited run. `singleton_isolated`
    # must mean "actually compared, no match," not "not part of this run."
    all_asmids = set()

    with open(args.mash_dist, newline="") as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            query_path, ref_path, distance = parts[0], parts[1], parts[2]
            query = asmid_from_mash_id(query_path)
            ref = asmid_from_mash_id(ref_path)
            all_asmids.add(query)
            all_asmids.add(ref)
            if query == ref:
                continue
            try:
                ani_pct = (1.0 - float(distance)) * 100.0
            except ValueError:
                continue
            if ani_pct >= args.prefilter_ani:
                uf.union(query, ref)

    # Assign stable, size-ordered component ids.
    members = defaultdict(list)
    for asmid in sorted(all_asmids):
        members[uf.find(asmid)].append(asmid)

    ordered_roots = sorted(members, key=lambda r: (-len(members[r]), r))
    root_to_id = {root: f"comp_{i+1:05d}" for i, root in enumerate(ordered_roots)}

    components_path = os.path.join(args.outdir, "components.tsv")
    with open(components_path, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["asmid", "component_id", "component_size", "existing_species"])
        for root in ordered_roots:
            comp_id = root_to_id[root]
            asmids = members[root]
            for asmid in asmids:
                species = samples.get(asmid, {}).get("SPECIES", "")
                w.writerow([asmid, comp_id, len(asmids), species])

    n_multi = 0
    for root in ordered_roots:
        asmids = members[root]
        if len(asmids) < 2:
            continue
        n_multi += 1
        comp_id = root_to_id[root]
        genome_list_path = os.path.join(args.outdir, f"{comp_id}.genomes.tsv")
        with open(genome_list_path, "w", newline="") as out:
            w = csv.writer(out, delimiter="\t")
            for asmid in asmids:
                # Fungi_BFD_runs/input_clean_genomes convention: <ASMID>.fa.gz
                genome_path = os.path.join(args.genome_dir, f"{asmid}.fa.gz")
                w.writerow([asmid, genome_path])

    n_singleton = len(ordered_roots) - n_multi
    print(
        f"components: {len(ordered_roots)} total "
        f"({n_multi} multi-member -> skani refinement, "
        f"{n_singleton} singleton -> below mash prefilter ANI floor, "
        f"flagged for phylogenomic fallback)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
