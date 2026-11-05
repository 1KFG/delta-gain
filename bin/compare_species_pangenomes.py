#!/usr/bin/env python3
"""Cross-species comparison of pangenomes (follows bin/species_pangenome.py).

Inputs: the incidence matrix, samples.csv (taxonomy), and the table written
by species_pangenome.py (so the species set, --min-genomes filter and
--exclude-asmids choice are identical in both steps).

Outputs:
  species_cluster_sharing.tsv  one row per species pair that shares >= 1
      cluster: shared count, Jaccard of the two species' full cluster sets,
      and whether the pair is in the same genus / family.
  taxon_pangenome_summary.tsv  for each --ranks level (default GENUS, FAMILY,
      ORDER, CLASS, PHYLUM): n species, median alpha, median core fraction,
      median rarefied pangenome size, share of species with alpha < 1 (open).

A cluster present in more than --max-species-per-cluster species is too broad
to say anything about relatedness and would make pair counting quadratic. It
is left out of the pair counts and counted per species in n_broad_clusters
(reported in the sharing table's companion file species_broad_clusters.tsv),
so the cap is visible rather than silent.
"""
import argparse
import csv
import itertools
from statistics import median

from species_pangenome import load_incidence, load_exclude

RANKS = ["GENUS", "FAMILY", "ORDER", "CLASS", "PHYLUM"]


def load_taxonomy(samples_csv):
    tax, species_of = {}, {}
    with open(samples_csv) as fh:
        for row in csv.DictReader(fh):
            species_of[row["ASMID"]] = (row.get("SPECIES") or "").strip()
            tax[row["ASMID"]] = {r: (row.get(r) or "").strip() for r in RANKS}
    return species_of, tax


def species_cluster_sets(asmid_clusters, species_of, keep_species, exclude):
    sets = {}
    for asmid, clusters in asmid_clusters.items():
        sp = species_of.get(asmid, "")
        if asmid in exclude or sp not in keep_species:
            continue
        sets.setdefault(sp, set()).update(clusters)
    return sets


def sharing_pairs(sets, max_species):
    by_cluster = {}
    for sp, clusters in sets.items():
        for c in clusters:
            by_cluster.setdefault(c, []).append(sp)
    shared, broad = {}, {sp: 0 for sp in sets}
    for c, sps in by_cluster.items():
        if len(sps) > max_species:
            for sp in sps:
                broad[sp] += 1
            continue
        for a, b in itertools.combinations(sorted(sps), 2):
            shared[(a, b)] = shared.get((a, b), 0) + 1
    return shared, broad


def lineage_of_species(species_of, tax):
    lin = {}
    for asmid, sp in species_of.items():
        if sp and sp not in lin:
            lin[sp] = tax[asmid]
    return lin


def taxon_summary(rows, lineage, ranks):
    out = []
    for rank in ranks:
        groups = {}
        for r in rows:
            label = lineage.get(r["species"], {}).get(rank, "")
            if label:
                groups.setdefault(label, []).append(r)
        for label, rs in sorted(groups.items()):
            alphas = [float(r["alpha_mean"]) for r in rs if r["alpha_mean"] != ""]
            rare = [float(r["rarefied_pan_size"]) for r in rs if r["rarefied_pan_size"] != ""]
            out.append({
                "rank": rank, "taxon": label, "n_species": len(rs),
                "median_alpha": median(alphas) if alphas else "",
                "frac_species_open_alpha_lt_1": (sum(a < 1 for a in alphas) / len(alphas)) if alphas else "",
                "median_core_fraction": median(float(r["core_fraction"]) for r in rs),
                "median_rarefied_pan_size": median(rare) if rare else "",
            })
    return out


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--incidence-matrix", required=True)
    ap.add_argument("--samples-csv", required=True)
    ap.add_argument("--species-table", required=True, help="output of species_pangenome.py")
    ap.add_argument("--exclude-asmids", default="")
    ap.add_argument("--max-species-per-cluster", type=int, default=100)
    ap.add_argument("--ranks", default=",".join(RANKS))
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    with open(args.species_table) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    keep = {r["species"] for r in rows}
    species_of, tax = load_taxonomy(args.samples_csv)
    lineage = lineage_of_species(species_of, tax)
    sets = species_cluster_sets(load_incidence(args.incidence_matrix), species_of,
                                keep, load_exclude(args.exclude_asmids))
    missing = keep - set(sets)
    if missing:
        raise SystemExit(f"species in table but absent from incidence matrix: {sorted(missing)[:5]}")

    shared, broad = sharing_pairs(sets, args.max_species_per_cluster)
    pair_rows = []
    for (a, b), n in sorted(shared.items()):
        la, lb = lineage[a], lineage[b]
        pair_rows.append({
            "species_a": a, "species_b": b, "shared_clusters": n,
            "jaccard": n / (len(sets[a]) + len(sets[b]) - n),
            "same_genus": int(bool(la["GENUS"]) and la["GENUS"] == lb["GENUS"]),
            "same_family": int(bool(la["FAMILY"]) and la["FAMILY"] == lb["FAMILY"]),
        })
    write_tsv(f"{args.outdir}/species_cluster_sharing.tsv", pair_rows,
              ["species_a", "species_b", "shared_clusters", "jaccard", "same_genus", "same_family"])
    write_tsv(f"{args.outdir}/species_broad_clusters.tsv",
              [{"species": sp, "pan_size": len(sets[sp]), "n_broad_clusters": n}
               for sp, n in sorted(broad.items())],
              ["species", "pan_size", "n_broad_clusters"])
    ts = taxon_summary(rows, lineage, [r for r in args.ranks.split(",") if r])
    write_tsv(f"{args.outdir}/taxon_pangenome_summary.tsv", ts,
              ["rank", "taxon", "n_species", "median_alpha", "frac_species_open_alpha_lt_1",
               "median_core_fraction", "median_rarefied_pan_size"])
    print(f"species: {len(sets)}; pairs sharing >=1 cluster: {len(pair_rows)}; "
          f"taxon rows: {len(ts)}")


if __name__ == "__main__":
    main()
