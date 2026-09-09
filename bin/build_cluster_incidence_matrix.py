#!/usr/bin/env python3
"""Build the genome x cluster incidence matrix, cluster external-status
lookup, and per-genome metadata table (design doc: Components,
BUILD_CLUSTER_INCIDENCE_MATRIX; Genome universe and edge cases).

The permutation universe for bin/accumulation_curve.py is seeded from
bfd_proteins_provenance.tsv, never samples.csv directly -- a genome absent
from the real clustering input must never silently appear as if it
contributed zero (the same class of bug fixed once already in
find_components.py for the Stage 1 mash/skani pipeline, see learning L-8).
"""
import argparse
import csv

import duckdb
import pandas as pd

RANK_ORDER = ["FAMILY", "ORDER", "CLASS", "SUBCLASS", "SUBPHYLUM", "PHYLUM"]


def _read_clusters(path):
    """MMseqs2 createtsv format: repr_id<TAB>member_id, no header."""
    repr_of = {}
    with open(path) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            repr_id, member_id = row
            repr_of[member_id] = repr_id
    return repr_of


def _read_provenance(path):
    asmid_of = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        seen = set()
        for row in reader:
            new_id = row["new_id"]
            if new_id in seen:
                raise ValueError(f"duplicate protein id in provenance: {new_id}")
            seen.add(new_id)
            asmid_of[new_id] = row["asmid"]
    return asmid_of


def _read_external_status(path):
    status_of = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            status_of[row["query_id"]] = row["overall_status"]
    return status_of


def _read_taxonomy(path):
    tax = {}
    seen = set()
    with open(path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            asmid = row["ASMID"]
            if asmid in seen:
                raise ValueError(f"duplicate ASMID in samples.csv: {asmid}")
            seen.add(asmid)
            tax[asmid] = {
                rank: (row[rank] or None)
                for rank in ["PHYLUM", "SUBPHYLUM", "CLASS", "SUBCLASS",
                             "ORDER", "FAMILY", "GENUS", "SPECIES"]
            }
    return tax


def _read_genome_classification(path):
    info = {}
    seen = set()
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            asmid = row["asmid"]
            if asmid in seen:
                raise ValueError(f"duplicate ASMID in genome_classification.tsv: {asmid}")
            seen.add(asmid)
            info[asmid] = {
                "component_id": row["component_id"],
                "cluster_class": row["cluster_class"],
            }
    return info


def _assign_clade(asmid, tax_by_asmid, rank_counts, min_clade_n):
    """Density-adaptive clade assignment with an explicit roll-up label.

    A genome is tagged at the finest rank whose label has at least
    min_clade_n annotated genomes. When that means rolling UP past a finer
    rank that did have a value (but was too sparse to report on its own),
    the label is returned as e.g. "Agaricomycetes (other orders)" rather
    than the bare "Agaricomycetes" -- design doc, Components: a roll-up bar
    in the per-clade plot must never be misread as "the whole class/order"
    when it only represents the finer clades too sparse to report
    individually. A naturally-blank finer rank is NOT a roll-up: if the
    assigned rank is the genome's finest rank that had any value at all,
    the bare label is returned unchanged.
    """
    row = tax_by_asmid.get(asmid, {})
    finest_populated_failed = None  # finest rank that HAD a value but was too sparse
    for rank in RANK_ORDER:
        label = row.get(rank)
        if not label:
            continue
        if rank_counts[rank].get(label, 0) >= min_clade_n:
            if finest_populated_failed is None:
                return rank, label
            return rank, f"{label} (other {finest_populated_failed.lower()}s)"
        if finest_populated_failed is None:
            finest_populated_failed = rank
    phylum = row.get("PHYLUM") or "unclassified"
    if (finest_populated_failed is None or finest_populated_failed == "PHYLUM"
            or phylum == "unclassified"):
        return "PHYLUM", phylum
    return "PHYLUM", f"{phylum} (other {finest_populated_failed.lower()}s)"


def build_incidence_matrix(clusters_tsv, provenance_tsv, novelty_bins_tsv,
                            samples_csv, busco_parquet, asm_stats_parquet,
                            genome_classification_tsv, min_clade_n):
    repr_of = _read_clusters(clusters_tsv)
    asmid_of = _read_provenance(provenance_tsv)
    ext_status_by_protein = _read_external_status(novelty_bins_tsv)
    tax_by_asmid = _read_taxonomy(samples_csv)
    class_by_asmid = _read_genome_classification(genome_classification_tsv)

    # Genome universe = ASMIDs actually present in provenance (never samples.csv).
    pairs = set()
    n_proteins_by_asmid = {}
    cluster_reps_seen = set()
    for member_id, repr_id in repr_of.items():
        asmid = asmid_of.get(member_id)
        if asmid is None:
            continue  # protein not in provenance -- excluded, not a silent zero.
        pairs.add((asmid, repr_id))
        cluster_reps_seen.add(repr_id)
        n_proteins_by_asmid[asmid] = n_proteins_by_asmid.get(asmid, 0) + 1

    incidence = pd.DataFrame({
        "asmid": [p[0] for p in pairs],
        "cluster_id": [p[1] for p in pairs],
    })

    # Fail loudly rather than defaulting to "no_hit": bin_novelty_hits.py
    # emits a row for EVERY query in its input FASTA, no-hit ones included,
    # so a cluster representative missing here means the three inputs came
    # from mismatched runs -- a real data-mismatch bug. Defaulting would
    # silently inflate the external-novelty curve (same class of silent
    # false-negative as learning L-8).
    sorted_reps = sorted(cluster_reps_seen)
    missing_status = [c for c in sorted_reps if c not in ext_status_by_protein]
    if missing_status:
        raise ValueError(
            "cluster representative(s) present in the clustering but missing "
            f"from {novelty_bins_tsv}: {missing_status}"
        )

    external_status = pd.DataFrame({
        "cluster_id": sorted_reps,
        "external_status": [ext_status_by_protein[c] for c in sorted_reps],
    })

    # Per-clade genome counts at every rank, for the density-adaptive rollup.
    rank_counts = {rank: {} for rank in RANK_ORDER}
    all_asmids = sorted(n_proteins_by_asmid)
    for asmid in all_asmids:
        row = tax_by_asmid.get(asmid, {})
        for rank in RANK_ORDER:
            label = row.get(rank)
            if label:
                rank_counts[rank][label] = rank_counts[rank].get(label, 0) + 1

    con = duckdb.connect()
    busco = con.execute(f"SELECT ASMID, complete_pct FROM read_parquet('{busco_parquet}')").fetchall()
    asm_stats = con.execute(f"SELECT ASMID, N50_bp FROM read_parquet('{asm_stats_parquet}')").fetchall()
    complete_pct_by_asmid = dict(busco)
    n50_by_asmid = dict(asm_stats)

    genome_meta_rows = []
    for asmid in all_asmids:
        clade_rank, clade_label = _assign_clade(asmid, tax_by_asmid, rank_counts, min_clade_n)
        cls_info = class_by_asmid.get(asmid, {})
        tax_row = tax_by_asmid.get(asmid, {})
        genome_meta_rows.append({
            "asmid": asmid,
            "n_proteins": n_proteins_by_asmid[asmid],
            "complete_pct": complete_pct_by_asmid.get(asmid),
            "n50_bp": n50_by_asmid.get(asmid),
            "component_id": cls_info.get("component_id"),
            "cluster_class": cls_info.get("cluster_class"),
            "clade_rank": clade_rank,
            "clade_label": clade_label,
            **tax_row,
        })
    genome_metadata = pd.DataFrame(genome_meta_rows)

    return incidence, external_status, genome_metadata


if __name__ == "__main__":
    import pyarrow as pa
    import pyarrow.parquet as pq

    ap = argparse.ArgumentParser()
    ap.add_argument("--clusters-tsv", required=True)
    ap.add_argument("--provenance-tsv", required=True)
    ap.add_argument("--novelty-bins-tsv", required=True)
    ap.add_argument("--samples-csv", required=True)
    ap.add_argument("--busco-parquet", required=True)
    ap.add_argument("--asm-stats-parquet", required=True)
    ap.add_argument("--genome-classification-tsv", required=True)
    ap.add_argument("--min-clade-n", type=int, default=10)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    incidence, ext_status, genome_meta = build_incidence_matrix(
        args.clusters_tsv, args.provenance_tsv, args.novelty_bins_tsv,
        args.samples_csv, args.busco_parquet, args.asm_stats_parquet,
        args.genome_classification_tsv, args.min_clade_n,
    )
    pq.write_table(pa.Table.from_pandas(incidence, preserve_index=False),
                   f"{args.outdir}/incidence_matrix.parquet", compression="zstd")
    pq.write_table(pa.Table.from_pandas(ext_status, preserve_index=False),
                   f"{args.outdir}/cluster_external_status.parquet", compression="zstd")
    pq.write_table(pa.Table.from_pandas(genome_meta, preserve_index=False),
                   f"{args.outdir}/genome_metadata.parquet", compression="zstd")
