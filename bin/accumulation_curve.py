#!/usr/bin/env python3
"""CLI wrapping lib_accumulation.py: writes the accumulation curve summary,
pan-proteome power-law fit, per-clade Shapley contribution, and raw
per-permutation marginal records (design doc: Components, Pipeline
integration).

After the outputs are written it runs validation #4 (positive-control
cross-check with Stage 1) as an automated assertion, per the design doc's
explicit instruction to implement it here rather than only as a standalone
script; a violation exits nonzero."""
import argparse
import csv
import sys

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from lib_accumulation import (
    run_all_permutations,
    aggregate_clade_contribution,
    fit_power_law_per_permutation,
)
from validate_positive_control import check_positive_control


def _load_inputs(incidence_path, external_status_path, genome_metadata_path):
    incidence = pq.read_table(incidence_path).to_pylist()
    asmid_clusters = {}
    for row in incidence:
        asmid_clusters.setdefault(row["asmid"], set()).add(row["cluster_id"])

    ext_status_table = pq.read_table(external_status_path).to_pylist()
    external_status = {r["cluster_id"]: r["external_status"] for r in ext_status_table}

    genome_metadata = pq.read_table(genome_metadata_path).to_pylist()
    clade_of = {r["asmid"]: (r["clade_rank"], r["clade_label"]) for r in genome_metadata}

    # incidence_matrix.parquet and genome_metadata.parquet are independent
    # Nextflow path inputs; nothing guarantees they were staged from the same
    # run. A genome present in one but not the other would otherwise be
    # silently dropped from clade_contribution.tsv while still counted in the
    # curves -- a silent under-report, not a crash. Fail loudly instead (same
    # pattern as check_qc_confounds.py's missing-asmid guard, learning L-8).
    only_in_incidence = sorted(set(asmid_clusters) - set(clade_of))
    only_in_metadata = sorted(set(clade_of) - set(asmid_clusters))
    if only_in_incidence or only_in_metadata:
        raise ValueError(
            "genome universes disagree between incidence matrix and genome "
            f"metadata: in incidence matrix only={only_in_incidence}; "
            f"in genome metadata only={only_in_metadata}"
        )

    return asmid_clusters, external_status, clade_of


def _write_summary_tsv(walks, path, version_tag):
    n_genomes = len(walks[0])
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["N", "internal_mean", "internal_p2_5", "internal_p97_5",
                          "external_mean", "external_p2_5", "external_p97_5",
                          "version_tag"])
        for position in range(1, n_genomes + 1):
            internal_vals = np.array([w[position - 1]["internal_cumulative"] for w in walks], dtype=float)
            external_vals = np.array([w[position - 1]["external_cumulative"] for w in walks], dtype=float)
            writer.writerow([
                position,
                internal_vals.mean(), np.percentile(internal_vals, 2.5), np.percentile(internal_vals, 97.5),
                external_vals.mean(), np.percentile(external_vals, 2.5), np.percentile(external_vals, 97.5),
                version_tag,
            ])


def _write_fit_tsv(walks, path, version_tag):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["series", "alpha_mean", "alpha_ci_low", "alpha_ci_high",
                          "kappa_mean", "r_squared_mean", "version_tag"])
        for series, count_key in [("internal", "new_internal_count"), ("external", "new_external_count")]:
            fit = fit_power_law_per_permutation(walks, count_key)
            writer.writerow([series, fit["alpha_mean"], fit["alpha_ci_low"], fit["alpha_ci_high"],
                              fit["kappa_mean"], fit["r_squared_mean"], version_tag])


def _write_clade_tsv(walks, clade_of, path, version_tag):
    rows = aggregate_clade_contribution(walks, clade_of)
    for row in rows:
        row["version_tag"] = version_tag
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _write_raw_marginals_parquet(walks, permutation_ids, path, version_tag):
    records = []
    for perm_id, walk in zip(permutation_ids, walks):
        for record in walk:
            records.append({"permutation_id": perm_id, **record, "version_tag": version_tag})
    table = pa.table({key: [r[key] for r in records] for key in records[0]})
    pq.write_table(table, path, compression="zstd")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--incidence-matrix", required=True)
    ap.add_argument("--external-status", required=True)
    ap.add_argument("--genome-metadata", required=True)
    ap.add_argument("--n-permutations", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--version-tag", required=True)
    ap.add_argument("--max-marginal-frac", type=float, default=0.05,
                    help="validation #4 threshold (redundant_strain_max_marginal_frac)")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    asmid_clusters, external_status, clade_of = _load_inputs(
        args.incidence_matrix, args.external_status, args.genome_metadata
    )
    asmids = sorted(asmid_clusters)
    walks = run_all_permutations(asmids, asmid_clusters, external_status,
                                  args.n_permutations, args.seed)
    permutation_ids = list(range(args.n_permutations))

    _write_summary_tsv(walks, f"{args.outdir}/accumulation_curve_summary.tsv", args.version_tag)
    _write_fit_tsv(walks, f"{args.outdir}/pangenome_powerlaw_fit.tsv", args.version_tag)
    _write_clade_tsv(walks, clade_of, f"{args.outdir}/clade_contribution.tsv", args.version_tag)
    raw_marginals_path = f"{args.outdir}/permutation_marginals.parquet"
    _write_raw_marginals_parquet(walks, permutation_ids, raw_marginals_path, args.version_tag)

    # Validation #4 (design doc: "Implement as an automated assertion in
    # bin/accumulation_curve.py's output validation"). A violation means a
    # real bug in the ANI calls or the incidence matrix, so fail the process
    # rather than publishing outputs that look fine.
    result = check_positive_control(raw_marginals_path, args.genome_metadata,
                                     args.max_marginal_frac)
    if not result["passed"]:
        for v in result["violations"]:
            print(f"VIOLATION: {v['asmid']} conditioned median marginal "
                  f"{v['median_marginal']}/{v['n_proteins']} proteins "
                  f"({v['frac']:.1%}) >= threshold {args.max_marginal_frac:.1%}",
                  file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
