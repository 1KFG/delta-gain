#!/usr/bin/env python3
"""Validation #3: correlate per-genome marginal new-cluster contribution
against proteome size and BUSCO/N50 QC metrics, before attributing a
clade's high contribution to real phylogenetic novelty (design doc:
Confound check)."""
import argparse

import numpy as np
import pyarrow.parquet as pq


def compute_qc_confound_correlations(raw_marginals_path, genome_metadata_path):
    raw = pq.read_table(raw_marginals_path).to_pylist()
    genome_meta = {r["asmid"]: r for r in pq.read_table(genome_metadata_path).to_pylist()}

    mean_marginal = {}
    for record in raw:
        asmid = record["asmid"]
        acc = mean_marginal.setdefault(asmid, {"internal": [], "external": []})
        acc["internal"].append(record["new_internal_count"])
        acc["external"].append(record["new_external_count"])

    asmids = sorted(mean_marginal)

    # Check for asmids in raw marginals but missing from genome metadata
    missing = set(asmids) - set(genome_meta.keys())
    if missing:
        raise ValueError(
            f"asmid(s) in raw marginals but missing from genome metadata: {sorted(missing)}"
        )

    internal_raw = np.array([np.mean(mean_marginal[a]["internal"]) for a in asmids])
    n_proteins = np.array([genome_meta[a]["n_proteins"] for a in asmids], dtype=float)
    complete_pct = np.array([genome_meta[a]["complete_pct"] for a in asmids], dtype=float)
    n50_bp = np.array([genome_meta[a]["n50_bp"] for a in asmids], dtype=float)
    internal_per_1000 = internal_raw / n_proteins * 1000

    def _corr(x, y):
        return float(np.corrcoef(x, y)[0, 1])

    return {
        "internal_raw_vs_n_proteins": _corr(internal_raw, n_proteins),
        "internal_raw_vs_complete_pct": _corr(internal_raw, complete_pct),
        "internal_raw_vs_n50_bp": _corr(internal_raw, n50_bp),
        "internal_per_1000_vs_n_proteins": _corr(internal_per_1000, n_proteins),
        "internal_per_1000_vs_complete_pct": _corr(internal_per_1000, complete_pct),
        "internal_per_1000_vs_n50_bp": _corr(internal_per_1000, n50_bp),
    }


if __name__ == "__main__":
    import json

    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-marginals", required=True)
    ap.add_argument("--genome-metadata", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    result = compute_qc_confound_correlations(args.raw_marginals, args.genome_metadata)
    json.dump(result, open(args.out, "w"), indent=2)
