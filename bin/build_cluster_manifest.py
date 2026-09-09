#!/usr/bin/env python3
"""Build the bfd_cluster_manifest.json provenance sidecar for one BFD
clustering run (design doc: BFD dataset versioning). Records genome count,
QC screen-outs (by reason, not just the final kept set), and clustering
params, so multiple clustering iterations with different parameter/QC
choices stay distinguishable."""
from version_tag import build_version_tag


def build_cluster_manifest(run_date: str, n_genomes_total: int,
                            n_genomes_annotated: int, mmseqs_version: str,
                            mmseqs_params: dict, qc_thresholds: dict,
                            screened_out: dict, zero_protein_asmids: list) -> dict:
    date_compact = run_date.replace("-", "")
    manifest = {
        "run_date": run_date,
        "n_genomes_total": n_genomes_total,
        "n_genomes_annotated": n_genomes_annotated,
        "mmseqs_version": mmseqs_version,
        "mmseqs_params": mmseqs_params,
        "qc_thresholds": qc_thresholds,
        "screened_out": screened_out,
        "zero_protein_asmids": zero_protein_asmids,
    }
    params_short = f"mmseqs{int(mmseqs_params['min_seq_id']*100)}c{int(mmseqs_params['c']*100)}"
    manifest["version_tag"] = build_version_tag(
        {"date": date_compact, "n": n_genomes_annotated, "params": params_short},
        prefix="bfd",
        manifest=manifest,
    )
    return manifest
