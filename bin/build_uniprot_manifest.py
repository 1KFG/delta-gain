#!/usr/bin/env python3
"""Build the uniprot_fungi_nr.manifest.json provenance sidecar (design doc:
UniProt reference versioning)."""
from version_tag import build_version_tag


def build_uniprot_manifest(release_info: dict, mmseqs_version: str,
                            mmseqs_params: dict, seq_count_in: int,
                            seq_count_out: int) -> dict:
    manifest = {
        "release": release_info["release"],
        "release_date": release_info["release_date"],
        "mmseqs_version": mmseqs_version,
        "mmseqs_params": mmseqs_params,
        "seq_count_in": seq_count_in,
        "seq_count_out": seq_count_out,
    }
    manifest["version_tag"] = build_version_tag(
        {"release": release_info["release"], "n": seq_count_out},
        prefix="uniprot-fungi",
        manifest=manifest,
    )
    return manifest


if __name__ == "__main__":
    import argparse
    import json
    from parse_uniprot_relnotes import parse_relnotes

    ap = argparse.ArgumentParser()
    ap.add_argument("--relnotes", required=True)
    ap.add_argument("--mmseqs-version", required=True)
    ap.add_argument("--min-seq-id", type=float, required=True)
    ap.add_argument("--cov", type=float, required=True)
    ap.add_argument("--cov-mode", type=int, required=True)
    ap.add_argument("--seq-count-in", type=int, required=True)
    ap.add_argument("--seq-count-out", type=int, required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    release_info = parse_relnotes(open(args.relnotes).read())
    manifest = build_uniprot_manifest(
        release_info, args.mmseqs_version,
        {"min_seq_id": args.min_seq_id, "c": args.cov, "cov_mode": args.cov_mode},
        args.seq_count_in, args.seq_count_out,
    )
    json.dump(manifest, open(args.out, "w"), indent=2)
