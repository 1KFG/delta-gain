from build_cluster_manifest import build_cluster_manifest

def test_manifest_records_screen_out_reasons_separately():
    manifest = build_cluster_manifest(
        run_date="2026-09-09",
        n_genomes_total=23683,
        n_genomes_annotated=10922,
        mmseqs_version="15.6f452",
        mmseqs_params={"min_seq_id": 0.95, "c": 0.9, "cov_mode": 0},
        qc_thresholds={"busco_complete_pct_min": 70.0, "n50_bp_min": 10000},
        screened_out={"busco_complete_pct_min": ["ASM1", "ASM2"], "n50_bp_min": ["ASM3"]},
        zero_protein_asmids=["ASM4"],
    )
    assert manifest["n_genomes_annotated"] == 10922
    assert manifest["screened_out"]["busco_complete_pct_min"] == ["ASM1", "ASM2"]
    assert manifest["zero_protein_asmids"] == ["ASM4"]
    assert "version_tag" in manifest
    assert manifest["version_tag"].startswith("bfd-v20260909-")
