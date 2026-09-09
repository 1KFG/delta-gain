from build_uniprot_manifest import build_uniprot_manifest


def test_manifest_has_required_fields_and_tag():
    release_info = {"release": "2026_03", "release_date": "2026-09-02"}
    manifest = build_uniprot_manifest(
        release_info=release_info,
        mmseqs_version="15.6f452",
        mmseqs_params={"min_seq_id": 0.95, "c": 0.9, "cov_mode": 0},
        seq_count_in=17_500_000,
        seq_count_out=15_089_757,
    )
    assert manifest["release"] == "2026_03"
    assert manifest["seq_count_out"] == 15_089_757
    assert "version_tag" in manifest
    assert manifest["version_tag"].startswith("uniprot-fungi-v2026_03-")
