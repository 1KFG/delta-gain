import pytest
import pyarrow as pa
import pyarrow.parquet as pq
from check_qc_confounds import compute_qc_confound_correlations

def test_detects_perfect_correlation_with_proteome_size(tmp_path):
    # 3 genomes, marginal contribution engineered to scale exactly with
    # n_proteins -- correlation must come out ~1.0.
    raw = pa.table({
        "permutation_id": [0, 0, 0],
        "asmid": ["A", "B", "C"],
        "new_internal_count": [10, 20, 30],
        "new_external_count": [1, 2, 3],
    })
    pq.write_table(raw, tmp_path / "permutation_marginals.parquet")
    genome_meta = pa.table({
        "asmid": ["A", "B", "C"],
        "n_proteins": [100, 200, 300],
        "complete_pct": [95.0, 80.0, 60.0],
        "n50_bp": [500000, 100000, 20000],
    })
    pq.write_table(genome_meta, tmp_path / "genome_metadata.parquet")

    result = compute_qc_confound_correlations(
        str(tmp_path / "permutation_marginals.parquet"),
        str(tmp_path / "genome_metadata.parquet"),
    )
    assert result["internal_raw_vs_n_proteins"] > 0.99
    assert result["internal_raw_vs_complete_pct"] < -0.99  # engineered inversely here
    assert "internal_per_1000_vs_n_proteins" in result


def test_raises_on_missing_asmid_in_metadata(tmp_path):
    # Test that a clear ValueError is raised if an asmid appears in raw marginals
    # but is missing from genome metadata (e.g. genome dropped from QC but still in results)
    raw = pa.table({
        "permutation_id": [0, 0, 0],
        "asmid": ["A", "B", "C"],  # C is present in marginals
        "new_internal_count": [10, 20, 30],
        "new_external_count": [1, 2, 3],
    })
    pq.write_table(raw, tmp_path / "permutation_marginals.parquet")
    genome_meta = pa.table({
        "asmid": ["A", "B"],  # C is missing from metadata
        "n_proteins": [100, 200],
        "complete_pct": [95.0, 80.0],
        "n50_bp": [500000, 100000],
    })
    pq.write_table(genome_meta, tmp_path / "genome_metadata.parquet")

    with pytest.raises(ValueError, match="asmid\\(s\\) in raw marginals but missing from genome metadata"):
        compute_qc_confound_correlations(
            str(tmp_path / "permutation_marginals.parquet"),
            str(tmp_path / "genome_metadata.parquet"),
        )
