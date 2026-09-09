import pyarrow as pa
import pyarrow.parquet as pq

from validate_positive_control import check_positive_control


def test_redundant_strain_conditioned_on_ani_partner_first(tmp_path):
    # A and B are the same ANI component (comp_1); A is redundant_strain.
    # Permutation 0: order A,B (A first -> A gets full credit for shared
    # cluster c2 -- must be EXCLUDED from A's conditioned average).
    # Permutation 1: order B,A (B first -> A's marginal on c2 is 0, the
    # correct conditioned observation).
    raw = pa.table({
        "permutation_id": [0, 0, 1, 1],
        "asmid": ["A", "B", "B", "A"],
        "position": [1, 2, 1, 2],
        "new_internal_count": [2, 1, 2, 0],
        "new_external_count": [0, 0, 0, 0],
    })
    pq.write_table(raw, tmp_path / "permutation_marginals.parquet")
    genome_meta = pa.table({
        "asmid": ["A", "B", "C"],
        "component_id": ["comp_1", "comp_1", "comp_2"],
        "cluster_class": ["redundant_strain", "redundant_strain", "singleton_isolated"],
        "n_proteins": [3, 3, 1],
    })
    pq.write_table(genome_meta, tmp_path / "genome_metadata.parquet")

    result = check_positive_control(
        str(tmp_path / "permutation_marginals.parquet"),
        str(tmp_path / "genome_metadata.parquet"),
        max_marginal_frac=0.05,
    )
    # A's only conditioned observation (permutation 1, B-before-A) is 0 new
    # clusters -> 0% of A's 3 proteins -> well under 5% -> passes.
    assert result["passed"] is True
    assert result["violations"] == []


def test_detects_a_real_violation():
    # B (A's ANI-component partner) is inserted first (position 1), so this
    # IS a valid conditioned observation for A -- A's 3/3-protein marginal
    # at position 2 is a genuine violation, not a missing-data non-event.
    raw = pa.table({
        "permutation_id": [0, 0],
        "asmid": ["B", "A"],
        "position": [1, 2],
        "new_internal_count": [2, 3],  # 3/3 proteins "new" even though conditioned -> a bug signal
        "new_external_count": [0, 0],
    })
    import tempfile, os
    tmp = tempfile.mkdtemp()
    pq.write_table(raw, os.path.join(tmp, "permutation_marginals.parquet"))
    genome_meta = pa.table({
        "asmid": ["A", "B"],
        "component_id": ["comp_1", "comp_1"],
        "cluster_class": ["redundant_strain", "redundant_strain"],
        "n_proteins": [3, 3],
    })
    pq.write_table(genome_meta, os.path.join(tmp, "genome_metadata.parquet"))
    result = check_positive_control(
        os.path.join(tmp, "permutation_marginals.parquet"),
        os.path.join(tmp, "genome_metadata.parquet"),
        max_marginal_frac=0.05,
    )
    assert result["passed"] is False
    assert result["violations"][0]["asmid"] == "A"
