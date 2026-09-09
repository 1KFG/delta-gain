import subprocess
import sys
import pyarrow.parquet as pq
import csv


def test_cli_runs_end_to_end_and_writes_expected_outputs(tmp_path):
    # Build the incidence-matrix inputs using Task 5's builder directly,
    # reusing its fixtures, then feed them into the CLI under test.
    from build_cluster_incidence_matrix import build_incidence_matrix
    import pyarrow as pa
    import pyarrow.parquet as pq2
    from pathlib import Path

    fixtures = Path("tests/data/accumulation")
    from tests.test_build_cluster_incidence_matrix import _write_qc_fixtures
    busco_pq, asm_pq = _write_qc_fixtures(tmp_path)
    incidence, ext_status, genome_meta = build_incidence_matrix(
        clusters_tsv=fixtures / "bfd_proteins_clusters.tsv",
        provenance_tsv=fixtures / "bfd_proteins_provenance.tsv",
        novelty_bins_tsv=fixtures / "bfd_novelty_bins.tsv",
        samples_csv=fixtures / "samples.csv",
        busco_parquet=busco_pq, asm_stats_parquet=asm_pq,
        genome_classification_tsv=fixtures / "genome_classification.tsv",
        min_clade_n=2,
    )
    # Task 5's build_incidence_matrix returns pandas DataFrames; convert to
    # pyarrow Tables before writing Parquet, matching the pattern used in
    # build_cluster_incidence_matrix.py's own __main__ block.
    pq2.write_table(pa.Table.from_pandas(incidence, preserve_index=False),
                     tmp_path / "incidence_matrix.parquet")
    pq2.write_table(pa.Table.from_pandas(ext_status, preserve_index=False),
                     tmp_path / "cluster_external_status.parquet")
    pq2.write_table(pa.Table.from_pandas(genome_meta, preserve_index=False),
                     tmp_path / "genome_metadata.parquet")

    result = subprocess.run(
        [sys.executable, "bin/accumulation_curve.py",
         "--incidence-matrix", str(tmp_path / "incidence_matrix.parquet"),
         "--external-status", str(tmp_path / "cluster_external_status.parquet"),
         "--genome-metadata", str(tmp_path / "genome_metadata.parquet"),
         "--n-permutations", "20", "--seed", "1",
         "--version-tag", "test-v1-abc1234",
         # The toy fixture's two comp_1 "redundant_strain" genomes are not
         # actually redundant at the protein level (GENOME_B carries its own
         # GENOME_B__p2 cluster), so the real 0.05 validation-#4 threshold
         # trips on it by construction. Relax it here so this test keeps
         # covering the happy path; the threshold itself is exercised by
         # test_positive_control_violation_exits_nonzero below.
         "--max-marginal-frac", "1.1",
         "--outdir", str(tmp_path)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr

    # Task 5's fixtures now contain 4 genomes (GENOME_A..D, added during
    # Task 5's review fix round) with 5 distinct clusters total
    # (GENOME_A__p1, GENOME_A__p2, GENOME_B__p2, GENOME_C__p1, GENOME_D__p1),
    # of which 3 are "no_hit" (GENOME_A__p2, GENOME_C__p1, GENOME_D__p1).
    with open(tmp_path / "accumulation_curve_summary.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 4  # N = 1, 2, 3, 4 (four genomes total)
    assert rows[-1]["internal_mean"] == "5.0"  # all 5 clusters seen by N=4
    assert rows[-1]["external_mean"] == "3.0"  # all 3 no_hit clusters seen by N=4
    assert all(r["version_tag"] == "test-v1-abc1234" for r in rows)

    with open(tmp_path / "pangenome_powerlaw_fit.tsv") as fh:
        fit_rows = list(csv.DictReader(fh, delimiter="\t"))
    assert {r["series"] for r in fit_rows} == {"internal", "external"}
    assert all(r["version_tag"] == "test-v1-abc1234" for r in fit_rows)

    with open(tmp_path / "clade_contribution.tsv") as fh:
        clade_rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(clade_rows) >= 1
    assert all(r["version_tag"] == "test-v1-abc1234" for r in clade_rows)

    raw = pq.read_table(tmp_path / "permutation_marginals.parquet")
    assert raw.num_rows == 20 * 4  # n_permutations * n_genomes
    assert set(raw.column("version_tag").to_pylist()) == {"test-v1-abc1234"}


def _write_synthetic_inputs(tmp_path, incidence_rows, genome_rows):
    """Write the three Parquet inputs the CLI consumes, hand-crafted rather
    than derived from build_incidence_matrix, so a specific mismatch or
    positive-control violation can be constructed deliberately."""
    import pyarrow as pa
    import pyarrow.parquet as pq2

    clusters = sorted({c for _, c in incidence_rows})
    pq2.write_table(
        pa.table({"asmid": [a for a, _ in incidence_rows],
                  "cluster_id": [c for _, c in incidence_rows]}),
        tmp_path / "incidence_matrix.parquet")
    pq2.write_table(
        pa.table({"cluster_id": clusters,
                  "external_status": ["no_hit"] * len(clusters)}),
        tmp_path / "cluster_external_status.parquet")
    pq2.write_table(
        pa.table({key: [r[key] for r in genome_rows] for key in genome_rows[0]}),
        tmp_path / "genome_metadata.parquet")


def _run_cli(tmp_path, extra_args=()):
    return subprocess.run(
        [sys.executable, "bin/accumulation_curve.py",
         "--incidence-matrix", str(tmp_path / "incidence_matrix.parquet"),
         "--external-status", str(tmp_path / "cluster_external_status.parquet"),
         "--genome-metadata", str(tmp_path / "genome_metadata.parquet"),
         "--n-permutations", "10", "--seed", "1",
         "--version-tag", "test-v1-abc1234",
         "--outdir", str(tmp_path), *extra_args],
        capture_output=True, text=True,
    )


def _genome_row(asmid, component, cluster_class, n_proteins):
    return {"asmid": asmid, "n_proteins": n_proteins, "complete_pct": 95.0,
            "n50_bp": 100000, "component_id": component,
            "cluster_class": cluster_class, "clade_rank": "PHYLUM",
            "clade_label": "Ascomycota"}


def test_positive_control_violation_exits_nonzero(tmp_path):
    """Validation #4 is wired into accumulation_curve.py as an automated
    assertion (design doc: "Implement as an automated assertion in
    bin/accumulation_curve.py's output validation"). Two same-component
    redundant_strain genomes with completely disjoint protein content are
    exactly the impossible case that assertion exists to catch."""
    incidence_rows = ([("GENOME_A", f"cA{i}") for i in range(10)]
                      + [("GENOME_B", f"cB{i}") for i in range(10)])
    genome_rows = [_genome_row("GENOME_A", "comp_1", "redundant_strain", 10),
                   _genome_row("GENOME_B", "comp_1", "redundant_strain", 10)]
    _write_synthetic_inputs(tmp_path, incidence_rows, genome_rows)

    result = _run_cli(tmp_path)
    assert result.returncode == 1, result.stdout
    assert "VIOLATION" in result.stderr
    assert "GENOME_" in result.stderr


def test_positive_control_passes_for_genuinely_redundant_pair(tmp_path):
    """Same setup, but the two same-component genomes actually share all of
    their protein content -- the later one contributes nothing new, so the
    assertion passes and the CLI exits 0. A third, unrelated genome with its
    own content is included so the per-permutation power-law fit has more
    than one distinct marginal value to work with."""
    incidence_rows = ([(g, f"c{i}") for g in ("GENOME_A", "GENOME_B") for i in range(10)]
                      + [("GENOME_C", f"cC{i}") for i in range(10)])
    genome_rows = [_genome_row("GENOME_A", "comp_1", "redundant_strain", 10),
                   _genome_row("GENOME_B", "comp_1", "redundant_strain", 10),
                   _genome_row("GENOME_C", "comp_2", "singleton_isolated", 10)]
    _write_synthetic_inputs(tmp_path, incidence_rows, genome_rows)

    result = _run_cli(tmp_path)
    assert result.returncode == 0, result.stderr


def test_mismatched_genome_universes_raise(tmp_path):
    """incidence_matrix.parquet and genome_metadata.parquet are independent
    Nextflow path inputs; a genome in one but not the other would otherwise
    be silently dropped from clade_contribution.tsv while still counted in
    the curves."""
    incidence_rows = [(g, f"c{i}") for g in ("GENOME_A", "GENOME_B") for i in range(3)]
    genome_rows = [_genome_row("GENOME_A", "comp_1", "singleton_isolated", 3)]
    _write_synthetic_inputs(tmp_path, incidence_rows, genome_rows)

    result = _run_cli(tmp_path)
    assert result.returncode != 0
    assert "genome universes disagree" in result.stderr
    assert "GENOME_B" in result.stderr
