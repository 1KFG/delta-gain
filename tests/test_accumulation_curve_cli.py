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

    with open(tmp_path / "pangenome_powerlaw_fit.tsv") as fh:
        fit_rows = list(csv.DictReader(fh, delimiter="\t"))
    assert {r["series"] for r in fit_rows} == {"internal", "external"}

    with open(tmp_path / "clade_contribution.tsv") as fh:
        clade_rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(clade_rows) >= 1

    raw = pq.read_table(tmp_path / "permutation_marginals.parquet")
    assert raw.num_rows == 20 * 4  # n_permutations * n_genomes
