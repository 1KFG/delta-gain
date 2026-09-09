from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from build_cluster_incidence_matrix import build_incidence_matrix

FIXTURES = Path("tests/data/accumulation")


def _write_qc_fixtures(tmp_path):
    busco = pa.table({"ASMID": ["GENOME_A", "GENOME_B", "GENOME_C", "GENOME_D"],
                       "complete_pct": [98.5, 97.0, 60.0, 95.0]})
    asm = pa.table({"ASMID": ["GENOME_A", "GENOME_B", "GENOME_C", "GENOME_D"],
                     "N50_bp": [500000, 480000, 20000, 300000]})
    pq.write_table(busco, tmp_path / "busco_genome.parquet")
    pq.write_table(asm, tmp_path / "asm_stats.parquet")
    return tmp_path / "busco_genome.parquet", tmp_path / "asm_stats.parquet"


def test_incidence_matrix_shape_and_membership(tmp_path):
    busco_pq, asm_pq = _write_qc_fixtures(tmp_path)
    incidence, ext_status, genome_meta = build_incidence_matrix(
        clusters_tsv=FIXTURES / "bfd_proteins_clusters.tsv",
        provenance_tsv=FIXTURES / "bfd_proteins_provenance.tsv",
        novelty_bins_tsv=FIXTURES / "bfd_novelty_bins.tsv",
        samples_csv=FIXTURES / "samples.csv",
        busco_parquet=busco_pq,
        asm_stats_parquet=asm_pq,
        genome_classification_tsv=FIXTURES / "genome_classification.tsv",
        min_clade_n=2,
    )
    # GENOME_A__p1 is the representative for both GENOME_A's and GENOME_B's
    # copy of that protein -> both genomes incident to cluster GENOME_A__p1.
    pairs = set(zip(incidence["asmid"], incidence["cluster_id"]))
    assert ("GENOME_A", "GENOME_A__p1") in pairs
    assert ("GENOME_B", "GENOME_A__p1") in pairs
    assert ("GENOME_A", "GENOME_A__p2") in pairs
    assert ("GENOME_B", "GENOME_B__p2") in pairs
    assert ("GENOME_C", "GENOME_C__p1") in pairs
    assert ("GENOME_D", "GENOME_D__p1") in pairs
    assert len(pairs) == 6


def test_external_status_lookup_keyed_by_cluster_representative(tmp_path):
    busco_pq, asm_pq = _write_qc_fixtures(tmp_path)
    _, ext_status, _ = build_incidence_matrix(
        clusters_tsv=FIXTURES / "bfd_proteins_clusters.tsv",
        provenance_tsv=FIXTURES / "bfd_proteins_provenance.tsv",
        novelty_bins_tsv=FIXTURES / "bfd_novelty_bins.tsv",
        samples_csv=FIXTURES / "samples.csv",
        busco_parquet=busco_pq, asm_stats_parquet=asm_pq,
        genome_classification_tsv=FIXTURES / "genome_classification.tsv",
        min_clade_n=2,
    )
    status = dict(zip(ext_status["cluster_id"], ext_status["external_status"]))
    assert status["GENOME_A__p2"] == "no_hit"
    assert status["GENOME_C__p1"] == "no_hit"


def test_genome_metadata_joins_qc_and_taxonomy_and_rolls_up_clade(tmp_path):
    busco_pq, asm_pq = _write_qc_fixtures(tmp_path)
    _, _, genome_meta = build_incidence_matrix(
        clusters_tsv=FIXTURES / "bfd_proteins_clusters.tsv",
        provenance_tsv=FIXTURES / "bfd_proteins_provenance.tsv",
        novelty_bins_tsv=FIXTURES / "bfd_novelty_bins.tsv",
        samples_csv=FIXTURES / "samples.csv",
        busco_parquet=busco_pq, asm_stats_parquet=asm_pq,
        genome_classification_tsv=FIXTURES / "genome_classification.tsv",
        min_clade_n=2,
    )
    meta = genome_meta.set_index("asmid")
    assert meta.loc["GENOME_A", "complete_pct"] == 98.5
    assert meta.loc["GENOME_A", "n50_bp"] == 500000
    assert meta.loc["GENOME_A", "component_id"] == "comp_1"
    assert meta.loc["GENOME_A", "cluster_class"] == "redundant_strain"
    # ORDER=Pleosporales has 2 genomes (>= min_clade_n=2) -> that's the clade
    # rank/label. FAMILY is blank in the fixture, so ORDER is this genome's
    # finest POPULATED rank: no roll-up happened, label stays bare.
    assert meta.loc["GENOME_A", "clade_rank"] == "ORDER"
    assert meta.loc["GENOME_A", "clade_label"] == "Pleosporales"
    # ORDER=Agaricales has only 1 genome (< min_clade_n), but CLASS=
    # Agaricomycetes has 2 (GENOME_C + GENOME_D) -> rolls up to CLASS. ORDER
    # WAS populated and failed the density check, so this is a genuine
    # roll-up and the label must say so (design doc: Components).
    assert meta.loc["GENOME_C", "clade_rank"] == "CLASS"
    assert meta.loc["GENOME_C", "clade_label"] == "Agaricomycetes (other orders)"
    # GENOME_D: ORDER=Auriculariales has only 1 genome -> same CLASS rollup.
    assert meta.loc["GENOME_D", "clade_rank"] == "CLASS"
    assert meta.loc["GENOME_D", "clade_label"] == "Agaricomycetes (other orders)"


def test_duplicate_asmid_in_samples_csv_raises(tmp_path):
    import pytest
    bad_samples = tmp_path / "samples.csv"
    bad_samples.write_text(
        "ASMID,SPECIES_IN,STRAIN,BIOPROJECT,NCBI_TAXONID,BUSCO_LINEAGE,"
        "PHYLUM,SUBPHYLUM,CLASS,SUBCLASS,ORDER,FAMILY,GENUS,SPECIES,"
        "TRANSL_TABLE,LOCUSTAG\n"
        "GENOME_A,Foo bar,S1,PRJ1,1,dikarya,Ascomycota,,Dothideomycetes,,"
        "Pleosporales,,Foo,Foo bar,1,X1\n"
        "GENOME_A,Foo bar,S1,PRJ1,1,dikarya,Ascomycota,,Dothideomycetes,,"
        "Pleosporales,,Foo,Foo bar,1,X1\n"  # exact duplicate ASMID row
    )
    busco_pq, asm_pq = _write_qc_fixtures(tmp_path)
    with pytest.raises(ValueError, match="duplicate"):
        build_incidence_matrix(
            clusters_tsv=FIXTURES / "bfd_proteins_clusters.tsv",
            provenance_tsv=FIXTURES / "bfd_proteins_provenance.tsv",
            novelty_bins_tsv=FIXTURES / "bfd_novelty_bins.tsv",
            samples_csv=bad_samples,
            busco_parquet=busco_pq, asm_stats_parquet=asm_pq,
            genome_classification_tsv=FIXTURES / "genome_classification.tsv",
            min_clade_n=2,
        )


def test_duplicate_asmid_in_provenance_raises(tmp_path):
    import pytest
    bad_provenance = tmp_path / "bfd_proteins_provenance.tsv"
    bad_provenance.write_text(
        "new_id\tasmid\toriginal_id\n"
        "GENOME_A__p1\tGENOME_A\tp1\n"
        "GENOME_A__p1\tGENOME_A\tp1\n"  # exact duplicate row -> still a real dup
    )
    busco_pq, asm_pq = _write_qc_fixtures(tmp_path)
    with pytest.raises(ValueError, match="duplicate"):
        build_incidence_matrix(
            clusters_tsv=FIXTURES / "bfd_proteins_clusters.tsv",
            provenance_tsv=bad_provenance,
            novelty_bins_tsv=FIXTURES / "bfd_novelty_bins.tsv",
            samples_csv=FIXTURES / "samples.csv",
            busco_parquet=busco_pq, asm_stats_parquet=asm_pq,
            genome_classification_tsv=FIXTURES / "genome_classification.tsv",
            min_clade_n=2,
        )


def test_cluster_missing_from_novelty_bins_raises(tmp_path):
    """A cluster representative absent from bfd_novelty_bins.tsv means the
    three inputs came from mismatched runs -- bin_novelty_hits.py emits a row
    for every query including no-hit ones. Must fail loudly rather than
    silently defaulting to the novelty-inflating "no_hit"."""
    import pytest
    truncated_bins = tmp_path / "bfd_novelty_bins.tsv"
    truncated_bins.write_text(
        "query_id\toverall_status\n"
        "GENOME_A__p1\tstrong_hit\n"
        "GENOME_A__p2\tno_hit\n"
        "GENOME_B__p2\tweak_hit\n"
        # GENOME_C__p1 and GENOME_D__p1 deliberately omitted.
    )
    busco_pq, asm_pq = _write_qc_fixtures(tmp_path)
    with pytest.raises(ValueError, match="GENOME_C__p1"):
        build_incidence_matrix(
            clusters_tsv=FIXTURES / "bfd_proteins_clusters.tsv",
            provenance_tsv=FIXTURES / "bfd_proteins_provenance.tsv",
            novelty_bins_tsv=truncated_bins,
            samples_csv=FIXTURES / "samples.csv",
            busco_parquet=busco_pq, asm_stats_parquet=asm_pq,
            genome_classification_tsv=FIXTURES / "genome_classification.tsv",
            min_clade_n=2,
        )
