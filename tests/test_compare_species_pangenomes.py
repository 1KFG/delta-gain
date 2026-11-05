import csv, subprocess, sys
from pathlib import Path
import pyarrow as pa
import pyarrow.parquet as pq

BIN = Path(__file__).resolve().parent.parent / "bin"

# Genus G1 has Sp1 and Sp2; genus G2 has Sp3. Each species: 2 genomes.
# Sp1 clusters {a,b,s1}, Sp2 {a,b,s2}, Sp3 {a,z}. Cluster "a" is in all three.
# Sp1/Sp2 share {a,b}: shared=2, union=4 -> jaccard 0.5. Sp1/Sp3 share {a}: 1/(3+2-1)=0.25.
def test_sharing_and_taxon_summary(tmp_path):
    spec = {"Sp1": ("G1", ["a", "b", "s1"]), "Sp2": ("G1", ["a", "b", "s2"]),
            "Sp3": ("G2", ["a", "z"])}
    pairs, samples = [], [["ASMID", "SPECIES", "GENUS", "FAMILY", "ORDER", "CLASS", "PHYLUM"]]
    for sp, (genus, clusters) in spec.items():
        for g in (1, 2):
            asm = f"{sp}_{g}"
            samples.append([asm, sp, genus, "F", "O", "C", "P"])
            pairs += [(asm, c) for c in clusters]
    pq.write_table(pa.table({"asmid": [p[0] for p in pairs], "cluster_id": [p[1] for p in pairs]}),
                   tmp_path / "inc.parquet")
    csv.writer(open(tmp_path / "samples.csv", "w", newline="")).writerows(samples)
    sp_tsv = tmp_path / "sp.tsv"
    subprocess.run([sys.executable, str(BIN / "species_pangenome.py"),
                    "--incidence-matrix", str(tmp_path / "inc.parquet"),
                    "--samples-csv", str(tmp_path / "samples.csv"),
                    "--min-genomes", "2", "--rarefy-n", "2", "--n-permutations", "20",
                    "--out", str(sp_tsv)], check=True)
    subprocess.run([sys.executable, str(BIN / "compare_species_pangenomes.py"),
                    "--incidence-matrix", str(tmp_path / "inc.parquet"),
                    "--samples-csv", str(tmp_path / "samples.csv"),
                    "--species-table", str(sp_tsv), "--outdir", str(tmp_path)], check=True)
    pr = {(r["species_a"], r["species_b"]): r for r in
          csv.DictReader(open(tmp_path / "species_cluster_sharing.tsv"), delimiter="\t")}
    assert pr[("Sp1", "Sp2")]["shared_clusters"] == "2"
    assert abs(float(pr[("Sp1", "Sp2")]["jaccard"]) - 0.5) < 1e-9
    assert pr[("Sp1", "Sp2")]["same_genus"] == "1"
    assert abs(float(pr[("Sp1", "Sp3")]["jaccard"]) - 0.25) < 1e-9
    assert pr[("Sp1", "Sp3")]["same_genus"] == "0"
    tx = {(r["rank"], r["taxon"]): r for r in
          csv.DictReader(open(tmp_path / "taxon_pangenome_summary.tsv"), delimiter="\t")}
    assert tx[("GENUS", "G1")]["n_species"] == "2"
    assert tx[("PHYLUM", "P")]["n_species"] == "3"
