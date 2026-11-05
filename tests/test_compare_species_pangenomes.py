import csv
import subprocess
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

BIN = Path(__file__).resolve().parent.parent / "bin"

# Aa one {a,b,s1}, Aa two {a,b,s2} (genus Aa); Bb one {a,z} (genus Bb, 3 genomes, one excluded).
SPEC = {"Aa one": ("Aa", ["a", "b", "s1"], 2), "Aa two": ("Aa", ["a", "b", "s2"], 2),
        "Bb one": ("Bb", ["a", "z"], 3)}


def _setup(tmp):
    pairs, samples = [], [["ASMID", "SPECIES", "GENUS", "FAMILY", "ORDER", "CLASS", "PHYLUM"]]
    for sp, (genus, clusters, n) in SPEC.items():
        for g in range(1, n + 1):
            asm = f"{sp}_{g}"
            samples.append([asm, sp, genus, "F", "O", "C", "P"])
            pairs += [(asm, c) for c in clusters]
    pq.write_table(pa.table({"asmid": [p[0] for p in pairs], "cluster_id": [p[1] for p in pairs]}),
                   tmp / "inc.parquet")
    with open(tmp / "samples.csv", "w", newline="") as fh:
        csv.writer(fh).writerows(samples)
    (tmp / "excl.txt").write_text("Bb one_3\n")
    subprocess.run([sys.executable, str(BIN / "species_pangenome.py"),
                    "--incidence-matrix", str(tmp / "inc.parquet"), "--samples-csv", str(tmp / "samples.csv"),
                    "--exclude-asmids", str(tmp / "excl.txt"), "--min-genomes", "2", "--rarefy-n", "2",
                    "--n-permutations", "20", "--out", str(tmp / "sp.tsv")], check=True, capture_output=True)


def _compare(tmp, *extra):
    return subprocess.run([sys.executable, str(BIN / "compare_species_pangenomes.py"),
                           "--incidence-matrix", str(tmp / "inc.parquet"), "--samples-csv", str(tmp / "samples.csv"),
                           "--species-table", str(tmp / "sp.tsv"), "--outdir", str(tmp), *extra],
                          capture_output=True, text=True)


def test_exact_and_narrow_sharing_and_taxon_summary(tmp_path):
    _setup(tmp_path)
    r = _compare(tmp_path, "--exclude-asmids", str(tmp_path / "excl.txt"), "--max-species-per-cluster", "2")
    assert r.returncode == 0, r.stderr
    pr = {(x["species_a"], x["species_b"]): x for x in
          csv.DictReader(open(tmp_path / "species_cluster_sharing.tsv"), delimiter="\t")}
    p = pr[("Aa one", "Aa two")]       # exact: {a,b} shared, union 4; narrow drops cluster a (3 species > cap 2)
    assert p["shared_clusters"] == "2" and abs(float(p["jaccard"]) - 0.5) < 1e-9
    assert p["shared_clusters_narrow"] == "1" and abs(float(p["jaccard_narrow"]) - 1 / 3) < 1e-6
    assert p["same_genus"] == "1"
    q = pr[("Aa one", "Bb one")]       # only the broad cluster a is shared
    assert q["shared_clusters"] == "1" and abs(float(q["jaccard"]) - 0.25) < 1e-9
    assert q["shared_clusters_narrow"] == "0" and q["jaccard_narrow"] == "0" and q["same_genus"] == "0"
    tx = {(x["rank"], x["taxon"]): x for x in
          csv.DictReader(open(tmp_path / "taxon_pangenome_summary.tsv"), delimiter="\t")}
    g = tx[("GENUS", "Aa")]
    assert (g["n_species"], g["n_closed"], g["n_open"]) == ("2", "2", "0")
    assert float(g["frac_open"]) == 0.0            # closed species stay in the denominator
    assert tx[("PHYLUM", "P")]["n_species"] == "3"


def test_fails_when_exclusions_differ_between_steps(tmp_path):
    _setup(tmp_path)                                # species table built WITH the exclusion
    r = _compare(tmp_path)                          # compare run WITHOUT it
    assert r.returncode != 0 and "differ from the species table" in r.stderr
