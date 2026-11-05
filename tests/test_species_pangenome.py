import subprocess, sys, csv
from pathlib import Path
import pyarrow as pa
import pyarrow.parquet as pq

BIN = Path(__file__).resolve().parent.parent / "bin" / "species_pangenome.py"

# Species S1: 5 genomes. c0 in all (core); c1..c4 each in one genome (singletons).
# Species S2: 5 genomes sharing the same 2 clusters (closed pangenome, all core).
# Genome X1: unknown species (must be skipped, not crash). R1 is excluded.
def _write(tmp_path):
    pairs = []
    for i in range(5):
        pairs.append((f"A{i}", "c0"))
        if i > 0:
            pairs.append((f"A{i}", f"c{i}"))
    for i in range(5):
        pairs += [(f"B{i}", "d0"), (f"B{i}", "d1")]
    pairs += [("X1", "c0"), ("R1", "c0")]
    pq.write_table(pa.table({"asmid": [p[0] for p in pairs],
                              "cluster_id": [p[1] for p in pairs]}),
                   tmp_path / "inc.parquet")
    with open(tmp_path / "samples.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["ASMID", "SPECIES"])
        for i in range(5):
            w.writerow([f"A{i}", "Sp one"]); w.writerow([f"B{i}", "Sp two"])
        w.writerow(["R1", "Sp one"]); w.writerow(["X1", ""])
    (tmp_path / "excl.txt").write_text("R1\n")

def test_species_summary(tmp_path):
    _write(tmp_path)
    out = tmp_path / "out.tsv"
    subprocess.run([sys.executable, str(BIN),
                    "--incidence-matrix", str(tmp_path / "inc.parquet"),
                    "--samples-csv", str(tmp_path / "samples.csv"),
                    "--exclude-asmids", str(tmp_path / "excl.txt"),
                    "--min-genomes", "5", "--rarefy-n", "5",
                    "--n-permutations", "50", "--out", str(out)], check=True)
    rows = {r["species"]: r for r in csv.DictReader(open(out), delimiter="\t")}
    assert set(rows) == {"Sp one", "Sp two"}
    s1, s2 = rows["Sp one"], rows["Sp two"]
    assert (s1["n_genomes"], s1["pan_size"], s1["core"], s1["singleton_clusters"]) == ("5", "5", "1", "4")
    assert (s2["pan_size"], s2["core"]) == ("2", "2")
    assert float(s1["rarefied_pan_size"]) == 5.0   # all 5 genomes -> full pan
    assert float(s2["core_fraction"]) == 1.0
