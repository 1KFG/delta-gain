import csv
import subprocess
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

BIN = Path(__file__).resolve().parent.parent / "bin" / "species_pangenome.py"


def _write(tmp, genomes, exclude=()):
    """genomes: {asmid: (species, [clusters])}"""
    pairs = [(a, c) for a, (_, cl) in genomes.items() for c in cl]
    pq.write_table(pa.table({"asmid": [p[0] for p in pairs], "cluster_id": [p[1] for p in pairs]}),
                   tmp / "inc.parquet")
    with open(tmp / "samples.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["ASMID", "SPECIES"])
        for a, (sp, _) in genomes.items():
            w.writerow([a, sp])
    (tmp / "excl.txt").write_text("".join(f"{a}\n" for a in exclude))


def _run(tmp, *extra):
    out = tmp / "out.tsv"
    subprocess.run([sys.executable, str(BIN), "--incidence-matrix", str(tmp / "inc.parquet"),
                    "--samples-csv", str(tmp / "samples.csv"), "--n-permutations", "50",
                    "--out", str(out), *extra], check=True, capture_output=True)
    return {r["species"]: r for r in csv.DictReader(open(out), delimiter="\t")}


def test_small_n_classes_closed_species_and_placeholders(tmp_path):
    g = {}
    for i in range(5):   # c0 core; c1..c4 each in one genome
        g[f"A{i}"] = ("Genus alpha", ["c0"] + ([f"c{i}"] if i else []))
        g[f"B{i}"] = ("Genus beta", ["d0", "d1"])          # identical genomes -> closed
        g[f"F{i}"] = ("Fusarium sp.", ["f0", f"f{i}"])      # placeholder -> skipped
    g["X1"] = ("", ["c0"])                                   # no species -> skipped
    g["R1"] = ("Genus alpha", ["c0"])                        # excluded
    _write(tmp_path, g, exclude=["R1"])
    rows = _run(tmp_path, "--min-genomes", "5", "--rarefy-n", "5", "--exclude-asmids", str(tmp_path / "excl.txt"))
    assert set(rows) == {"Genus alpha", "Genus beta"}
    a, b = rows["Genus alpha"], rows["Genus beta"]
    assert (a["n_genomes"], a["pan_size"], a["core"], a["soft_core"], a["cloud"], a["shell"]) == \
        ("5", "5", "1", "", "4", "0")                        # singletons are cloud even at n=5
    assert a["singleton_clusters"] == "4"
    assert float(a["rarefied_pan_size"]) == 5.0 and float(a["rarefied_core"]) == 1.0
    assert (b["pan_size"], b["core"], b["fit_status"], b["openness"]) == ("2", "2", "no_new_after_first", "closed")


def test_soft_core_and_shell_at_n20(tmp_path):
    g = {}
    for i in range(20):
        cl = ["c0"]
        if i < 19: cl.append("x")        # 19/20 -> soft-core
        if i < 10: cl.append("z")        # 10/20 -> shell
        if i < 3: cl.append("y")         # 3/20 -> cloud (floor(0.15*20) = 3)
        g[f"G{i}"] = ("Genus gamma", cl)
    _write(tmp_path, g)
    r = _run(tmp_path, "--min-genomes", "20", "--rarefy-n", "10")["Genus gamma"]
    assert (r["pan_size"], r["core"], r["soft_core"], r["shell"], r["cloud"]) == ("4", "1", "1", "1", "1")


def test_alpha_excludes_first_genome_and_reports_points(tmp_path):
    g = {f"H{i}": ("Genus delta", ["shared"] + [f"p{i}_{k}" for k in range(3)]) for i in range(6)}
    _write(tmp_path, g)
    r = _run(tmp_path, "--min-genomes", "6", "--rarefy-n", "6")["Genus delta"]
    # every genome after the first adds exactly 3 clusters (the first adds 4): flat -> alpha 0.
    # Including N=1 would give a clearly positive alpha.
    assert r["fit_status"] == "ok" and r["alpha_n_points"] == "5"
    assert abs(float(r["alpha"])) < 1e-9 and r["openness"] == "open"
