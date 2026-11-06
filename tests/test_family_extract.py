"""Hand-checked regions for bin/family_extract.py on a tiny synthetic BFD."""
import csv
import gzip
import os
import sys

import duckdb

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
import family_extract  # noqa: E402

SEQ = "ACDEFGHIKLMNPQRSTVWY" * 3  # 60 aa, position-identifiable


def _pq(con, path, cols, rows):
    con.execute(f"create or replace table t ({', '.join(cols)})")
    for r in rows:
        con.execute("insert into t values (" + ",".join("?" * len(r)) + ")", r)
    con.execute(f"copy t to '{path}' (format parquet)")


def build(tmp):
    t = tmp / "tables"
    t.mkdir()
    con = duckdb.connect()
    _pq(con, f"{t}/species.parquet", ["ASMID varchar", "LOCUSTAG varchar", "SPECIES varchar"],
        [("A1", "PA1", "sp1"), ("N1", "PN1", "sp2"), ("X", "PX", "sp3")])
    with open(tmp / "classes.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["ASMID", "accession", "class", "ncbi_status", "pep_status"])
        for a, c in (("A1", "baseline"), ("N1", "new"), ("X", "annotated_pending")):
            w.writerow([a, a, c, "s", "matched"])
    _pq(con, f"{t}/pfam.parquet",
        ["protein_id varchar", "pfam_id varchar", "pfam_acc varchar", "env_from bigint", "env_to bigint",
         "domain_score double", "domain_i_evalue double", "hmm_from bigint", "hmm_to bigint", "pfam_len bigint"],
        [("PA1_1-T1", "X1", "PF00001.5", 5, 50, 80.0, 1e-20, 1, 40, 50),
         ("PA1_1-T1", "X1", "PF00001.5", 52, 60, 20.0, 1e-3, 1, 8, 50),   # < MIN_REGION -> dropped
         ("PN1_1-T1", "X1", "PF00001.5", 1, 40, 70.0, 1e-18, 1, 40, 50),
         ("PX_1-T1", "X1", "PF00001.5", 1, 40, 70.0, 1e-18, 1, 40, 50)])  # excluded genome
    _pq(con, f"{t}/cazy.cazymes_hmm.parquet",
        ["protein_id varchar", "HMM_id varchar", "s_start bigint", "s_end bigint", "evalue double",
         "coverage double"],
        [("PA1_1-T1", "GH5_7", 2, 40, 1e-30, 0.9), ("PA1_1-T1", "GH5", 30, 58, 1e-20, 0.8),
         ("PN1_1-T1", "GH5", 1, 45, 1e-25, 0.95)])
    _pq(con, f"{t}/cazy.overview.parquet",
        ["protein_id varchar", "EC varchar", "sub_fam varchar", "substrate varchar"],
        [("PA1_1-T1", "3.2.1.4", "GH5_e1", "cellulose")])
    _pq(con, f"{t}/merops.parquet",
        ["protein_id varchar", "merops_id varchar", "q_start bigint", "q_end bigint", "evalue double",
         "percent_identity double", "bitscore double"],
        [("PA1_1-T1", "MER2", 1, 35, 1e-40, 70.0, 100.0), ("PA1_1-T1", "MER2", 20, 58, 1e-10, 40.0, 50.0),
         ("PN1_1-T1", "MER1", 1, 40, 1e-9, 33.0, 40.0)])
    (tmp / "mmap.tab").write_text("MER1\tS01A\tS01.001\nMER2\tS08A\tS08.001\n")
    man = tmp / "manifest.tsv"
    with open(man, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["ASMID", "tag", "pep_path", "status"])
        for asm, pre in (("A1", "PA1"), ("N1", "PN1"), ("X", "PX")):
            p = tmp / f"{pre}.faa.gz"
            with gzip.open(p, "wt") as g:
                g.write(f">{pre}_1-T1 some description\n{SEQ[:30]}\n{SEQ[30:]}*\n>{pre}_2-T1\n{SEQ}\n")
            w.writerow([asm, asm, str(p), "matched"])
    return t, man


def members(path):
    with open(path) as fh:
        return {r["region_id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def test_extract(tmp_path):
    t, man = build(tmp_path)
    out = tmp_path / "out"
    family_extract.main(["--tables-dir", str(t), "--classes", str(tmp_path / "classes.tsv"),
                         "--manifest", str(man), "--merops-map", str(tmp_path / "mmap.tab"),
                         "--outdir", str(out), "--threads", "2",
                         "--families", "pfam:PF00001", "cazy_family:GH5", "merops_family:S08"])
    m = members(out / "pfam_PF00001.members.tsv")
    assert sorted(m) == ["PA1_1-T1|5-50", "PN1_1-T1|1-40"]  # short domain + excluded genome gone
    r = m["PA1_1-T1|5-50"]
    assert (r["asmid"], r["class"], r["region_len"], r["protein_len"], r["model_cov"]) == \
        ("A1", "baseline", "46", "60", "0.8")
    fa = (out / "pfam_PF00001.faa").read_text().split("\n")
    assert fa[fa.index(">PA1_1-T1|5-50") + 1] == SEQ[4:50]  # 1-based inclusive slice, '*' stripped

    g = members(out / "cazy_family_GH5.members.tsv")
    assert "PA1_1-T1|2-58" in g and "PA1_1-T1|2-40" not in g  # overlapping hits merged
    a = g["PA1_1-T1|2-58"]
    assert (a["sub"], a["ec"], a["efam"], a["substrate"]) == ("GH5_7", "3.2.1.4", "GH5_e1", "cellulose")
    assert g["PN1_1-T1|1-45"]["ec"] == "NA"

    s = members(out / "merops_family_S08.members.tsv")
    assert list(s) == ["PA1_1-T1|1-58"] and s["PA1_1-T1|1-58"]["pident"] == "70.0"
    assert not (out / "merops_family_S01.faa").exists()  # PN1's S01 hit is not an S08 request

    summ = {r["family"]: r for r in csv.DictReader(open(out / "regions_summary.tsv"), delimiter="\t")}
    assert summ["pfam_PF00001"]["genomes_baseline"] == "1" and summ["pfam_PF00001"]["genomes_new"] == "1"
