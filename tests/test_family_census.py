"""Hand-checked numbers for bin/family_census.py on a tiny synthetic BFD."""
import csv
import os
import sys

import duckdb

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
import family_census  # noqa: E402


def _pq(con, path, cols, rows):
    con.execute(f"create or replace table t ({', '.join(cols)})")
    for r in rows:
        con.execute("insert into t values (" + ",".join("?" * len(r)) + ")", r)
    con.execute(f"copy t to '{path}' (format parquet)")


def build(tmp):
    t = tmp / "tables"
    t.mkdir()
    con = duckdb.connect()
    # genomes: A1,A2 baseline sp1 ; N1 new sp1 ; N2,N3 new sp2 (new taxon) ; X excluded
    _pq(con, f"{t}/species.parquet", ["ASMID varchar", "LOCUSTAG varchar", "SPECIES varchar"],
        [("A1", "PA1", "sp1"), ("A2", "PA2", "sp1"), ("N1", "PN1", "sp1"),
         ("N2", "PN2", "sp2"), ("N3", "PN3", "sp2"), ("X", "PX", "sp3")])
    with open(tmp / "classes.tsv", "w", newline="") as fh:  # CRLF like the real file
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["ASMID", "accession", "class", "ncbi_status", "pep_status"])
        for a, c, p in (("A1", "baseline", "matched"), ("A2", "baseline", "matched"),
                        ("N1", "new", "matched"), ("N2", "new", "matched"),
                        ("N3", "new", "matched"), ("X", "annotated_pending", "matched")):
            w.writerow([a, a, c, "s", p])
    _pq(con, f"{t}/pfam.parquet", ["protein_id varchar", "pfam_id varchar", "pfam_acc varchar"],
        [("PA1_1-T1", "X1", "PF00001.5"), ("PA1_1-T1", "X1", "PF00001.5"), ("PA1_2-T1", "X1", "PF00001.5"),
         ("PA2_1-T1", "X1", "PF00001.5"), ("PN1_1-T1", "X1", "PF00001.5"),
         ("PN2_1-T1", "X2", "PF00002.3"), ("PN3_1-T1", "X2", "PF00002.3"),
         ("PX_1-T1", "X3", "PF00003.1")])
    _pq(con, f"{t}/cazy.cazymes_hmm.parquet", ["protein_id varchar", "HMM_id varchar"],
        [("PA1_1-T1", "GH5_7"), ("PN1_1-T1", "GH5_7"), ("PN2_1-T1", "GH5")])
    _pq(con, f"{t}/cazy.overview.parquet", ["protein_id varchar", "sub_fam varchar"],
        [("PN1_1-T1", "CBM35_e7+GH27_e9"), ("PN2_1-T1", "-")])
    _pq(con, f"{t}/merops.parquet",
        ["protein_id varchar", "merops_id varchar", "percent_identity double", "evalue double"],
        [("PA1_1-T1", "MER1", 30.0, 1e-20), ("PA1_1-T1", "MER1", 50.0, 1e-30),
         ("PN2_1-T1", "MER1", 35.0, 1e-12), ("PN2_1-T1", "MER2", 80.0, 1e-50)])
    (tmp / "mmap.tab").write_text("MER1\tS01A\tS01.001\nMER2\tS08A\tS08.001\n")
    return t


def rows(path):
    with open(path) as fh:
        return {(r["level"], r["family"]): r for r in csv.DictReader(fh, delimiter="\t")}


def test_census(tmp_path):
    t = build(tmp_path)
    out = tmp_path / "out"
    family_census.main(["--tables-dir", str(t), "--classes", str(tmp_path / "classes.tsv"),
                        "--merops-map", str(tmp_path / "mmap.tab"), "--outdir", str(out)])
    tot = dict(l.split("\t") for l in (out / "census_totals.tsv").read_text().split("\n")[1:] if l)
    assert tot == {"genomes_baseline": "2", "genomes_new": "3", "species_with_baseline": "1",
                   "species_with_new": "2", "species_new_taxon": "1"}

    p = rows(out / "pfam_census.tsv")
    assert ("pfam", "PF00003") not in p  # excluded genome never counted
    f1 = p[("pfam", "PF00001")]
    assert (f1["g_base"], f1["g_new"], f1["prot_base"], f1["inst_base"]) == ("2", "1", "3", "4")
    assert f1["copies_base"] == "1.5" and f1["new_only"] == "0"
    f2 = p[("pfam", "PF00002")]
    assert (f2["g_base"], f2["g_new"], f2["sp_new"], f2["sp_newtaxon"], f2["sp_gain_known"]) == \
        ("0", "2", "1", "1", "0")
    assert f2["new_only"] == "0"  # needs >= 3 new genomes

    c = rows(out / "cazy_census.tsv")
    assert (c[("cazy_family", "GH5")]["g_base"], c[("cazy_family", "GH5")]["g_new"]) == ("1", "2")
    assert c[("cazy_family", "GH5")]["sp_gain_known"] == "0"  # sp1 has it in base already
    assert (c[("cazy_subfamily", "GH5_7")]["g_base"], c[("cazy_subfamily", "GH5_7")]["g_new"]) == ("1", "1")
    assert c[("cazy_efam", "GH27_e9")]["g_new"] == "1" and c[("cazy_efam", "CBM35_e7")]["g_base"] == "0"
    assert ("cazy_efam", "-") not in c

    m = rows(out / "merops_census.tsv")
    s01 = m[("merops_family", "S01")]
    assert (s01["g_base"], s01["g_new"], s01["lowid_base"], s01["lowid_new"]) == ("1", "1", "0", "1")
    assert m[("merops_family", "S08")]["g_new"] == "1"
    assert m[("merops_subfamily", "S01A")]["g_base"] == "1"


def test_bh_monotone():
    q = family_census.bh([0.01, None, 0.04, 0.03])
    assert q[1] is None
    assert q[0] <= q[3] <= q[2] and all(0 <= x <= 1 for x in q if x is not None)
