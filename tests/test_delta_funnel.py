import csv, subprocess, sys, gzip
from pathlib import Path

BIN = Path(__file__).resolve().parent.parent / "bin"
PY = sys.executable


def run(script, *args):
    r = subprocess.run([PY, str(BIN / script), *map(str, args)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r.stderr


def write_status(p, rows):
    with open(p, "w") as f:
        f.write("accession\tfolder\tannotation_provider\tannotation_release\tprotein_coding\tprotein_file_ok\tstatus\tassembly_release\tsource_database\n")
        for acc, folder, ok, st in rows:
            f.write(f"{acc}\t{folder}\t\t\t\t{ok}\t{st}\t\t\n")


def test_classify_and_build(tmp_path):
    st = tmp_path / "status.tsv"
    write_status(st, [
        ("GCA_1.1", "GCA_1.1_b", 1, "baseline_ok"),
        ("GCA_2.1", "GCA_2.1_n", 0, "no_prior_annotation_present"),
        ("GCA_3.1", "GCA_3.1_n", 0, "no_prior_annotation_present"),
        ("GCA_4.1", "GCA_4.1_p", 0, "annotated_folder_missing_protein"),
    ])
    sib = tmp_path / "sib.txt"; sib.write_text("GCA_3.1\n")
    pep = {}
    for n in ["b", "n", "n3", "p"]:
        f = tmp_path / f"{n}.fa"; f.write_text(f">{n}1 d\nMKVLAAGIVGL\n>{n}2\nMKV\n>{n}1\nMKVLAAGIVGL\n"); pep[n] = f
    man = tmp_path / "man.tsv"
    man.write_text("ASMID\ttag\tpep_path\tstatus\n"
                   f"GCA_1.1_b\tb\t{pep['b']}\tmatched\n"
                   f"GCA_2.1_n\tn\t{pep['n']}\tmatched\n"
                   f"GCA_3.1_n\tn3\t{pep['n3']}\tmatched\n"
                   f"GCA_4.1_p\tp\t{pep['p']}\tmatched\n"
                   f"PRIVATE_x\tx\t{pep['p']}\tmatched\n")
    err = run("classify_genomes_by_prior_annotation.py", "--pep-manifest", man, "--status", st,
              "--sibling-flags", sib, "--out-classes", tmp_path / "cls.tsv",
              "--out-new-manifest", tmp_path / "new.tsv")
    cls = {r["ASMID"]: r["class"] for r in csv.DictReader(open(tmp_path / "cls.tsv"), delimiter="\t")}
    assert cls == {"GCA_1.1_b": "baseline", "GCA_2.1_n": "new", "GCA_3.1_n": "new_sibling_flagged",
                   "GCA_4.1_p": "annotated_pending", "PRIVATE_x": "not_in_ncbi_table"}
    matched = [r["ASMID"] for r in csv.DictReader(open(tmp_path / "new.tsv"), delimiter="\t") if r["status"] == "matched"]
    assert matched == ["GCA_2.1_n"]
    # baseline protein file in the NCBI layout
    asm = tmp_path / "asm" / "GCA_1.1_b"; asm.mkdir(parents=True)
    with gzip.open(asm / "GCA_1.1_b_protein.faa.gz", "wt") as g:
        g.write(">b1 x\nMKVLAAGIVGL*\n>short\nMK\n")
    run("build_funnel_input.py", "--status", st, "--asm-dir", tmp_path / "asm", "--new-manifest", tmp_path / "new.tsv",
        "--out-fasta", tmp_path / "all.faa", "--out-prov", tmp_path / "prov.tsv", "--min-len", 5)
    ids = [l[1:].strip() for l in open(tmp_path / "all.faa") if l.startswith(">")]
    assert ids == ["BASE__GCA_1.1_b__b1", "NEW__GCA_2.1_n__n1"]   # short + duplicate dropped, '*' stripped


def test_summary_nested_tiers(tmp_path):
    prov = tmp_path / "prov.tsv"
    prov.write_text("id\tclass\tgenome\n"
                    "B1\tBASE\tgb\n"
                    "N1\tNEW\tg1\nN2\tNEW\tg1\nN3\tNEW\tg2\nN4\tNEW\tg2\nN5\tNEW\tg2\n")
    # tier0 clusters: {B1,N1}->known ; {N2}; {N3,N4}; {N5}
    t0 = tmp_path / "t0.tsv"
    t0.write_text("B1\tB1\nB1\tN1\nN2\tN2\nN3\tN3\nN3\tN4\nN5\tN5\n")
    # tier1 over reps {B1,N2,N3,N5}: N2 joins B1 (becomes known), N3+N5 together, no base
    t1 = tmp_path / "t1.tsv"
    t1.write_text("B1\tB1\nB1\tN2\nN3\tN3\nN3\tN5\n")
    out = tmp_path / "out"
    run("delta_funnel_summary.py", "--prov", prov, "--tier-clusters", t0, t1,
        "--tier-labels", "95", "70", "--outdir", out, "--threads", 2, "--memory-limit", "1GB")
    ts = {r["tier"]: r for r in csv.DictReader(open(out / "tier_summary.tsv"), delimiter="\t")}
    assert (ts["95"]["new_proteins"], ts["95"]["novel_proteins"]) == ("5", "4")   # only N1 known
    assert (ts["70"]["new_proteins"], ts["70"]["novel_proteins"]) == ("5", "3")   # N1,N2 known
    g = {r["genome"]: r for r in csv.DictReader(open(out / "genome_tier.tsv"), delimiter="\t")}
    assert g["g1"]["novel_95"] == "1" and g["g1"]["novel_70"] == "0"
    assert g["g2"]["novel_95"] == "3" and g["g2"]["novel_70"] == "3"
    sv = list(csv.DictReader(open(out / "survivors_clusters.tsv"), delimiter="\t"))
    assert len(sv) == 1 and sv[0]["rep"] == "N3" and sv[0]["n_new_proteins"] == "3" and sv[0]["n_genomes"] == "1"
