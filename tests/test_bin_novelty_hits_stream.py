import csv, random, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / "bin" / "bin_novelty_hits.py"
NEW = HERE.parent / "bin" / "bin_novelty_hits_stream.py"
TARGETS = ["nr_cluster", "uniprot_fungi"]


def _row(q, ev, bits, qc, sc):
    return "\t".join([q, "subj", "50.0", "100", "0", "0", "1", "100", "1", "100",
                      ev, str(bits), str(qc), str(sc), "some title"])


def _make(tmp, seed=1, n_chunks=3, per_chunk=40):
    rnd = random.Random(seed)
    chunks, all_ids = [], []
    for c in range(n_chunks):
        ids = [f"GCF_{c}.1_Sp|prot_{c}_{i}" for i in range(per_chunk)]
        all_ids += ids
        fa = tmp / f"c{c}.faa"
        fa.write_text("".join(f">{q} desc here\nMKV\n" for q in ids))
        hits = {}
        for t in TARGETS:
            rows = []
            for q in ids:
                if rnd.random() < 0.2:          # no hit in this target
                    continue
                for _ in range(rnd.randint(1, 4)):
                    ev = rnd.choice(["0", "1e-50", "1e-8", "1e-5", "5e-4", "2e-2", "1e-8"])
                    qc = rnd.choice([30, 60, 70, 85, 99.5])
                    sc = rnd.choice([30, 60, 70, 90, 100])
                    rows.append(_row(q, ev, rnd.choice([100, 200, 300]), qc, sc))
            rnd.shuffle(rows)                    # interleave rows of different queries
            p = tmp / f"c{c}.{t}.tsv"
            p.write_text("\n".join(rows) + ("\n" if rows else ""))
            hits[t] = p
        chunks.append((f"c{c}", fa, hits))
    full = tmp / "all.faa"
    full.write_text("".join(f">{q}\nMKV\n" for q in all_ids))
    return chunks, full


def _read(path):
    with open(path) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    out = []
    for r in rows[1:]:
        out.append([float(x) if x.replace(".", "").replace("-", "").replace("e", "").replace("+", "").isdigit()
                    else x for x in r])
    return rows[0], out


def _run_both(tmp, mode="best-evalue", seed=1):
    chunks, full = _make(tmp, seed)
    old_out, new_out = tmp / "old.tsv", tmp / "new.tsv"
    cmd = [sys.executable, str(OLD), "--queries", str(full), "--out", str(old_out)]
    for _, _, hits in chunks:
        for t, p in hits.items():
            cmd += ["--hits", f"{t}={p}"]
    subprocess.run(cmd, check=True, capture_output=True)
    man = tmp / "manifest.tsv"
    man.write_text("chunk\tqueries\t" + "\t".join(TARGETS) + "\n" +
                   "".join(f"{c}\t{fa}\t" + "\t".join(str(h[t]) for t in TARGETS) + "\n"
                           for c, fa, h in chunks))
    subprocess.run([sys.executable, str(NEW), "--manifest", str(man), "--queries", str(full),
                    "--classify", mode, "--out", str(new_out)], check=True, capture_output=True)
    return _read(old_out), _read(new_out)


def test_default_matches_original_script(tmp_path):
    for seed in (1, 2, 3):
        d = tmp_path / str(seed); d.mkdir()
        (oh, orows), (nh, nrows) = _run_both(d, seed=seed)
        assert oh == nh
        assert orows == nrows, f"seed {seed}: streaming output differs from original"


def test_best_class_unmasks_strong_hit(tmp_path):
    fa = tmp_path / "c.faa"; fa.write_text(">q1\nMKV\n>q2\nMKV\n")
    # q1: partial-coverage hit with the best E-value, plus a full-coverage hit at 1e-8.
    nr = tmp_path / "nr.tsv"
    nr.write_text(_row("q1", "1e-50", 400, 40, 40) + "\n" + _row("q1", "1e-8", 150, 90, 90) + "\n")
    up = tmp_path / "up.tsv"; up.write_text("")
    man = tmp_path / "m.tsv"
    man.write_text(f"chunk\tqueries\tnr_cluster\tuniprot_fungi\nc0\t{fa}\t{nr}\t{up}\n")
    for mode, expect in (("best-evalue", "weak_hit"), ("best-class", "strong_hit")):
        out = tmp_path / f"{mode}.tsv"
        subprocess.run([sys.executable, str(NEW), "--manifest", str(man), "--classify", mode,
                        "--out", str(out)], check=True, capture_output=True)
        rows = {r[0]: r for r in csv.reader(open(out), delimiter="\t")}
        assert rows["q1"][-1] == expect, (mode, rows["q1"])
        assert rows["q2"][-1] == "no_hit"


def test_fails_loudly_on_foreign_query_and_missing_file(tmp_path):
    fa = tmp_path / "c.faa"; fa.write_text(">q1\nMKV\n")
    bad = tmp_path / "bad.tsv"; bad.write_text(_row("other_chunk_q", "1e-9", 100, 90, 90) + "\n")
    empty = tmp_path / "e.tsv"; empty.write_text("")
    man = tmp_path / "m.tsv"
    man.write_text(f"chunk\tqueries\tnr_cluster\nc0\t{fa}\t{bad}\n")
    r = subprocess.run([sys.executable, str(NEW), "--manifest", str(man), "--out", str(tmp_path / "o.tsv")],
                       capture_output=True, text=True)
    assert r.returncode != 0 and "not in this chunk" in r.stderr and not (tmp_path / "o.tsv").exists()
    man.write_text(f"chunk\tqueries\tnr_cluster\nc0\t{fa}\t{tmp_path / 'missing.tsv'}\n")
    r = subprocess.run([sys.executable, str(NEW), "--manifest", str(man), "--out", str(tmp_path / "o.tsv")],
                       capture_output=True, text=True)
    assert r.returncode != 0 and not (tmp_path / "o.tsv").exists()
