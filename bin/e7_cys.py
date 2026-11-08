#!/usr/bin/env python3
"""E7 stage 2: cysteine counts for SSP candidate proteins by scanning each genome's BFD protein FASTA once (multiprocessing).

Reads cand.parquet (prefix, protein_id) and the pep manifest; writes cys.parquet (protein_id, cys, plen_seq).
"""
import argparse
import csv
import gzip
import multiprocessing as mp
import os
import sys
from collections import defaultdict

import duckdb

CAND = {}


def work(job):
    prefix, path, outdir = job
    ids = CAND.get(prefix)
    if not ids:
        return prefix, 0, 0
    got = 0
    opener = gzip.open if path.endswith(".gz") else open
    out = []
    cur = None
    seq = []
    try:
        with opener(path, "rt") as fh:
            for line in fh:
                if line[0] == ">":
                    if cur is not None:
                        s = "".join(seq)
                        out.append((cur, s.count("C"), len(s)))
                        got += 1
                    pid = line[1:].split(None, 1)[0]
                    cur = pid if pid in ids else None
                    seq = []
                elif cur is not None:
                    seq.append(line.strip())
            if cur is not None:
                s = "".join(seq)
                out.append((cur, s.count("C"), len(s)))
                got += 1
    except Exception as e:  # noqa: BLE001
        return prefix, -1, repr(e)[:100]
    with open(os.path.join(outdir, prefix + ".tsv"), "w") as fo:
        for pid, cc, n in out:
            fo.write(f"{pid}\t{cc}\t{n}\n")
    return prefix, got, len(ids)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--cand", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--workers", type=int, default=32)
    a = ap.parse_args(argv)
    shards = os.path.join(a.outdir, "cys_shards")
    os.makedirs(shards, exist_ok=True)
    con = duckdb.connect()
    for prefix, pid in con.execute(f"select prefix, protein_id from read_parquet('{a.cand}')").fetchall():
        CAND.setdefault(prefix, set()).add(pid)
    print("candidate proteins", sum(len(v) for v in CAND.values()), "genomes", len(CAND), flush=True)
    path_of = {}
    with open(a.manifest) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            path_of[r["ASMID"]] = r["pep_path"]
    jobs = []
    for asmid, prefix in con.execute(f"select asmid, prefix from read_parquet('{a.meta}')").fetchall():
        if prefix in CAND and asmid in path_of and not os.path.exists(os.path.join(shards, prefix + ".tsv")):
            jobs.append((prefix, path_of[asmid], shards))
    print("jobs", len(jobs), flush=True)
    bad = 0
    with mp.get_context("fork").Pool(a.workers) as pool:
        for i, (prefix, got, tot) in enumerate(pool.imap_unordered(work, jobs, chunksize=4)):
            if got < 0:
                bad += 1
                print("ERR", prefix, tot, flush=True)
            if i % 2000 == 0:
                print("done", i, flush=True)
    print("errors", bad, flush=True)
    con.execute(f"""copy (select column0 protein_id, column1 cys, column2 plen_seq
      from read_csv('{shards}/*.tsv', delim='\t', header=false, columns={{'column0':'VARCHAR','column1':'INTEGER','column2':'INTEGER'}}))
      to '{a.outdir}/cys.parquet' (format parquet)""")
    print("cys rows", con.execute(f"select count(*) from read_parquet('{a.outdir}/cys.parquet')").fetchone(), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
