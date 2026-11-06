#!/usr/bin/env python3
"""Write ONE FASTA holding the baseline proteins (NCBI, prior-annotated genomes)
and the query proteins (BFD proteins of genomes with no prior annotation).

IDs:  BASE__<ncbi folder>__<protein id>      NEW__<ASMID>__<original id>
Provenance TSV: id, class (BASE|NEW), genome.

Proteins shorter than --min-len are dropped (default 30 aa); duplicate IDs
within one genome are dropped and counted.
"""
import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from combine_bfd_proteins import iter_fasta  # noqa: E402


def emit(fa, prov, cls, genome, path, min_len, stats):
    seen = set()
    for header, seq in iter_fasta(path):
        pid = header.split()[0] if header else ""
        seq = seq.rstrip("*")
        if not pid or len(seq) < min_len:
            stats["short"] += 1
            continue
        if pid in seen:
            stats["dup_id"] += 1
            continue
        seen.add(pid)
        new_id = f"{cls}__{genome}__{pid}"
        fa.write(f">{new_id}\n{seq}\n")
        prov.write(f"{new_id}\t{cls}\t{genome}\n")
        stats[cls] += 1


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--status", required=True)
    ap.add_argument("--asm-dir", required=True)
    ap.add_argument("--new-manifest", required=True)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-prov", required=True)
    ap.add_argument("--min-len", type=int, default=30)
    a = ap.parse_args(argv)

    stats = {"BASE": 0, "NEW": 0, "short": 0, "dup_id": 0, "base_genomes": 0,
             "new_genomes": 0, "base_missing_file": 0}
    with open(a.out_fasta, "w") as fa, open(a.out_prov, "w") as prov:
        prov.write("id\tclass\tgenome\n")
        with open(a.status, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["status"] != "baseline_ok" or r["protein_file_ok"] != "1":
                    continue
                p = os.path.join(a.asm_dir, r["folder"].strip(), r["folder"].strip() + "_protein.faa.gz")
                if not os.path.exists(p):
                    stats["base_missing_file"] += 1
                    continue
                stats["base_genomes"] += 1
                emit(fa, prov, "BASE", r["folder"].strip(), p, a.min_len, stats)
        with open(a.new_manifest, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["status"] != "matched":
                    continue
                stats["new_genomes"] += 1
                emit(fa, prov, "NEW", r["ASMID"], r["pep_path"], a.min_len, stats)
    print("build_funnel_input: " + ", ".join(f"{k}={v}" for k, v in stats.items()),
          file=sys.stderr)


if __name__ == "__main__":
    main()
