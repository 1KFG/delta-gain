#!/usr/bin/env python3
"""Extract cluster-representative sequences from a `<representative>\\t<member>`
cluster-membership table (DIAMOND `cluster --out` or MMseqs2 `createtsv`
format) plus the original combined FASTA those IDs came from.

Written for DIAMOND_CLUSTER rather than using `diamond cluster --reps`
directly: real bug found 2026-09-08 -- combining `--reps` and `--out` in the
same `diamond cluster` invocation (DIAMOND 2.2.6) crashes ("Error opening
file .../input.tsv") even on trivial inputs, while `--out` alone works fine.
Deriving the representative FASTA from `--out`'s table ourselves sidesteps
the bug entirely.
"""
import argparse
import sys


def iter_fasta(path):
    header = None
    seq_lines = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        yield header, "".join(seq_lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clusters", required=True, help="rep<TAB>member table, no header")
    ap.add_argument("--fasta", required=True, help="original combined FASTA")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rep_ids = set()
    with open(args.clusters) as fh:
        for line in fh:
            if not line.strip():
                continue
            rep_ids.add(line.rstrip("\n").split("\t")[0])

    n_written = 0
    with open(args.out, "w") as out_fh:
        for header, seq in iter_fasta(args.fasta):
            seq_id = header.split()[0] if header else ""
            if seq_id in rep_ids:
                out_fh.write(f">{header}\n")
                for i in range(0, len(seq), 60):
                    out_fh.write(seq[i:i + 60] + "\n")
                n_written += 1

    print(
        f"extract_cluster_reps: {n_written}/{len(rep_ids)} representative "
        f"sequences written to {args.out}",
        file=sys.stderr,
    )
    if n_written != len(rep_ids):
        print(
            "WARNING: representative count mismatch -- some rep IDs in "
            "the cluster table were not found in the FASTA",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
