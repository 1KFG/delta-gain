#!/usr/bin/env python3
"""Concatenate the BFD proteins named `matched` in build_pep_manifest.py's
manifest into one combined FASTA for MMSEQS_LINCLUST, prefixing every
sequence ID with its ASMID so per-protein provenance survives dereplication
and the later DIAMOND search (a cluster representative's original ASMID
would otherwise be lost once IDs are just funannotate's per-genome locus
tags, which are not unique across genomes).

Dependency-free FASTA parsing (no biopython) -- input files are plain
funannotate protein FASTA, nothing exotic in the header/sequence format that
would need a real parser.
"""
import argparse
import csv
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
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-provenance", required=True)
    args = ap.parse_args()

    n_genomes = 0
    n_proteins = 0
    with open(args.manifest, newline="") as mf, \
         open(args.out_fasta, "w") as fasta_out, \
         open(args.out_provenance, "w", newline="") as prov_out:

        prov_w = csv.writer(prov_out, delimiter="\t")
        prov_w.writerow(["new_id", "asmid", "original_id"])

        for row in csv.DictReader(mf, delimiter="\t"):
            if row["status"] != "matched":
                continue
            asmid = row["ASMID"]
            n_genomes += 1
            for header, seq in iter_fasta(row["pep_path"]):
                orig_id = header.split()[0] if header else ""
                new_id = f"{asmid}__{orig_id}"
                fasta_out.write(f">{new_id}\n")
                for i in range(0, len(seq), 60):
                    fasta_out.write(seq[i:i + 60] + "\n")
                prov_w.writerow([new_id, asmid, orig_id])
                n_proteins += 1

    print(
        f"combine_bfd_proteins: {n_proteins} proteins from {n_genomes} "
        f"matched genomes -> {args.out_fasta}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
