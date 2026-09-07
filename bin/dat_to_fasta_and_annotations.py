#!/usr/bin/env python3
"""Convert a UniProtKB .dat(.gz) flatfile to FASTA + a cross-reference table.

UniProt's taxonomic-division bulk files (uniprot_{sprot,trembl}_fungi.dat.gz)
are only distributed as .dat/.xml, not FASTA -- this does the one-time local
conversion. Uses Bio.SwissProt.parse (not the plainer SeqIO "swiss" parser)
specifically because it preserves the full cross_references structure: every
entry's DR lines (InterPro, Pfam, GO, PROSITE, ...), which SeqIO's "swiss"
format collapses into dbxrefs and mostly discards. For TrEMBL entries these
InterPro/Pfam/GO calls come from UniProt's own automatic annotation pipeline
(InterProScan-based), so they're free, already-computed domain/GO evidence
for whatever a BFD protein hits against -- kept here rather than thrown away.

Output:
  <prefix>.fasta            standard header: ">ACCESSION description"
  <prefix>.xrefs.tsv        accession, db, id, extra1, extra2  (one row per
                             DR line; db is e.g. "InterPro", "Pfam", "GO")
"""
import argparse
import gzip
import sys

from Bio import SwissProt


def openmaybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dat", required=True)
    ap.add_argument("--prefix", required=True)
    args = ap.parse_args()

    n = 0
    with openmaybe_gz(args.dat) as fh, \
         open(f"{args.prefix}.fasta", "w") as fasta_out, \
         open(f"{args.prefix}.xrefs.tsv", "w") as xref_out:

        xref_out.write("accession\tdb\tid\textra1\textra2\n")

        for record in SwissProt.parse(fh):
            n += 1
            accession = record.accessions[0]
            fasta_out.write(f">{accession} {record.description}\n")
            seq = record.sequence
            for i in range(0, len(seq), 60):
                fasta_out.write(seq[i : i + 60] + "\n")

            for xref in record.cross_references:
                # xref is a tuple like ('InterPro', 'IPR000001', 'Kringle', '1')
                db = xref[0]
                rest = list(xref[1:]) + ["", ""]
                xref_out.write(f"{accession}\t{db}\t{rest[0]}\t{rest[1]}\t{rest[2]}\n")

    print(f"dat_to_fasta_and_annotations: converted {n} records from {args.dat}", file=sys.stderr)


if __name__ == "__main__":
    main()
