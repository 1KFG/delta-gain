#!/usr/bin/env python3
"""Map samples.csv ASMIDs to their annotated-protein FASTA in Fungi_BFD_runs's
input/pep/ (Stage 3 BFD-side input, HANDOFF.md item 1/3).

Fungi_BFD_runs/input/pep/*.proteins.fa is named by species+strain tag, not by
ASMID -- funannotate writes it using the exact same `makeSampleTag(SPECIES,
STRAIN)` convention as the annotation pipeline itself
(Fungi_BFD/nextflow/modules/common/utils.nf:58-114), so this ports that
Groovy logic verbatim rather than guessing a sanitization rule. IMPORTANT:
the tag is built from the `SPECIES` column (genus+species, no strain), not
`SPECIES_IN` (free-text, sometimes carries "sp."/"aff." qualifiers and
already-embedded strain text) -- confirmed from the callers
(funannotate.nf:110, BFD.nf:65,147, paralogoscope.nf:65,92,
setup_symlinks.nf:38), which all pass row.SPECIES.

Some tags collide across >1 ASMID (e.g. two accessions for the same
species+strain deposited separately) -- the pep filename alone can't say
which ASMID's proteins are actually in the file, so those rows are flagged
`ambiguous_tag` rather than guessed at.

Not every genome is annotated yet (54% had no matching pep file as of
2026-09-08) -- rows without a matching file are flagged `not_annotated`
rather than silently dropped, since "not yet annotated" is itself a relevant
bookkeeping fact (see HANDOFF.md item 1), not just a missing row.
"""
import argparse
import csv
import os
import re
import sys
from collections import Counter


def clean_strain(raw_strain):
    s = (raw_strain or "").strip()
    s = re.sub(r"""['"]""", "", s)
    s = s.split(";")[0].strip()
    s = s.replace(":", " ")
    s = re.sub(r"^\s*\*+", "", s)
    s = re.sub(r"\*+\s*$", "", s)
    s = re.sub(r"\s*\*+\s*", "-", s)
    return s.strip()


def make_sample_tag(raw_species, raw_strain):
    sp = (raw_species or "").strip()
    sp = re.sub(r"""['"]""", "", sp)
    st = clean_strain(raw_strain)
    tag = "_".join(part for part in (sp, st) if part)
    tag = re.sub(r"[\s/#\[\]?{}]+", "_", tag)
    return tag


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True, help="Fungi_BFD_runs/samples.csv")
    ap.add_argument("--pep-dir", required=True, help="Fungi_BFD_runs/input/pep")
    ap.add_argument("--n-test", type=int, default=0,
                     help="cap the number of matched rows used downstream "
                          "(0 = no cap); rows past the cap are still listed "
                          "in the manifest, marked matched_excluded_by_n_test, "
                          "so annotation coverage stays fully visible")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.samples, newline="") as fh:
        rows = list(csv.DictReader(fh))

    tags = [make_sample_tag(row["SPECIES"], row["STRAIN"]) for row in rows]
    tag_counts = Counter(tags)

    manifest_rows = []
    n_matched_used = 0
    for row, tag in zip(rows, tags):
        asmid = row["ASMID"].strip()
        pep_path = os.path.join(args.pep_dir, f"{tag}.proteins.fa")

        if tag_counts[tag] > 1:
            status = "ambiguous_tag"
        elif os.path.exists(pep_path):
            if args.n_test > 0 and n_matched_used >= args.n_test:
                status = "matched_excluded_by_n_test"
            else:
                status = "matched"
                n_matched_used += 1
        else:
            status = "not_annotated"
            pep_path = ""

        manifest_rows.append((asmid, tag, pep_path, status))

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["ASMID", "tag", "pep_path", "status"])
        w.writerows(manifest_rows)

    counts = Counter(r[3] for r in manifest_rows)
    print(
        f"build_pep_manifest: {len(manifest_rows)} samples.csv rows -> "
        f"matched={counts['matched']} "
        f"matched_excluded_by_n_test={counts['matched_excluded_by_n_test']} "
        f"ambiguous_tag={counts['ambiguous_tag']} "
        f"not_annotated={counts['not_annotated']}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
