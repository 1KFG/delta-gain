#!/usr/bin/env python3
"""Split the BFD pep manifest by NCBI prior-annotation status.

Inputs
  --pep-manifest  BUILD_PEP_MANIFEST output (ASMID, tag, pep_path, status)
  --status        prior_annotation_status.tsv (accession, folder, status, ...)
  --sibling-flags optional file, one NCBI accession per line: strains for which
                  a DIFFERENT assembly of the same strain carries NCBI
                  annotation (dropped by the project's species+strain dedup).

Classes written to --out-classes
  baseline            NCBI annotation present and protein file on disk
  new                 no prior NCBI annotation (BFD is the only annotation)
  new_sibling_flagged would be 'new' but a same-strain assembly is annotated
                      (excluded from the new set unless --include-sibling-flagged)
  annotated_pending   NCBI-annotated but proteins not on disk yet (excluded
                      from both sets so it can never be mis-counted as new)
  not_in_ncbi_table   ASMID not found in the status table (excluded, reported)

--out-new-manifest is the pep manifest with status 'matched' kept only for
class 'new' genomes that were 'matched' in the input, so that
combine_bfd_proteins.py / build_funnel_input.py read exactly the new set.
"""
import argparse
import csv
import re
import sys
from collections import Counter

ACC = re.compile(r"^(GC[AF]_\d+\.\d+)")


def load_status(path):
    by_folder, by_acc = {}, {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["folder"]:
                by_folder[r["folder"]] = r
            by_acc[r["accession"]] = r
    return by_folder, by_acc


def classify(rec, siblings, include_sibling):
    if rec is None:
        return "not_in_ncbi_table"
    st = rec["status"]
    if st == "baseline_ok":
        return "baseline"
    if st.startswith("annotated"):
        return "annotated_pending"
    # no_prior_annotation_*
    if rec["accession"] in siblings and not include_sibling:
        return "new_sibling_flagged"
    return "new"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--pep-manifest", required=True)
    ap.add_argument("--status", required=True)
    ap.add_argument("--sibling-flags")
    ap.add_argument("--include-sibling-flagged", action="store_true")
    ap.add_argument("--out-classes", required=True)
    ap.add_argument("--out-new-manifest", required=True)
    a = ap.parse_args(argv)

    by_folder, by_acc = load_status(a.status)
    siblings = set()
    if a.sibling_flags:
        siblings = {l.strip() for l in open(a.sibling_flags) if l.strip()}

    counts = Counter()
    with open(a.pep_manifest, newline="") as mf, \
         open(a.out_classes, "w", newline="") as oc, \
         open(a.out_new_manifest, "w", newline="") as on:
        rd = csv.DictReader(mf, delimiter="\t")
        wc = csv.writer(oc, delimiter="\t")
        wc.writerow(["ASMID", "accession", "class", "ncbi_status", "pep_status"])
        wn = csv.DictWriter(on, fieldnames=rd.fieldnames, delimiter="\t")
        wn.writeheader()
        for row in rd:
            asm = row["ASMID"]
            rec = by_folder.get(asm)
            if rec is None:
                m = ACC.match(asm)
                rec = by_acc.get(m.group(1)) if m else None
            cls = classify(rec, siblings, a.include_sibling_flagged)
            counts[(cls, row["status"])] += 1
            wc.writerow([asm, rec["accession"] if rec else "", cls,
                         rec["status"] if rec else "", row["status"]])
            out = dict(row)
            if not (cls == "new" and row["status"] == "matched"):
                out["status"] = "excluded_" + cls if row["status"] == "matched" else row["status"]
            wn.writerow(out)
    for (cls, pst), n in sorted(counts.items()):
        print(f"classify: {cls:22s} pep_status={pst:24s} {n}", file=sys.stderr)


if __name__ == "__main__":
    main()
