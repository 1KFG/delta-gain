#!/usr/bin/env python3
"""Alignment-based catalytic-residue check for GH18 (PF00704) and S08 (PF00082) regions.

Replaces the regex-only approximation of family_motif.py: every region is aligned to the Pfam profile HMM with hmmalign;
the catalytic positions are located once on the HMM consensus (motif regex on the consensus) and read from the
match columns of each region. For each group (baseline / new, known cluster / new, novel cluster at --tier):
  spans  = fraction of regions with a residue (not a gap) at every catalytic position   (fragment / deletion control)
  intact = fraction with the canonical residue at every catalytic position
  intact_given_span = intact among regions that span all sites   (pseudo-enzyme signal free of truncation)
"""
import argparse
import os
import re
import subprocess
import sys

import duckdb

# family key -> (hmm file, consensus motif regex, {name: (offset in motif match, required residues)})
SPEC = {
    "cazy_family_GH18": ("PF00704.hmm", r"[LIVMF]D[GA][LIVMF]D[LIVMF][DE][WF]E",
                         {"D1": (1, "D"), "D4": (4, "D"), "D6": (6, "DE"), "E8": (8, "E")}),
    "merops_family_S08": ("PF00082.hmm", None, {}),
}
# S08 has three separate motifs on the consensus
S08_MOTIFS = [("Asp", r"[LIVMF][LIVMF]D[ST]G", 2, "D"), ("His", r"H[GA][ST]H", 3, "H"), ("Ser", r"G[TS]S[MA][AT][ASTG]P", 2, "S")]


def run(cmd, **kw):
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, **kw)


def consensus(hmm):
    out = subprocess.run(["hmmemit", "-c", hmm], check=True, capture_output=True, text=True).stdout
    return "".join(l.strip() for l in out.split("\n") if not l.startswith(">")).upper()


def sites(fam, cons):
    """-> list of (name, hmm_node_1based, allowed residues)."""
    out = []
    if fam == "cazy_family_GH18":
        m = re.search(SPEC[fam][1], cons)
        assert m, "GH18 motif not found on consensus"
        for n, (off, res) in SPEC[fam][2].items():
            out.append((n, m.start() + off + 1, res))
    else:
        for n, rx, off, res in S08_MOTIFS:
            m = re.search(rx, cons)
            assert m, f"S08 {n} motif not found on consensus"
            out.append((n, m.start() + off + 1, res))
    return out


def read_pfam_aln(path, cols_needed):
    """Parse hmmalign --outformat pfam (one block). cols_needed: list of 1-based HMM nodes. -> {name: residues}"""
    rf = None
    seqs = {}
    with open(path) as fh:
        for l in fh:
            if l.startswith("#=GC RF"):
                rf = l.split(None, 2)[2].strip()
            elif l and l[0] not in "#/\n":
                p = l.split()
                if len(p) == 2:
                    seqs[p[0]] = p[1]
    match_idx = [i for i, ch in enumerate(rf) if ch != "."]  # alignment column of each HMM node (RF 'x' marks match)
    cols = [match_idx[k - 1] for k in cols_needed]
    return {n: "".join(s[c] for c in cols) for n, s in seqs.items()}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--regions-dir", required=True)
    ap.add_argument("--clusters-dir", required=True)
    ap.add_argument("--hmm-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tier", default="c70")
    ap.add_argument("--families", nargs="+", default=list(SPEC))
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    con = duckdb.connect()
    summ = open(os.path.join(a.outdir, "catalytic_summary.tsv"), "w")
    summ.write("family\tgroup\tregions\tspans\tintact\tintact_given_span\t" + "per_site_identity\n")
    for fam in a.families:
        hmm = os.path.join(a.hmm_dir, SPEC[fam][0])
        cons = consensus(hmm)
        st = sites(fam, cons)
        print(fam, "sites on consensus:", [(n, k, cons[k - 1]) for n, k, _ in st], flush=True)
        aln = os.path.join(a.outdir, fam + ".aln.pfam")
        if not os.path.exists(aln):
            run(["hmmalign", "--amino", "--trim", "--outformat", "pfam", "-o", aln, hmm,
                 os.path.join(a.regions_dir, fam + ".faa")])
        res = read_pfam_aln(aln, [k for _, k, _ in st])
        rows = []
        for rid, r in res.items():
            span = all(c not in "-." for c in r)
            ok = all(c in allowed for c, (_, _, allowed) in zip(r, st))
            rows.append((rid, r, int(span), int(ok)))
        con.execute("create or replace table cat(region_id varchar, res varchar, span int, ok int)")
        con.executemany("insert into cat values (?,?,?,?)", rows)
        con.execute(f"""create or replace table j as select m.region_id, replace(m."class", chr(13), '') cls, k.{a.tier} cid
          from read_csv('{a.regions_dir}/{fam}.members.tsv', delim='\t', header=true, all_varchar=true) m
          join read_csv('{a.clusters_dir}/{fam}.clusters.tsv', delim='\t', header=true, all_varchar=true) k using (region_id)""")
        con.execute("create or replace table cl as select cid, max((cls='baseline')::int) hb from j group by 1")
        con.execute("""create or replace table g as select j.region_id, cat.res, cat.span, cat.ok,
            case when j.cls='baseline' then 'baseline' when cl.hb=0 then 'new_novel' else 'new_known' end grp
            from j join cl using (cid) join cat using (region_id)""")
        for grp, n, sp, ok, ig in con.execute("""select grp, count(*), avg(span), avg(ok),
              sum(ok)/nullif(sum(span), 0) from g group by 1 order by 1""").fetchall():
            per = []
            for i, (nm, _, allowed) in enumerate(st):
                f = con.execute(f"select avg((substr(res, {i + 1}, 1) in ({','.join(repr(c) for c in allowed)}))::int) "
                                f"from g where grp = '{grp}' and span = 1").fetchone()[0]
                per.append(f"{nm}={f:.3f}")
            summ.write(f"{fam}\t{grp}\t{n}\t{sp:.4f}\t{ok:.4f}\t{ig:.4f}\t{' '.join(per)}\n")
        summ.flush()
        # most common non-canonical residue strings among new_novel regions that span all sites
        con.execute(f"""copy (select grp, res, count(*) n from g where span = 1 and ok = 0 group by 1, 2
                         order by 3 desc limit 60) to '{a.outdir}/{fam}.noncanonical.tsv' (delimiter '\t', header)""")
    summ.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
