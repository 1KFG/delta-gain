#!/usr/bin/env python3
"""F2 within-family clustering of extracted regions + baseline/new novelty summary.

Input  (from bin/family_extract.py):  <regions-dir>/<key>.faa  and  <key>.members.tsv
Output (per family, in --outdir):
  <key>.clusters.tsv      region_id, c90, c70, c50, c30   (representative id at each nested tier)
  <key>.tier_summary.tsv  one row per tier (see SUMMARY_COLS)
  tier_summary_all.tsv    all families concatenated

Nested tiers: 90 % by `mmseqs linclust`, later tiers by `mmseqs cluster` on the previous tier's
representatives (80 % bidirectional coverage, cov-mode 0, like the funnel and the metagenome-families paper).
A region is "novel at tier t" if its tier-t cluster contains NO baseline member.  Counted per member
region (never per representative).  Reported overall and split by whether the region's species has any
baseline genome ("within known species" = strict test: sequence not matched by baseline strains of the
same species; "new taxon" = species absent from the baseline).  Species-collapsed numbers count distinct
(species, cluster) pairs of new species, so heavily sequenced species do not dominate.
"""
import argparse
import csv
import os
import shutil
import subprocess
import sys

import duckdb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import family_census  # noqa: E402

SUMMARY_COLS = [
    "family", "tier", "regions_base", "regions_new", "clusters", "clusters_with_base",
    "clusters_new_only", "clusters_base_only", "new_regions_in_new_only", "frac_new_novel",
    "new_known_species", "novel_known_species", "frac_novel_known_species",
    "new_new_taxon", "novel_new_taxon", "frac_novel_new_taxon",
    "pairs_new_species", "pairs_novel", "frac_pairs_novel",
    "pairs_known_species", "pairs_novel_known_species", "frac_pairs_novel_known_species"]


def run(cmd):
    subprocess.run(cmd, shell=True, check=True)


def read_fasta(path):
    seqs, cur, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf)
                cur, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if cur is not None:
        seqs[cur] = "".join(buf)
    return seqs


def read_pairs(path):
    m = {}
    with open(path) as fh:
        for line in fh:
            rep, mem = line.rstrip("\n").split("\t")[:2]
            m[mem] = rep
    return m


def compose(first, later):
    """first: member->rep0 ; later: [rep_{k-1}->rep_k, ...] -> {member: [rep0, rep1, ...]}"""
    out = {}
    for mem, r0 in first.items():
        chain, cur = [r0], r0
        for step in later:
            cur = step.get(cur, cur)
            chain.append(cur)
        out[mem] = chain
    return out


def cluster_family(key, regions_dir, outdir, a):
    faa = os.path.join(regions_dir, key + ".faa")
    work = os.path.join(a.tmp_dir, key)
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    mm, th = a.mmseqs, a.threads
    run(f"{mm} createdb {faa} {work}/DB -v 1")
    ident = a.tiers[0]
    run(f"{mm} linclust {work}/DB {work}/c0 {work}/tmp0 --min-seq-id {ident} -c {a.cov0} --cov-mode 1 "
        f"--threads {th} -v 1")
    run(f"{mm} createtsv {work}/DB {work}/DB {work}/c0 {work}/c0.tsv --threads {th} -v 1")
    first = read_pairs(os.path.join(work, "c0.tsv"))
    seqs = read_fasta(faa)
    later, cur_ids = [], sorted(set(first.values()))
    for k, ident in enumerate(a.tiers[1:], start=1):
        rf = os.path.join(work, f"reps{k}.faa")
        with open(rf, "w") as fh:
            for i in cur_ids:
                fh.write(f">{i}\n{seqs[i]}\n")
        run(f"{mm} createdb {rf} {work}/R{k} -v 1")
        run(f"{mm} cluster {work}/R{k} {work}/c{k} {work}/tmp{k} --min-seq-id {ident} -c {a.cov} "
            f"--cov-mode 0 -s {a.sens} --threads {th} -v 1")
        run(f"{mm} createtsv {work}/R{k} {work}/R{k} {work}/c{k} {work}/c{k}.tsv --threads {th} -v 1")
        step = read_pairs(os.path.join(work, f"c{k}.tsv"))
        later.append(step)
        cur_ids = sorted(set(step.values()))
    chains = compose(first, later)
    labels = [f"c{int(round(t * 100))}" for t in a.tiers]
    asg = os.path.join(outdir, key + ".clusters.tsv")
    with open(asg, "w") as fh:
        fh.write("region_id\t" + "\t".join(labels) + "\n")
        for rid, ch in chains.items():
            fh.write(rid + "\t" + "\t".join(ch) + "\n")
    shutil.rmtree(work, ignore_errors=True)
    return asg, labels


def summarize(key, members_tsv, assign_tsv, labels, has_base):
    """has_base: {species: bool}. Returns list of rows (SUMMARY_COLS)."""
    con = duckdb.connect()
    con.execute("create table spt (species varchar, has_base boolean)")
    con.executemany("insert into spt values (?, ?)", list(has_base.items()))
    con.execute(f"""create table mem as select region_id, "class" cls, species from
        read_csv('{members_tsv}', delim='\t', header=true, all_varchar=true)""")
    con.execute(f"create table asg as select * from read_csv('{assign_tsv}', delim='\t', header=true, all_varchar=true)")
    rows = []
    for lab in labels:
        con.execute(f"""create or replace table m as
            select a.{lab} c, mem.cls, mem.species, coalesce(s.has_base, false) kb
            from asg a join mem using (region_id) left join spt s on s.species = mem.species""")
        con.execute("""create or replace table cl as
            select c, count(*) filter (cls = 'baseline') nb, count(*) filter (cls = 'new') nn from m group by c""")
        (rb, rn, ncl, nwb, nno, nbo, nin, nk, nok, nt, notx) = con.execute("""
            select (select count(*) from m where cls='baseline'), (select count(*) from m where cls='new'),
              (select count(*) from cl), (select count(*) from cl where nb>0),
              (select count(*) from cl where nb=0), (select count(*) from cl where nn=0),
              (select count(*) from m join cl using (c) where m.cls='new' and cl.nb=0),
              (select count(*) from m where cls='new' and kb),
              (select count(*) from m join cl using (c) where m.cls='new' and kb and cl.nb=0),
              (select count(*) from m where cls='new' and not kb),
              (select count(*) from m join cl using (c) where m.cls='new' and not kb and cl.nb=0)""").fetchone()
        (np_, npn, npk, npkn) = con.execute("""
            with p as (select distinct m.species, m.c, m.kb, cl.nb from m join cl using (c) where m.cls='new')
            select count(*), count(*) filter (nb=0), count(*) filter (kb), count(*) filter (kb and nb=0) from p"""
                                              ).fetchone()
        f = lambda x, y: (x / y if y else None)  # noqa: E731
        rows.append([key, lab, rb, rn, ncl, nwb, nno, nbo, nin, f(nin, rn), nk, nok, f(nok, nk),
                     nt, notx, f(notx, nt), np_, npn, f(npn, np_), npk, npkn, f(npkn, npk)])
    return rows


def fmt(v):
    return "NA" if v is None else (f"{v:.6g}" if isinstance(v, float) else str(v))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--regions-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tables-dir", required=True)
    ap.add_argument("--classes", required=True)
    ap.add_argument("--families", nargs="*", default=None, help="file keys; default: all in regions-dir")
    ap.add_argument("--tiers", type=float, nargs="+", default=[0.90, 0.70, 0.50, 0.30])
    ap.add_argument("--cov0", type=float, default=0.90, help="tier-0 coverage (target-covers, mode 1)")
    ap.add_argument("--cov", type=float, default=0.80, help="later tiers, bidirectional (mode 0)")
    ap.add_argument("--sens", default="6")
    ap.add_argument("--mmseqs", default="mmseqs")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--tmp-dir", default=None)
    ap.add_argument("--min-regions", type=int, default=20)
    ap.add_argument("--aggregate", action="store_true",
                    help="only concatenate existing <key>.tier_summary.tsv into tier_summary_all.tsv")
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    if a.aggregate:
        files = sorted(f for f in os.listdir(a.outdir) if f.endswith(".tier_summary.tsv"))
        with open(os.path.join(a.outdir, "tier_summary_all.tsv"), "w") as out:
            out.write("\t".join(SUMMARY_COLS) + "\n")
            for f in files:
                with open(os.path.join(a.outdir, f)) as fh:
                    next(fh)
                    out.write(fh.read())
        print("aggregated", len(files), "families")
        return 0
    a.tmp_dir = a.tmp_dir or os.path.join(a.outdir, "_tmp")
    os.makedirs(a.tmp_dir, exist_ok=True)

    con = duckdb.connect()
    a.max_genomes = 0
    family_census.build_genomes(con, a)
    has_base = dict(con.execute("select species, has_base from sp").fetchall())
    con.close()

    keys = a.families or sorted(f[:-4] for f in os.listdir(a.regions_dir) if f.endswith(".faa"))
    allrows = []
    for key in keys:
        mt = os.path.join(a.regions_dir, key + ".members.tsv")
        n = sum(1 for _ in open(mt)) - 1
        if n < a.min_regions:
            print(f"skip {key}: {n} regions", flush=True)
            continue
        print(f"cluster {key}: {n} regions", flush=True)
        asg, labels = cluster_family(key, a.regions_dir, a.outdir, a)
        rows = summarize(key, mt, asg, labels, has_base)
        with open(os.path.join(a.outdir, key + ".tier_summary.tsv"), "w") as fh:
            fh.write("\t".join(SUMMARY_COLS) + "\n")
            for r in rows:
                fh.write("\t".join(fmt(v) for v in r) + "\n")
        allrows += rows
        for r in rows:
            print("  ", r[1], "clusters", r[4], "new-only", r[6], "frac_new_novel", fmt(r[9]),
                  "known-sp", fmt(r[12]), flush=True)
    if a.families is None:
        with open(os.path.join(a.outdir, "tier_summary_all.tsv"), "w") as fh:
            fh.write("\t".join(SUMMARY_COLS) + "\n")
            for r in allrows:
                fh.write("\t".join(fmt(v) for v in r) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
