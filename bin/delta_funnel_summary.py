#!/usr/bin/env python3
"""Summarise a nested MMseqs2 identity funnel over baseline + new proteins.

Tier 0 clusters ALL proteins (BASE_ and NEW_); each later tier clusters the
representatives of the previous tier.  A cluster is 'known' if any BASE_
protein sits in it at that tier (propagated through the nesting).  A NEW_
protein is novel at tier k when its tier-k cluster holds no BASE_ protein.
Counting is per protein (member-weighted), not per representative.

Outputs (in --outdir)
  tier_summary.tsv       tier label, new proteins, novel proteins, fraction
  genome_tier.tsv        per new genome: n proteins, novel count and fraction per tier
  survivors_clusters.tsv last-tier clusters with no BASE_ member:
                         rep, n_new_proteins, n_genomes
"""
import argparse
import os
import sys

import duckdb


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--prov", required=True)
    ap.add_argument("--tier-clusters", nargs="+", required=True,
                    help="rep<TAB>member TSVs, tier 0 first")
    ap.add_argument("--tier-labels", nargs="+", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--memory-limit", default="64GB")
    ap.add_argument("--tmp-dir", default=None)
    a = ap.parse_args(argv)
    if len(a.tier_clusters) != len(a.tier_labels):
        sys.exit("need one label per tier cluster file")
    os.makedirs(a.outdir, exist_ok=True)

    c = duckdb.connect()
    c.execute(f"PRAGMA threads={a.threads}")
    c.execute(f"PRAGMA memory_limit='{a.memory_limit}'")
    if a.tmp_dir:
        c.execute(f"PRAGMA temp_directory='{a.tmp_dir}'")

    c.execute(f"""create table prov as select id, class cls, genome from
        read_csv('{a.prov}', delim='\t', header=true, columns={{'id':'VARCHAR','class':'VARCHAR','genome':'VARCHAR'}})""")

    def read_tier(path):
        return (f"read_csv('{path}', delim='\t', header=false, "
                f"columns={{'rep':'VARCHAR','member':'VARCHAR'}})")

    # tier 0: known flag per rep from BASE members; np = NEW proteins -> current rep
    c.execute(f"""create table flag as select t.rep, bool_or(p.cls='BASE') known
                  from {read_tier(a.tier_clusters[0])} t join prov p on p.id=t.member group by 1""")
    c.execute(f"""create table np as select p.id, p.genome, t.rep cur
                  from {read_tier(a.tier_clusters[0])} t join prov p on p.id=t.member
                  where p.cls='NEW'""")
    results = []  # (label, per-genome table name)

    def record(k):
        c.execute(f"""create table g{k} as select np.genome, count(*) n,
                      sum((not f.known)::int) novel from np join flag f on f.rep=np.cur
                      group by 1""")

    record(0)
    for k in range(1, len(a.tier_clusters)):
        tc = read_tier(a.tier_clusters[k])
        c.execute(f"""create or replace table flag as select t.rep, bool_or(f.known) known
                      from {tc} t join flag f on f.rep=t.member group by 1""")
        c.execute(f"""create or replace table np as select np.id, np.genome, t.rep cur
                      from np join {tc} t on t.member=np.cur""")
        record(k)

    labels = a.tier_labels
    # tier summary
    with open(os.path.join(a.outdir, "tier_summary.tsv"), "w") as fh:
        fh.write("tier\tnew_proteins\tnovel_proteins\tfraction_novel\n")
        for k, lab in enumerate(labels):
            n, nov = c.execute(f"select sum(n), sum(novel) from g{k}").fetchone()
            fh.write(f"{lab}\t{n}\t{nov}\t{(nov / n) if n else float('nan'):.6f}\n")
    # genome table
    sel = ", ".join(f"g{k}.novel novel_{lab}, g{k}.novel*1.0/g0.n frac_{lab}"
                    for k, lab in enumerate(labels))
    joins = " ".join(f"join g{k} using(genome)" for k in range(1, len(labels)))
    c.execute(f"""copy (select g0.genome, g0.n new_proteins, {sel}
                  from g0 {joins} order by g0.genome)
                  to '{os.path.join(a.outdir, 'genome_tier.tsv')}' (delimiter '\t', header)""")
    # survivors at last tier
    last = len(labels) - 1
    c.execute(f"""copy (select f.rep, count(*) n_new_proteins, count(distinct np.genome) n_genomes
                  from np join flag f on f.rep=np.cur where not f.known group by 1
                  order by n_new_proteins desc)
                  to '{os.path.join(a.outdir, 'survivors_clusters.tsv')}' (delimiter '\t', header)""")
    print(f"delta_funnel_summary: wrote {a.outdir} ({len(labels)} tiers, last={labels[last]})",
          file=sys.stderr)


if __name__ == "__main__":
    main()
