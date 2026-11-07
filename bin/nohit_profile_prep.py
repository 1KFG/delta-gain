#!/usr/bin/env python3
"""Profile search prep: expand the 'no_hit' survivor clusters (sensitive sequence search found no baseline homolog) to all member
proteins, so MMseqs2 can build one profile per cluster and search it against baseline-known representatives.
Writes clusters.tsv (rep<TAB>member, rep included as its own member), members.faa, reps.txt."""
import argparse
import duckdb

ap = argparse.ArgumentParser()
ap.add_argument("--funnel-dir", required=True); ap.add_argument("--best", required=True)
ap.add_argument("--outdir", required=True); ap.add_argument("--threads", type=int, default=16)
ap.add_argument("--memory-limit", default="96GB"); ap.add_argument("--tmp-dir", default=None)
ap.add_argument("--max-clusters", type=int, default=0)
a = ap.parse_args()
F = a.funnel_dir
c = duckdb.connect()
c.execute(f"pragma threads={a.threads}"); c.execute(f"pragma memory_limit='{a.memory_limit}'")
if a.tmp_dir:
    c.execute(f"pragma temp_directory='{a.tmp_dir}'")
rd = lambda p: (f"read_csv('{p}', delim='\t', quote='', escape='', strict_mode=false, header=false, columns={{'rep':'VARCHAR','member':'VARCHAR'}})")
lim = f"limit {a.max_clusters}" if a.max_clusters else ""
c.execute(f"create table s as select rep from read_csv('{a.best}', delim='\t', header=true) where cls = 'no_hit' order by rep {lim}")
c.execute(f"create table l1 as select t.rep srep, t.member m from {rd(F + '/tier_50_clusters.tsv')} t semi join s on s.rep = t.rep")
c.execute(f"create table l2 as select l1.srep, t.member m from {rd(F + '/tier_70_clusters.tsv')} t join l1 on l1.m = t.rep")
c.execute(f"create table l3 as select l2.srep, t.member m from {rd(F + '/tier_95_clusters.tsv')} t join l2 on l2.m = t.rep where t.member like 'NEW__%'")
c.execute(f"copy l3 to '{a.outdir}/clusters.tsv' (delimiter '\t', header false)")
c.execute("create table ids as select distinct m from l3")
n = c.execute("select count(*) from ids").fetchone()[0]
print("no_hit clusters", c.execute("select count(*) from s").fetchone()[0], "member proteins", n, flush=True)
ids = {r[0] for r in c.execute("select m from ids").fetchall()}
c.close()
out = None
kept = 0
with open(f"{F}/funnel_input.faa") as fh, open(f"{a.outdir}/members.faa", "w") as fo:
    for line in fh:
        if line[0] == ">":
            keep = line[1:].split()[0] in ids
            if keep:
                kept += 1
        if keep:
            fo.write(line)
print("sequences written", kept, flush=True)
