#!/usr/bin/env python3
"""Summarise E9: new-minus-null novel fraction (F3, genome-level) per family x variant x tier, at the larger matched genome count.

Writes results/e9/e9_summary.tsv and prints a markdown table. A 'flip' is a variant whose diff has the opposite sign to the default
variant's (|default diff| > 0.003 and |variant diff| > 0.003), a 'size change' a variant diff differing from the default by > 50 % and > 0.01."""
import os
import sys

import duckdb

E = "/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/results/e9"
VARS = ["env_cov80_m0", "env_cov60_m0", "env_cov90_m0", "env_cov80_m1", "ali_cov80_m0"]
TIERS = ["c70", "c50", "c30"]
c = duckdb.connect()
rows = []
for v in VARS:
    f = f"{E}/null_{v}/null_all.tsv"
    if not os.path.exists(f):
        print("missing", f, file=sys.stderr)
        continue
    q = c.execute(f"""
      select '{v}' variant, family, tier, try_cast(genomes_H as int) gH, try_cast(stat_new as double) new, try_cast(stat_null as double) nul,
             try_cast(diff_mean as double) diff, try_cast(diff_lo as double) lo, try_cast(diff_hi as double) hi
      from read_csv('{f}', delim=chr(9), header=true, all_varchar=true)
      where variant = 'genome' and tier in ({','.join(repr(t) for t in TIERS)})""").fetchall()
    rows += q
c.execute("create table r(variant varchar, family varchar, tier varchar, gH int, new double, nul double, diff double, lo double, hi double)")
c.executemany("insert into r values (?,?,?,?,?,?,?,?,?)", rows)
c.execute("""create table big as select r.* from r join (select family, variant, max(gH) m from r group by 1, 2) x
             on x.family = r.family and x.variant = r.variant and x.m = r.gH""")
c.execute(f"copy (select * from big order by family, tier, variant) to '{E}/e9_summary.tsv' (delimiter '\t', header)")
fams = [r[0] for r in c.execute("select distinct family from big order by 1").fetchall()]
print("| family | tier | " + " | ".join(VARS) + " | flips | size changes |")
print("|---|---|" + "---|" * (len(VARS) + 2))
for fam in fams:
    for t in TIERS:
        cells, d = [], {}
        for v in VARS:
            r = c.execute("select diff, lo, hi, gH from big where family=? and tier=? and variant=?", [fam, t, v]).fetchone()
            if r is None or r[0] is None:
                cells.append("NA")
                continue
            d[v] = r[0]
            cells.append(f"{r[0]:+.3f} [{r[1]:+.3f},{r[2]:+.3f}]")
        base = d.get("env_cov80_m0")
        flips = sizes = 0
        if base is not None:
            for v, x in d.items():
                if v == "env_cov80_m0":
                    continue
                if abs(base) > 0.003 and abs(x) > 0.003 and base * x < 0:
                    flips += 1
                if abs(x - base) > 0.01 and abs(x - base) > 0.5 * abs(base):
                    sizes += 1
        print(f"| {fam} | {t} | " + " | ".join(cells) + f" | {flips} | {sizes} |")
print("genomes_H used:", c.execute("select family, max(gH) from big group by 1").fetchall())
