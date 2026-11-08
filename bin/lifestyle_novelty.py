#!/usr/bin/env python3
"""Lifestyle correlates of novelty (FungalTraits 1.2 genus table, Polme et al. 2020).

Unit = genus (new genomes of a genus averaged), so strains of one genus do not inflate n. Per lifestyle: mean novelty,
and the within-order deviation (genus mean minus the mean of its order, orders with >=2 lifestyle classes) so clade
composition is controlled. Significance by permuting lifestyle labels among genera inside each order.
"""
import argparse
import os
import random
import sys
from collections import defaultdict

import duckdb


def collapse(p):
    p = (p or "").strip()
    if not p:
        return "unknown"
    if p.endswith("saprotroph"):
        return "saprotroph"
    return {"plant_pathogen": "plant_pathogen", "animal_parasite": "animal_assoc", "animal_endosymbiont": "animal_assoc",
            "ectomycorrhizal": "mycorrhizal", "arbuscular_mycorrhizal": "mycorrhizal", "ericoid_mycorrhizal": "mycorrhizal",
            "foliar_endophyte": "endophyte", "root_endophyte": "endophyte", "lichenized": "lichen",
            "mycoparasite": "mycoparasite", "unspecified": "unknown"}.get(p, "other")


def load(stratum, ft, min_genomes, strata=None):
    con = duckdb.connect()
    con.execute(f"""create table ft as select lower(trim(GENUS)) genus, primary_lifestyle pl, Growth_form_template gf
        from read_csv('{ft}', header=true, all_varchar=true, ignore_errors=true, strict_mode=false)""")
    where = ("where left(s.stratum, 1) in (" + ",".join("'%s'" % c for c in strata.split(",")) + ")") if strata else ""
    rows = con.execute(f"""
      select lower(trim(s.genus)) genus, any_value(s."ord") ord, any_value(s.stratum) stratum, count(*) n,
             avg(frac_95), avg(frac_70), avg(frac_50), any_value(ft.pl), any_value(ft.gf)
      from read_csv('{stratum}', delim='\t', header=true) s left join ft on ft.genus = lower(trim(s.genus))
      {where}
      group by lower(trim(s.genus)) having count(*) >= {min_genomes}""").fetchall()
    con.close()
    return rows


def dev(items, key, labels):
    """within-order deviation of `key` per label (orders holding >= 2 labels only)."""
    byo = defaultdict(list)
    for g, l in zip(items, labels):
        byo[g["ord"]].append((g[key], l))
    s, c = defaultdict(float), defaultdict(int)
    for v in byo.values():
        if len({l for _, l in v}) < 2:
            continue
        m = sum(x for x, _ in v) / len(v)
        for x, l in v:
            s[l] += x - m
            c[l] += 1
    return {l: s[l] / c[l] for l in s}, c


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--stratum", required=True)
    ap.add_argument("--traits", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-genomes", type=int, default=1)
    ap.add_argument("--strata", default=None, help="comma list of stratum digits to keep, e.g. 2,3")
    ap.add_argument("--perm", type=int, default=2000)
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    rows = load(a.stratum, a.traits, a.min_genomes, a.strata)
    G = [dict(genus=r[0], ord=r[1] or "NA", n=r[3], f95=r[4], f70=r[5], f50=r[6], life=collapse(r[7]),
              yeast=("yeast" in (r[8] or ""))) for r in rows]
    known = [g for g in G if g["life"] != "unknown"]
    print("genera", len(G), "with lifestyle", len(known), "genomes", sum(g["n"] for g in G),
          "covered", sum(g["n"] for g in known))

    out = open(os.path.join(a.outdir, "lifestyle_novelty.tsv"), "w")
    out.write("metric\tlifestyle\tgenera\tgenomes\tmean\twithin_order_dev\tn_dev\tperm_p\n")
    rnd = random.Random(7)
    byo = defaultdict(list)
    for i, g in enumerate(known):
        byo[g["ord"]].append(i)
    for key in ("f95", "f70", "f50"):
        obs, cnt = dev(known, key, [g["life"] for g in known])
        ge = {l: 0 for l in obs}
        for _ in range(a.perm):
            lab = [None] * len(known)
            for idx in byo.values():
                ls = [known[i]["life"] for i in idx]
                rnd.shuffle(ls)
                for i, l in zip(idx, ls):
                    lab[i] = l
            d, _c = dev(known, key, lab)
            for l in obs:
                if abs(d.get(l, 0)) >= abs(obs[l]):
                    ge[l] += 1
        for l in sorted(obs):
            sel = [g for g in known if g["life"] == l]
            out.write(f"{key}\t{l}\t{len(sel)}\t{sum(g['n'] for g in sel)}\t{sum(g[key] for g in sel) / len(sel):.4f}\t"
                      f"{obs[l]:.4f}\t{cnt[l]}\t{(ge[l] + 1) / (a.perm + 1):.4f}\n")
    obs, cnt = dev(G, "f50", ["yeast" if g["yeast"] else "non_yeast" for g in G])
    for l in sorted(obs):
        n = sum(1 for g in G if g["yeast"] == (l == "yeast"))
        out.write(f"f50\tgrowth:{l}\t{n}\t\t\t{obs[l]:.4f}\t{cnt[l]}\tNA\n")
    out.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
