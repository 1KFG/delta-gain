#!/usr/bin/env python3
"""E4: species accumulation of tier-50 clusters (baseline-first / new-first / fully random orderings).

Stage A (--stage pairs): distinct (genome, tier-50 cluster) pairs by expanding tier-95 members through the nested tiers
(95 -> 70 -> 50). Cluster ids are 64-bit hashes of the tier-50 representative name.
Stage B (--stage curves): units = species within group (baseline genomes of a species, new genomes of a species). A unit is
'N_known' if its species also has a baseline genome, else 'N_newsp'. For each permutation, a cluster is credited to the first unit
that contains it; credits give the marginal clusters per unit and the cumulative curve.
  order A: baseline units (random) then new units (random)   -> marginal gain of new units after the baseline is exhausted
  order B: new units first then baseline units                -> baseline depth vs new-unit novelty contrast
  order C: all units randomly interleaved                     -> overall species accumulation curve
"""
import argparse
import os
import sys

import numpy as np


def stage_pairs(a):
    import duckdb
    F = a.funnel_dir
    c = duckdb.connect()
    c.execute(f"pragma threads={a.threads}")
    c.execute(f"pragma memory_limit='{a.memory_limit}'")
    if a.tmp_dir:
        c.execute(f"pragma temp_directory='{a.tmp_dir}'")
    rd = lambda p: (f"read_csv('{p}', delim='\t', quote='', escape='', strict_mode=false, header=false, "
                    f"columns={{'rep':'VARCHAR','member':'VARCHAR'}})")
    c.execute(f"create table t50 as select rep r50, member m70 from {rd(F + '/tier_50_clusters.tsv')}")
    c.execute(f"create table j as select t.member m95, hash(t50.r50) h50 from {rd(F + '/tier_70_clusters.tsv')} t join t50 on t50.m70 = t.rep")
    print("j rows", c.execute("select count(*) from j").fetchone()[0], flush=True)
    c.execute(f"""create table pairs as
      select distinct split_part(t.member, '__', 2) asm, j.h50 from {rd(F + '/tier_95_clusters.tsv')} t join j on j.m95 = t.rep""")
    n = c.execute("select count(*), count(distinct asm), count(distinct h50) from pairs").fetchone()
    print("pairs", n, flush=True)
    c.execute(f"copy pairs to '{a.outdir}/pairs.parquet' (format parquet)")


def stage_curves(a):
    import duckdb
    c = duckdb.connect()
    c.execute(f"pragma threads={a.threads}")
    pr = c.execute(f"select asm, h50 from read_parquet('{a.outdir}/pairs.parquet')").fetchnumpy()
    asm_names = pr["asm"]
    h = pr["h50"].astype(np.uint64)
    uniq_asm, asm_idx = np.unique(asm_names, return_inverse=True)
    print("pairs", len(h), "genomes", len(uniq_asm), flush=True)
    # genome -> species, class
    meta = c.execute(f"""select cls.ASMID asmid, cls.class cl, coalesce(sp.species, cls.ASMID) species
        from read_csv('{a.classes}', delim='\t', header=true) cls
        left join (select ASMID, any_value(SPECIES) species from read_parquet('{a.tables_dir}/species.parquet') group by 1) sp on sp.ASMID = cls.ASMID""").fetchall()
    m = {r[0]: (r[1], r[2]) for r in meta}
    base_species = {sp for (cl, sp) in m.values() if cl == "baseline"}
    unit_key = {}
    unit_of_genome = np.full(len(uniq_asm), -1, dtype=np.int64)
    kinds = []
    for i, g in enumerate(uniq_asm):
        if g not in m:
            continue
        cl, sp = m[g]
        key = (cl, sp)
        if key not in unit_key:
            unit_key[key] = len(unit_key)
            kinds.append("B" if cl == "baseline" else ("N_known" if sp in base_species else "N_newsp"))
        unit_of_genome[i] = unit_key[key]
    kinds = np.array(kinds)
    nu = len(unit_key)
    print("units", nu, {k: int((kinds == k).sum()) for k in ("B", "N_known", "N_newsp")}, flush=True)
    unit_of_pair = unit_of_genome[asm_idx]
    keep = unit_of_pair >= 0
    h, unit_of_pair = h[keep], unit_of_pair[keep]
    order = np.argsort(h, kind="stable")
    h, unit_of_pair = h[order], unit_of_pair[order]
    starts = np.flatnonzero(np.r_[True, h[1:] != h[:-1]])
    ncl = len(starts)
    print("clusters", ncl, flush=True)
    rng = np.random.default_rng(a.seed)
    out_credit = open(os.path.join(a.outdir, "accum_credit.tsv"), "w")
    out_credit.write("order\tperm\tkind\tunits\tclusters_credited\tmean_per_unit\n")
    curves = {o: [] for o in "ABC"}
    isB = kinds == "B"
    for o in "ABC":
        for p in range(a.perms):
            if o == "C":
                perm = rng.permutation(nu)
            else:
                b = rng.permutation(np.flatnonzero(isB))
                n = rng.permutation(np.flatnonzero(~isB))
                perm = np.r_[b, n] if o == "A" else np.r_[n, b]
            rank = np.empty(nu, dtype=np.int64)
            rank[perm] = np.arange(nu)
            first = np.minimum.reduceat(rank[unit_of_pair], starts)
            cnt = np.bincount(first, minlength=nu)  # clusters first covered at position r
            curves[o].append(np.cumsum(cnt))
            credit_by_unit = np.zeros(nu, dtype=np.int64)
            credit_by_unit[perm] = cnt
            for k in ("B", "N_known", "N_newsp"):
                sel = kinds == k
                out_credit.write(f"{o}\t{p}\t{k}\t{int(sel.sum())}\t{int(credit_by_unit[sel].sum())}\t{credit_by_unit[sel].mean():.1f}\n")
        out_credit.flush()
    out_credit.close()
    # mean cumulative curve (position = number of units added) with sd, thinned
    pos = np.unique(np.r_[np.arange(0, nu, 25), nu - 1, int(isB.sum()) - 1, int(isB.sum())])
    pos = pos[(pos >= 0) & (pos < nu)]
    with open(os.path.join(a.outdir, "accum_curves.tsv"), "w") as fh:
        fh.write("order\tunits_added\tmean_clusters\tsd_clusters\n")
        for o in "ABC":
            arr = np.vstack(curves[o])
            for p_ in pos:
                fh.write(f"{o}\t{p_ + 1}\t{arr[:, p_].mean():.0f}\t{arr[:, p_].std():.0f}\n")
    with open(os.path.join(a.outdir, "accum_units.tsv"), "w") as fh:
        fh.write("kind\tunits\n")
        for k in ("B", "N_known", "N_newsp"):
            fh.write(f"{k}\t{int((kinds == k).sum())}\n")
        fh.write(f"total_clusters\t{ncl}\n")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["pairs", "curves"], required=True)
    ap.add_argument("--funnel-dir")
    ap.add_argument("--tables-dir")
    ap.add_argument("--classes")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--memory-limit", default="150GB")
    ap.add_argument("--tmp-dir", default=None)
    ap.add_argument("--perms", type=int, default=20)
    ap.add_argument("--seed", type=int, default=11)
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    (stage_pairs if a.stage == "pairs" else stage_curves)(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
