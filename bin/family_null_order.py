#!/usr/bin/env python3
"""F3b leave-one-order-out null for within-family novelty (docs/FAMILY_DIVERSITY_PLAN.md, stage F3b).

F3 compares new vs baseline genomes only inside orders that have both. It cannot say what novelty is
expected when a whole taxon is absent from the reference. F3b answers that:

  For each order o: reference R_o = all baseline genomes NOT in order o (so o is entirely unseen).
    group `base_out`: K baseline genomes drawn from o (needs >= --min-base baseline genomes in o)
    group `new`     : K new genomes drawn from o (needs >= --min-new new genomes in o)
  stat = fraction of the drawn genomes' regions (genome variant) or distinct (species, cluster) pairs
  (species variant) whose tier-t cluster has no member in R_o.
Fixed K genomes per draw makes stats comparable across orders. Draws are generated once and reused for
every family. Orders with nb == 0 are the real "new-only" orders: their `new` stat is compared with the
`base_out` distribution (how novel a baseline order looks when it is held out of the baseline).
`--aggregate` writes order_null_summary.tsv (per family/tier/variant) from the per-family files.
"""
import argparse
import os
import sys

import duckdb
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # python -I drops the script dir
import family_null as fn  # noqa: E402

COLS = ["family", "tier", "variant", "order", "group", "n_base", "n_new", "k", "draws",
        "stat_mean", "stat_lo", "stat_hi", "cls"]


def make_order_draws(g_cls, g_ord, g_sp, K, n_draws, collapse, rng, min_base, min_new, g_cl=None):
    """-> list of (order_idx, group, k, [index arrays])."""
    out = []
    for o in np.unique(g_ord):
        b = np.where((g_ord == o) & (g_cls == 0))[0]
        n = np.where((g_ord == o) & (g_cls == 1))[0]
        for group, pool, need in (("base_out", b, min_base), ("new", n, min_new)):
            if len(pool) < need or len(pool) == 0:
                continue
            avail = len(np.unique(g_sp[pool])) if collapse else len(pool)
            k = min(K, avail)
            ds = []
            for _ in range(n_draws):
                p = fn.one_per_species(rng, pool, g_sp) if collapse else pool
                ds.append(rng.choice(p, k, replace=False))
            out.append((int(o), group, k, ds))
    # baseline-less orders are too small (a few genomes each) for K draws: pool them as one pseudo-order
    # (reference = the whole baseline; sentinel order index -1)
    base_ord = set(g_ord[g_cls == 0].tolist())
    n_only = [i for i in np.where(g_cls == 1)[0] if g_ord[i] not in base_ord]
    pool = np.array(n_only, dtype=int)
    avail = len(np.unique(g_sp[pool])) if collapse and len(pool) else len(pool)
    if len(pool) >= min_new and avail:
        k = min(K, avail)
        ds = []
        for _ in range(n_draws):
            p = fn.one_per_species(rng, pool, g_sp) if collapse else pool
            ds.append(rng.choice(p, k, replace=False))
        out.append((-1, "new_only_pool", k, ds))
    # F3c: baseline-less orders grouped by taxonomic class (reference = whole baseline; sentinel -(class index + 2))
    if g_cl is not None:
        base_cl = set(g_cl[g_cls == 0].tolist())
        for cl in sorted(set(g_cl[pool].tolist()) if len(pool) else []):
            if cl not in base_cl:
                continue
            pc = pool[g_cl[pool] == cl]
            avail = len(np.unique(g_sp[pc])) if collapse else len(pc)
            if len(pc) >= min_new and avail:
                k = min(K, avail)
                ds = []
                for _ in range(n_draws):
                    p = fn.one_per_species(rng, pc, g_sp) if collapse else pc
                    ds.append(rng.choice(p, k, replace=False))
                out.append((-int(cl) - 2, "new_only_class", k, ds))
    return out


def eval_order_family(cl, g_r, sp_g, g_cls, g_ord, order_draws, collapse):
    """Returns {(order, group): array of per-draw stats}."""
    ncl = int(cl.max()) + 1
    ng = len(g_cls)
    srt = np.argsort(g_r, kind="stable")
    off = np.searchsorted(g_r[srt], np.arange(ng + 1))
    sp_r = sp_g[g_r]
    base_r = g_cls[g_r] == 0
    ord_r = g_ord[g_r]
    res = {}
    cache = {}
    for o, group, k, ds in order_draws:
        if o not in cache:
            present = np.zeros(ncl, bool)
            present[cl[base_r & (ord_r != o)]] = True
            cache = {o: present}  # draws are grouped by order: keep only the current one
        nov = ~cache[o][cl]
        vals = []
        for idx in ds:
            ridx = np.concatenate([srt[off[g]:off[g + 1]] for g in idx]) if len(idx) else np.array([], int)
            if collapse:
                pair, first = np.unique(sp_r[ridx].astype(np.int64) * ncl + cl[ridx], return_index=True)
                vals.append(nov[ridx[first]].mean() if len(first) else np.nan)
            else:
                vals.append(nov[ridx].mean() if len(ridx) else np.nan)
        res[(o, group)] = np.array(vals)
    return res


def run_order_null(g_rows, fam_inputs, K, n_draws, seed, min_base, min_new):
    asm, g_cls, g_ord, g_sp, n_sp = fn.genome_arrays(g_rows)
    ord_names = np.unique([r[2] for r in g_rows])
    cl_names, g_cl = np.unique([r[4] if len(r) > 4 else "NA" for r in g_rows], return_inverse=True)
    ord_cl = {}
    for oi, cli in zip(g_ord, g_cl):
        ord_cl.setdefault(int(oi), int(cli))
    out = []
    base_ord = set(g_ord[g_cls == 0].tolist())
    n_pool = int(sum(1 for i in np.where(g_cls == 1)[0] if g_ord[i] not in base_ord))
    for variant in ("genome", "species"):
        collapse = variant == "species"
        rng = np.random.default_rng(seed)
        od = make_order_draws(g_cls, g_ord, g_sp, K, n_draws, collapse, rng, min_base, min_new, g_cl)
        for key, (labels, g_r, cls) in fam_inputs.items():
            for lab in labels:
                res = eval_order_family(cls[lab], g_r, g_sp, g_cls, g_ord, od, collapse)
                for o, group, k, _ in od:
                    v = res[(o, group)]
                    if o <= -2:
                        ci = -o - 2
                        nb, oname, ocl = 0, "CLASS:" + str(cl_names[ci]), str(cl_names[ci])
                        nn = int(sum(1 for i in np.where(g_cls == 1)[0] if g_ord[i] not in base_ord and g_cl[i] == ci))
                    elif o < 0:
                        nb, nn, oname, ocl = 0, n_pool, "NEW_ONLY_POOL", "NA"
                    else:
                        nb = int(((g_ord == o) & (g_cls == 0)).sum())
                        nn = int(((g_ord == o) & (g_cls == 1)).sum())
                        oname = ord_names[o]
                        ocl = str(cl_names[ord_cl[o]])
                    out.append([key, lab, variant, oname, group, nb, nn, k, n_draws,
                                float(np.nanmean(v)), fn.q(v, 2.5), fn.q(v, 97.5), ocl])
    return out


SUMMARY_SQL = """
with t as (select * from read_csv('{f}', delim='\t', header=true, all_varchar=true)),
x as (select family, tier, variant, "order", "group", cast(n_base as int) nb, cast(n_new as int) nn,
             cast(k as int) k, try_cast(stat_mean as double) s from t where stat_mean <> 'NA' and k = '{K}'),
b as (select family, tier, variant, count(*) n_base_orders, quantile_cont(s, 0.5) base_out_median,
             quantile_cont(s, 0.1) base_out_p10, quantile_cont(s, 0.9) base_out_p90
      from x where "group" = 'base_out' group by all),
nr as (select family, tier, variant, count(*) n_new_rich, quantile_cont(s, 0.5) new_in_base_orders_median
       from x where "group" = 'new' and nb >= {mb} group by all),
n0 as (select family, tier, variant, count(*) n_new_only, quantile_cont(s, 0.5) new_only_orders_median
       from x where "group" = 'new_only_pool' group by all)
select b.*, nr.n_new_rich, nr.new_in_base_orders_median, n0.n_new_only, n0.new_only_orders_median,
       n0.new_only_orders_median - b.base_out_median as new_only_minus_base_out
from b left join nr using (family, tier, variant) left join n0 using (family, tier, variant)
order by family, tier, variant
"""


def aggregate(outdir, K, min_base):
    files = sorted(f for f in os.listdir(outdir) if f.endswith(".ordnull.tsv"))
    allf = os.path.join(outdir, "ordnull_all.tsv")
    with open(allf, "w") as out:
        out.write("\t".join(COLS) + "\n")
        for f in files:
            with open(os.path.join(outdir, f)) as fh:
                next(fh)
                out.write(fh.read())
    con = duckdb.connect()
    con.execute(f"copy ({SUMMARY_SQL.format(f=allf, K=K, mb=min_base)}) to "
                f"'{os.path.join(outdir, 'order_null_summary.tsv')}' (delimiter '\t', header)")
    con.execute(f"""copy (
      with b as (select family, tier, variant, cls, count(*) n_base_orders, median(try_cast(stat_mean as double)) base_med
                 from read_csv('{allf}', delim='\t', header=true, all_varchar=true)
                 where "group" = 'base_out' and k = '{K}' and stat_mean <> 'NA' group by all),
           n as (select family, tier, variant, cls, try_cast(stat_mean as double) new_stat, n_new
                 from read_csv('{allf}', delim='\t', header=true, all_varchar=true)
                 where "group" = 'new_only_class' and stat_mean <> 'NA')
      select n.family, n.tier, n.variant, n.cls, n.n_new, n.new_stat, b.base_med, b.n_base_orders, n.new_stat - b.base_med diff
      from n join b using (family, tier, variant, cls) order by 1, 2, 3, 4)
      to '{os.path.join(outdir, 'order_null_class.tsv')}' (delimiter '\t', header)""")
    con.close()
    print("aggregated", len(files))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--members-dir", required=True)
    ap.add_argument("--clusters-dir", required=True)
    ap.add_argument("--tables-dir", required=True)
    ap.add_argument("--classes", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--families", nargs="*", default=None)
    ap.add_argument("-K", type=int, default=10, help="genomes per draw")
    ap.add_argument("--min-base", type=int, default=20)
    ap.add_argument("--min-new", type=int, default=10)
    ap.add_argument("--draws", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--aggregate", action="store_true")
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    if a.aggregate:
        aggregate(a.outdir, a.K, a.min_base)
        return 0
    g_rows = fn.load_genomes(a)
    asm, *_ = fn.genome_arrays(g_rows)
    keys = a.families or sorted(f[:-13] for f in os.listdir(a.clusters_dir) if f.endswith(".clusters.tsv"))
    for key in keys:
        labels, g_r, cls = fn.family_arrays(key, a, asm)
        rows = run_order_null(g_rows, {key: (labels, g_r, cls)}, a.K, a.draws, a.seed, a.min_base, a.min_new)
        with open(os.path.join(a.outdir, key + ".ordnull.tsv"), "w") as fh:
            fh.write("\t".join(COLS) + "\n")
            for r in rows:
                fh.write("\t".join(fn.fmt(v) for v in r) + "\n")
        print("wrote", key, len(rows), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
