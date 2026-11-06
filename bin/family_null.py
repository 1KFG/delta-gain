#!/usr/bin/env python3
"""F3 baseline-vs-baseline null for within-family novelty (see docs/FAMILY_DIVERSITY_PLAN.md, stage F3).

Question: do NEW genomes add more family diversity than the same number of extra BASELINE genomes
from the same taxa would?

Each draw (order-matched, paired):
  H = baseline genomes held out, N = equally many NEW genomes, drawn per taxonomic order
      (m_o = max(1, round(f * min(n_base_o, n_new_o))) over orders present in both sets).
  Reference R = all baseline genomes not in H.
  stat_new  = fraction of N's regions (or distinct species-cluster pairs) whose tier-t cluster has no member in R
  stat_null = the same for H (so H and N are measured against the SAME reference)
  diff      = stat_new - stat_null  (paired per draw)
Variants: `genome` (regions of all drawn genomes) and `species` (one random genome per species per pool and
distinct (species, cluster) pairs: heavily sequenced species cannot dominate).
Each also reported for the "known-species" stratum: only regions/pairs whose species has a genome in R.
Draws are generated once per (fraction, variant) and reused for every family.
"""
import argparse
import csv
import os
import sys

import duckdb
import numpy as np

COLS = ["family", "tier", "variant", "frac", "draws", "genomes_H", "genomes_N", "orders",
        "stat_new", "stat_null", "null_lo", "null_hi", "diff_mean", "diff_lo", "diff_hi", "p_diff_le0",
        "known_new", "known_null", "known_diff_mean", "known_diff_lo", "known_diff_hi", "known_p_diff_le0"]


def one_per_species(rng, idx, g_sp):
    sh = rng.permutation(idx)
    _, first = np.unique(g_sp[sh], return_index=True)
    return sh[first]


def make_draws(g_cls, g_ord, g_sp, f, collapse, n_draws, rng):
    pools = []
    for o in np.unique(g_ord):
        b = np.where((g_ord == o) & (g_cls == 0))[0]
        n = np.where((g_ord == o) & (g_cls == 1))[0]
        if len(b) and len(n):
            pools.append((b, n))
    draws = []
    for _ in range(n_draws):
        hs, ns = [], []
        for b, n in pools:
            if collapse:
                b, n = one_per_species(rng, b, g_sp), one_per_species(rng, n, g_sp)
            m = min(max(1, int(round(f * min(len(b), len(n))))), len(b), len(n))
            hs.append(rng.choice(b, m, replace=False))
            ns.append(rng.choice(n, m, replace=False))
        draws.append((np.concatenate(hs), np.concatenate(ns)))
    return draws, len(pools)


def eval_family(cl, g_r, sp_g, g_cls, draws, collapse, n_sp):
    """cl: cluster idx per region; g_r: genome idx per region. Returns arrays over draws."""
    ncl = int(cl.max()) + 1
    sp_r = sp_g[g_r]
    base_r = g_cls[g_r] == 0
    ng = len(g_cls)
    out = {k: [] for k in ("new", "null", "knew", "knull")}
    for H, N in draws:
        inH = np.zeros(ng, bool)
        inH[H] = True
        inN = np.zeros(ng, bool)
        inN[N] = True
        ref = base_r & ~inH[g_r]
        present = np.zeros(ncl, bool)
        present[cl[ref]] = True
        spref = np.zeros(n_sp, bool)
        refg = (g_cls == 0) & ~inH
        spref[sp_g[refg]] = True
        res = []
        for m in (inN[g_r], inH[g_r]):
            if collapse:
                p = np.unique(sp_r[m].astype(np.int64) * ncl + cl[m])
                c, s = p % ncl, p // ncl
            else:
                c, s = cl[m], sp_r[m]
            nov = ~present[c]
            kn = spref[s]
            res.append((nov.mean() if len(nov) else np.nan,
                        nov[kn].mean() if kn.any() else np.nan))
        (a, ka), (b, kb) = res
        out["new"].append(a); out["knew"].append(ka); out["null"].append(b); out["knull"].append(kb)
    return {k: np.array(v) for k, v in out.items()}


def q(a, p):
    a = a[~np.isnan(a)]
    return float(np.percentile(a, p)) if len(a) else None


def summarise(res):
    d = res["new"] - res["null"]
    kd = res["knew"] - res["knull"]
    ok, kok = ~np.isnan(d), ~np.isnan(kd)
    return dict(
        stat_new=float(np.nanmean(res["new"])), stat_null=float(np.nanmean(res["null"])),
        null_lo=q(res["null"], 2.5), null_hi=q(res["null"], 97.5),
        diff_mean=float(np.nanmean(d)), diff_lo=q(d, 2.5), diff_hi=q(d, 97.5),
        p=float((d[ok] <= 0).mean()) if ok.any() else None,
        known_new=float(np.nanmean(res["knew"])) if kok.any() else None,
        known_null=float(np.nanmean(res["knull"])) if kok.any() else None,
        kdm=float(np.nanmean(kd)) if kok.any() else None, kdl=q(kd, 2.5), kdh=q(kd, 97.5),
        kp=float((kd[kok] <= 0).mean()) if kok.any() else None)


def load_genomes(a):
    con = duckdb.connect()
    rows = con.execute(f"""
      select s.ASMID, replace(c."class", chr(13), ''), coalesce(nullif(s."ORDER", ''), 'NA'),
             coalesce(nullif(s.SPECIES, ''), s.ASMID)
      from read_parquet('{a.tables_dir}/species.parquet') s
      join read_csv('{a.classes}', delim='\t', header=true, all_varchar=true) c on c.ASMID = s.ASMID
      where replace(c."class", chr(13), '') in ('baseline', 'new')
        and replace(c.pep_status, chr(13), '') = 'matched'""").fetchall()
    con.close()
    return rows


def genome_arrays(rows):
    asm = {r[0]: i for i, r in enumerate(rows)}
    g_cls = np.array([0 if r[1] == "baseline" else 1 for r in rows])
    _, g_ord = np.unique([r[2] for r in rows], return_inverse=True)
    sp_u, g_sp = np.unique([r[3] for r in rows], return_inverse=True)
    return asm, g_cls, g_ord, g_sp, len(sp_u)


def family_arrays(key, a, asm):
    con = duckdb.connect()
    mt = os.path.join(a.members_dir, key + ".members.tsv")
    ct = os.path.join(a.clusters_dir, key + ".clusters.tsv")
    cols = con.execute(f"describe select * from read_csv('{ct}', delim='\t', header=true, all_varchar=true)"
                       ).fetchall()
    labels = [c[0] for c in cols if c[0] != "region_id"]
    sel = ", ".join(f"a.{l}" for l in labels)
    rows = con.execute(f"""select m.asmid, {sel} from
        read_csv('{mt}', delim='\t', header=true, all_varchar=true) m
        join read_csv('{ct}', delim='\t', header=true, all_varchar=true) a using (region_id)""").fetchall()
    con.close()
    g_r = np.array([asm[r[0]] for r in rows])
    cls = {}
    for j, l in enumerate(labels, start=1):
        _, inv = np.unique([r[j] for r in rows], return_inverse=True)
        cls[l] = inv
    return labels, g_r, cls


def run_null(g_rows, fam_inputs, fractions, n_draws, seed, tiers=None):
    """fam_inputs: {key: (labels, g_r, {label: cluster idx})}. Returns list of output rows."""
    asm, g_cls, g_ord, g_sp, n_sp = genome_arrays(g_rows)
    out = []
    for f in fractions:
        for variant in ("genome", "species"):
            collapse = variant == "species"
            rng = np.random.default_rng(seed)
            draws, n_ord = make_draws(g_cls, g_ord, g_sp, f, collapse, n_draws, rng)
            kH = float(np.mean([len(h) for h, _ in draws]))
            kN = float(np.mean([len(n) for _, n in draws]))
            for key, (labels, g_r, cls) in fam_inputs.items():
                for lab in labels:
                    s = summarise(eval_family(cls[lab], g_r, g_sp, g_cls, draws, collapse, n_sp))
                    out.append([key, lab, variant, f, n_draws, round(kH, 1), round(kN, 1), n_ord,
                                s["stat_new"], s["stat_null"], s["null_lo"], s["null_hi"], s["diff_mean"],
                                s["diff_lo"], s["diff_hi"], s["p"], s["known_new"], s["known_null"],
                                s["kdm"], s["kdl"], s["kdh"], s["kp"]])
    return out


def fmt(v):
    return "NA" if v is None or (isinstance(v, float) and np.isnan(v)) else (f"{v:.6g}" if isinstance(v, float) else str(v))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--members-dir", required=True)
    ap.add_argument("--clusters-dir", required=True)
    ap.add_argument("--tables-dir", required=True)
    ap.add_argument("--classes", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--families", nargs="*", default=None)
    ap.add_argument("--fractions", type=float, nargs="+", default=[0.05, 0.2])
    ap.add_argument("--draws", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--aggregate", action="store_true")
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    if a.aggregate:
        files = sorted(f for f in os.listdir(a.outdir) if f.endswith(".null.tsv"))
        with open(os.path.join(a.outdir, "null_all.tsv"), "w") as out:
            out.write("\t".join(COLS) + "\n")
            for f in files:
                with open(os.path.join(a.outdir, f)) as fh:
                    next(fh)
                    out.write(fh.read())
        print("aggregated", len(files))
        return 0
    asm, *_ = genome_arrays(load_genomes(a))
    keys = a.families or sorted(f[:-13] for f in os.listdir(a.clusters_dir) if f.endswith(".clusters.tsv"))
    g_rows = load_genomes(a)
    for key in keys:
        labels, g_r, cls = family_arrays(key, a, asm)
        rows = run_null(g_rows, {key: (labels, g_r, cls)}, a.fractions, a.draws, a.seed)
        with open(os.path.join(a.outdir, key + ".null.tsv"), "w") as fh:
            fh.write("\t".join(COLS) + "\n")
            for r in rows:
                fh.write("\t".join(fmt(v) for v in r) + "\n")
        print("wrote", key, len(rows), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
