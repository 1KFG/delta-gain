#!/usr/bin/env python3
"""E7 stage 3: SSP / signal-peptide fractions and CAZy / MEROPS copy number, new vs baseline within order (species level).

Units are (species, class): all genomes of a species in one class are averaged, so strains do not inflate n. Within-order
contrast = sum_o w_o (mean_new_o - mean_base_o) / sum_o w_o, w_o = nb*nn/(nb+nn) over orders with >= --min-units units of each class.
Significance: SSP/SP by permuting class labels of units inside each order; family copy numbers by bootstrapping units inside
(order, class) (CI, two-sided p, BH q). 'N_known' = new-class units whose species also has a baseline genome, 'N_newsp' = the rest.
Survivors: observed vs genome-rate expected signal-peptide / SSP counts for survivor representatives, by genome-count class.
"""
import argparse
import os
import sys
from collections import defaultdict

import duckdb
import numpy as np


def wdiff(vals, lab, ords, groups, a_lab, b_lab):
    """weighted within-order mean(a) - mean(b); vals: (n,) or (n, f); lab: int labels; ords: int order index per unit."""
    num, den = 0.0, 0.0
    for idx in groups:
        la = lab[idx]
        ia, ib = idx[la == a_lab], idx[la == b_lab]
        if len(ia) < 1 or len(ib) < 1:
            continue
        w = len(ia) * len(ib) / (len(ia) + len(ib))
        num = num + w * (vals[ia].mean(0) - vals[ib].mean(0))
        den += w
    return num / den if den else np.nan


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, help="results/e7 (stage 1/2 outputs)")
    ap.add_argument("--survivors", required=True)
    ap.add_argument("--perms", type=int, default=5000)
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--min-units", type=int, default=3)
    ap.add_argument("--min-genomes-fam", type=int, default=50)
    ap.add_argument("--seed", type=int, default=3)
    a = ap.parse_args(argv)
    d = a.dir
    c = duckdb.connect()
    c.execute("pragma threads=16")
    rng = np.random.default_rng(a.seed)
    # per-genome SSP counts
    c.execute(f"""create table ssp as
      select split_part(cand.protein_id, '_', 1) prefix, count(*) n_ssp
      from read_parquet('{d}/cand.parquet') cand join read_parquet('{d}/cys.parquet') y using (protein_id) where y.cys >= 4 group by 1""")
    c.execute(f"""create table g as
      select m.asmid, m.prefix, m.cl, m.species, m.ord, n.n_prot, n.n_sp, n.n_cand, coalesce(s.n_ssp, 0) n_ssp
      from read_parquet('{d}/genome_meta.parquet') m join read_parquet('{d}/genome_counts.parquet') n using (prefix)
      left join ssp s using (prefix) where n.n_prot > 0""")
    c.execute(f"copy g to '{d}/e7_genome_counts.tsv' (delimiter '\t', header)")
    out = open(os.path.join(d, "e7_summary.txt"), "w")

    def P(*x):
        s = " ".join(str(v) for v in x)
        print(s, flush=True)
        out.write(s + "\n")

    tot = c.execute("select cl, count(*), sum(n_prot), sum(n_sp), sum(n_cand), sum(n_ssp) from g group by 1 order by 1").fetchall()
    for cl, n, npr, nsp, ncd, nss in tot:
        P(f"genomes[{cl}]={n} proteins={npr} signal_peptide={nsp} ({nsp / npr:.4f}) candidates_pre_cys={ncd} SSP={nss} ({nss / npr:.5f})")
    # ---------------- species-level units
    c.execute("""create table u as
      select species, cl, any_value(ord) ord, count(*) ng, avg(n_sp::double / n_prot) sp_frac, avg(n_ssp::double / n_prot) ssp_frac,
             avg(n_ssp::double / nullif(n_sp, 0)) ssp_of_sp
      from g where ord is not null group by species, cl""")
    c.execute("""create table bsp as select distinct species from u where cl = 'baseline'""")
    rows = c.execute("""select u.species, u.cl, u.ord, u.ng, u.sp_frac, u.ssp_frac, coalesce(u.ssp_of_sp, 0),
        case when u.cl = 'baseline' then 'B' when u.species in (select species from bsp) then 'N_known' else 'N_newsp' end kind
        from u""").fetchall()
    kinds = np.array([r[7] for r in rows])
    ords = np.array([r[2] for r in rows])
    ord_ids = {o: i for i, o in enumerate(sorted(set(ords)))}
    oi = np.array([ord_ids[o] for o in ords])
    metrics = {"sp_frac": 4, "ssp_frac": 5, "ssp_of_sp": 6}
    V = {k: np.array([r[i] for r in rows], dtype=float) for k, i in metrics.items()}
    P("units", {k: int((kinds == k).sum()) for k in ("B", "N_known", "N_newsp")})
    with open(os.path.join(d, "e7_ssp_within_order.tsv"), "w") as fo:
        fo.write("metric\tcontrast\torders\tunits_a\tunits_b\tmean_a\tmean_b\twithin_order_diff\tci_lo\tci_hi\tperm_p\n")
        for contrast, ka in (("new_vs_base", ("N_known", "N_newsp")), ("N_known_vs_base", ("N_known",)), ("N_newsp_vs_base", ("N_newsp",))):
            sel = np.isin(kinds, list(ka) + ["B"])
            lab = np.where(np.isin(kinds, list(ka)), 1, np.where(kinds == "B", 0, -1))
            groups = []
            for o in range(len(ord_ids)):
                idx = np.flatnonzero((oi == o) & sel)
                if (lab[idx] == 1).sum() >= a.min_units and (lab[idx] == 0).sum() >= a.min_units:
                    groups.append(idx)
            for mname, vals in V.items():
                obs = wdiff(vals, lab, oi, groups, 1, 0)
                # permutation
                ge = 0
                for _ in range(a.perms):
                    pl = lab.copy()
                    for idx in groups:
                        pl[idx] = rng.permutation(lab[idx])
                    if abs(wdiff(vals, pl, oi, groups, 1, 0)) >= abs(obs):
                        ge += 1
                # bootstrap CI over units inside (order, class)
                bs = []
                for _ in range(a.boot):
                    num = den = 0.0
                    for idx in groups:
                        ia, ib = idx[lab[idx] == 1], idx[lab[idx] == 0]
                        ra, rb = rng.choice(ia, len(ia)), rng.choice(ib, len(ib))
                        w = len(ia) * len(ib) / (len(ia) + len(ib))
                        num += w * (vals[ra].mean() - vals[rb].mean())
                        den += w
                    bs.append(num / den)
                lo, hi = np.percentile(bs, [2.5, 97.5])
                ina = np.concatenate([g[lab[g] == 1] for g in groups])
                inb = np.concatenate([g[lab[g] == 0] for g in groups])
                fo.write(f"{mname}\t{contrast}\t{len(groups)}\t{len(ina)}\t{len(inb)}\t{vals[ina].mean():.5f}\t{vals[inb].mean():.5f}\t"
                         f"{obs:.5f}\t{lo:.5f}\t{hi:.5f}\t{(ge + 1) / (a.perms + 1):.4f}\n")
    P(open(os.path.join(d, "e7_ssp_within_order.tsv")).read())
    # ---------------- survivors
    c.execute(f"""create table sv as
      select split_part(rep, '__', 2) asmid, split_part(rep, '__', 3) pid, n_genomes,
        case when n_genomes = 1 then '1' when n_genomes < 5 then '2-4' when n_genomes < 20 then '5-19' else '>=20' end cls
      from read_csv('{a.survivors}', delim='\t', header=true)""")
    c.execute(f"""create table svf as
      select sv.cls, sv.asmid, (f.protein_id is not null)::int is_sp, coalesce(f.cand, 0) cand, (y.cys >= 4)::int ssp4
      from sv left join read_parquet('{d}/sp_flag.parquet') f on f.protein_id = sv.pid
      left join read_parquet('{d}/cys.parquet') y on y.protein_id = sv.pid""")
    c.execute("""create table gr as select asmid, n_sp::double / n_prot sp_rate, n_ssp::double / n_prot ssp_rate from g""")
    res = c.execute("""
      select svf.cls, count(*) n, sum(is_sp) sp, sum(coalesce(ssp4, 0) * cand) ssp, sum(gr.sp_rate) exp_sp, sum(gr.ssp_rate) exp_ssp
      from svf join gr using (asmid) group by 1 order by 1""").fetchall()
    newt = c.execute("select sum(n_prot), sum(n_sp), sum(n_ssp) from g where cl = 'new'").fetchone()
    P("survivor representatives vs genome-rate expectation (same genome's overall rate), and vs all new proteins")
    P(f"all new-genome proteins: n={newt[0]} SP frac={newt[1] / newt[0]:.4f} SSP frac={newt[2] / newt[0]:.5f}")
    with open(os.path.join(d, "e7_survivor_ssp.tsv"), "w") as fo:
        fo.write("genome_class\treps\tsp\tsp_frac\texp_sp_frac\tsp_obs_over_exp\tssp\tssp_frac\texp_ssp_frac\tssp_obs_over_exp\n")
        for cls, n, sp, ssp, esp, essp in res:
            fo.write(f"{cls}\t{n}\t{sp}\t{sp / n:.4f}\t{esp / n:.4f}\t{sp / esp:.3f}\t{ssp}\t{ssp / n:.5f}\t{essp / n:.5f}\t{ssp / essp:.3f}\n")
    P(open(os.path.join(d, "e7_survivor_ssp.tsv")).read())
    # ---------------- family copy number
    fam = c.execute(f"""select f.prefix, f.kind || ':' || f.fam fam, f.n from read_parquet('{d}/fam_counts.parquet') f""").fetchall()
    c.execute(f"""create table gg as select prefix, asmid, species, cl, ord, n_prot from g where ord is not null""")
    grows = c.execute("select prefix, species, cl, ord, n_prot from gg").fetchall()
    gidx = {r[0]: i for i, r in enumerate(grows)}
    fams = defaultdict(int)
    for p, fm, n in fam:
        if p in gidx:
            fams[fm] += 1
    famlist = sorted(f for f, k in fams.items() if k >= a.min_genomes_fam)
    fi = {f: i for i, f in enumerate(famlist)}
    Mraw = np.zeros((len(grows), len(famlist)), dtype=np.float32)
    for p, fm, n in fam:
        if p in gidx and fm in fi:
            Mraw[gidx[p], fi[fm]] = n
    nprot = np.array([r[4] for r in grows], dtype=np.float32)
    Mnorm = Mraw / nprot[:, None] * 1000.0
    ukey = {}
    uidx = np.zeros(len(grows), dtype=int)
    for i, r in enumerate(grows):
        k = (r[1], r[2])
        if k not in ukey:
            ukey[k] = len(ukey)
        uidx[i] = ukey[k]
    nU = len(ukey)

    def unitmean(M):
        out_ = np.zeros((nU, M.shape[1]), dtype=np.float64)
        np.add.at(out_, uidx, M)
        cnt = np.bincount(uidx, minlength=nU)[:, None]
        return out_ / cnt

    ulist = sorted(ukey, key=lambda k: ukey[k])
    ucl = np.array([k[1] for k in ulist])
    first = {}
    for i, r in enumerate(grows):
        first.setdefault(uidx[i], r[3])
    uord = np.array([first[i] for i in range(nU)])
    ulab = np.where(ucl == "new", 1, 0)
    uoi = np.array([ord_ids.get(o, -1) for o in uord])
    # orders in the copy-number analysis: >= min-units units of each class
    cgroups = []
    for o in sorted(set(uoi)):
        if o < 0:
            continue
        idx = np.flatnonzero(uoi == o)
        if (ulab[idx] == 1).sum() >= a.min_units and (ulab[idx] == 0).sum() >= a.min_units:
            cgroups.append(idx)
    P("copy-number orders", len(cgroups), "families", len(famlist), "units", nU)
    res_rows = []
    for tag, M in (("raw", Mraw), ("per1000", Mnorm)):
        Um = unitmean(M)
        obs = wdiff(Um, ulab, uoi, cgroups, 1, 0)
        # weighted base / new means for relative effect
        def wmeans(U, groups, lab_):
            nb, nn, den = 0.0, 0.0, 0.0
            for idx in groups:
                ia, ib = idx[lab_[idx] == 1], idx[lab_[idx] == 0]
                w = len(ia) * len(ib) / (len(ia) + len(ib))
                nb = nb + w * U[ib].mean(0)
                nn = nn + w * U[ia].mean(0)
                den += w
            return nb / den, nn / den
        mb, mn = wmeans(Um, cgroups, ulab)
        bs = np.zeros((a.boot, len(famlist)))
        for b in range(a.boot):
            num = 0.0
            den = 0.0
            for idx in cgroups:
                ia, ib = idx[ulab[idx] == 1], idx[ulab[idx] == 0]
                ra, rb = rng.choice(ia, len(ia)), rng.choice(ib, len(ib))
                w = len(ia) * len(ib) / (len(ia) + len(ib))
                num = num + w * (Um[ra].mean(0) - Um[rb].mean(0))
                den += w
            bs[b] = num / den
        lo, hi = np.percentile(bs, [2.5, 97.5], axis=0)
        p = 2 * np.minimum((bs <= 0).mean(0), (bs >= 0).mean(0))
        p = np.maximum(p, 1.0 / (a.boot + 1))
        # BH
        order = np.argsort(p)
        q = np.empty_like(p)
        m_ = len(p)
        prev = 1.0
        for rank, i in zip(range(m_, 0, -1), order[::-1]):
            prev = min(prev, p[i] * m_ / rank)
            q[i] = prev
        for j, f in enumerate(famlist):
            res_rows.append((tag, f, fams[f], mb[j], mn[j], obs[j], lo[j], hi[j], p[j], q[j]))
    with open(os.path.join(d, "e7_family_copy_number.tsv"), "w") as fo:
        fo.write("scale\tfamily\tgenomes_with_hits\tbase_mean\tnew_mean\twithin_order_diff\tci_lo\tci_hi\tp\tq_bh\tdirection\n")
        for r in res_rows:
            direction = "expanded" if r[9] < 0.05 and r[5] > 0 else ("contracted" if r[9] < 0.05 and r[5] < 0 else "ns")
            fo.write("\t".join(str(x) if isinstance(x, str) or isinstance(x, int) else f"{x:.5g}" for x in r) + f"\t{direction}\n")
    for tag in ("raw", "per1000"):
        sub = [r for r in res_rows if r[0] == tag]
        exp_ = [r for r in sub if r[9] < 0.05 and r[5] > 0]
        con_ = [r for r in sub if r[9] < 0.05 and r[5] < 0]
        P(f"copy number [{tag}]: families tested={len(sub)} expanded(q<0.05)={len(exp_)} contracted={len(con_)}")
        for kind in ("cazy", "merops"):
            s2 = [r for r in sub if r[1].startswith(kind)]
            P(f"   {kind}: tested={len(s2)} expanded={sum(1 for r in s2 if r[9] < 0.05 and r[5] > 0)} contracted={sum(1 for r in s2 if r[9] < 0.05 and r[5] < 0)}")
        for r in sorted(exp_, key=lambda r: -(r[5] / max(r[3], 1e-9)))[:12]:
            P(f"   top expanded {r[1]}: base={r[3]:.3f} new={r[4]:.3f} diff={r[5]:.3f} [{r[6]:.3f},{r[7]:.3f}] q={r[9]:.3g}")
        for r in sorted(con_, key=lambda r: r[5] / max(r[3], 1e-9))[:8]:
            P(f"   top contracted {r[1]}: base={r[3]:.3f} new={r[4]:.3f} diff={r[5]:.3f} [{r[6]:.3f},{r[7]:.3f}] q={r[9]:.3g}")
    out.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
