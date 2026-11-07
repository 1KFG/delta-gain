#!/usr/bin/env python3
"""F4 artefact / quality diagnostics for within-family novelty (docs/FAMILY_DIVERSITY_PLAN.md, stage F4).

For each family and tier t, a NEW region is "novel" if its tier-t cluster has no baseline member (same
definition as F2/F3). We ask whether novelty is concentrated where artefacts live:
  * low genome quality (BUSCO complete %, assembly N50)        -> qc_by_quality.tsv
  * a single species/genome, or clusters private to one species -> family_qc.tsv
  * truncated / partial regions (start == 1 or end == protein_len), low model coverage -> family_qc.tsv
  * taxonomy (phylum / order of the novel regions)             -> family_qc.tsv columns top_phylum*
Per-family metrics are then averaged within a group (enriched Pfam, control Pfam, CAZy, MEROPS).
"""
import argparse
import os
import sys

import duckdb

SQL_LOAD = """
create or replace temp table m as
select regexp_replace(regexp_extract(filename, '([^/]+)\\.members\\.tsv$', 1), '^\\s+|\\s+$', '') as family,
       region_id, asmid, replace("class", chr(13), '') cls, species,
       try_cast(start as int) s, try_cast("end" as int) e, try_cast(protein_len as int) plen,
       try_cast(model_cov as double) mcov, try_cast(region_len as int) rlen
from read_csv({mfiles}, delim='\t', header=true, all_varchar=true, filename=true, union_by_name=true)
"""


def run(a):
    con = duckdb.connect()
    con.execute(f"PRAGMA threads={a.threads}")
    fams = [l.strip() for l in open(a.groups).read().split("\n") if l.strip()]
    grp = {}
    for l in fams:
        f, g = l.split("\t")
        grp[f] = g
    mfiles = [os.path.join(a.members_dir, f + ".members.tsv") for f in grp
              if os.path.exists(os.path.join(a.members_dir, f + ".members.tsv"))
              and os.path.exists(os.path.join(a.clusters_dir, f + ".clusters.tsv"))]
    cfiles = [os.path.join(a.clusters_dir, os.path.basename(m).replace(".members.tsv", ".clusters.tsv")) for m in mfiles]
    con.execute(SQL_LOAD.format(mfiles=repr(mfiles)))
    con.execute("create or replace temp table g(family varchar, grp varchar)")
    con.executemany("insert into g values (?, ?)", list(grp.items()))
    tier_cols = a.tiers
    sel = ", ".join(tier_cols)
    con.execute(f"""
      create or replace temp table c as
      select regexp_extract(filename, '([^/]+)\\.clusters\\.tsv$', 1) family, region_id, {sel}
      from read_csv({cfiles!r}, delim='\t', header=true, all_varchar=true, filename=true, union_by_name=true)""")
    con.execute(f"""
      create or replace temp table q as
      select s.ASMID asmid, any_value(b.complete_pct) busco, any_value(x.N50_bp) n50, any_value(x.contig_count) contigs,
             any_value(s.PHYLUM) phylum, any_value(s."ORDER") ord
      from read_parquet('{a.tables_dir}/species.parquet') s
      left join read_parquet('{a.tables_dir}/busco_genome.parquet') b on b.ASMID = s.ASMID
      left join read_parquet('{a.tables_dir}/asm_stats.parquet') x on x.ASMID = s.ASMID
      group by s.ASMID""")
    con.execute("""
      create or replace temp table j as
      select m.*, c.* exclude (family, region_id), q.busco, q.n50, q.contigs, q.phylum, q.ord,
             (s = 1 or e >= plen) trunc
      from m join c using (family, region_id) left join q using (asmid)""")
    fam_rows, bin_rows, ord_rows = [], [], []
    for t in tier_cols:
        con.execute(f"""
          create or replace temp table cl as
          select family, {t} cid, max((cls = 'baseline')::int) hb, count(distinct asmid) ng, count(distinct species) nsp
          from j group by 1, 2""")
        con.execute(f"""
          create or replace temp table n as
          select j.*, (cl.hb = 0) nov, cl.nsp cl_nsp, cl.ng cl_ng
          from j join cl on cl.family = j.family and cl.cid = j.{t}
          where j.cls = 'new'""")
        fam_rows += con.execute(f"""
          with sp as (select family, species, count(*) k from n where nov group by 1, 2),
               top as (select family, max(k) kmax, sum(k) ktot from sp group by 1),
               ph as (select family, phylum, count(*) k, row_number() over (partition by family order by count(*) desc) r
                      from n where nov group by 1, 2)
          select '{t}' tier, n.family, any_value(g.grp) grp, count(*) n_new, avg(nov::int) novel_rate,
            sum(nov::int) n_novel, count(distinct case when nov then species end) novel_species,
            any_value(top.kmax) / nullif(any_value(top.ktot), 0) top_species_share,
            avg(case when nov and cl_nsp = 1 then 1.0 when nov then 0.0 end) frac_novel_private_species,
            avg(case when nov then trunc::int end) trunc_novel, avg(case when not nov then trunc::int end) trunc_known,
            avg(case when nov then mcov end) mcov_novel, avg(case when not nov then mcov end) mcov_known,
            avg(case when nov then rlen end) rlen_novel, avg(case when not nov then rlen end) rlen_known,
            any_value(ph.phylum) top_phylum_novel
          from n join g using (family) left join top using (family) left join ph on ph.family = n.family and ph.r = 1
          group by n.family""").fetchall()
        # within-order contrast: novelty in low-quality vs high-quality genomes of the SAME order (clade-matched)
        ord_rows += con.execute(f"""
          with b as (select family, ord, nov, busco < 95 lowb, n50 < 50000 lown from n where busco is not null and n50 is not null),
          r as (select family, ord, 'busco' metric, lowb low, count(*) k, avg(nov::int) rate from b group by all
                union all select family, ord, 'n50', lown, count(*), avg(nov::int) from b group by all),
          w as (select family, ord, metric,
                  max(case when low then rate end) rl, max(case when not low then rate end) rh,
                  max(case when low then k end) kl, max(case when not low then k end) kh from r group by all)
          select '{t}' tier, g.grp, w.metric, w.family, count(*) n_orders, avg(rl - rh) mean_diff_low_minus_high
          from w join g using (family) where kl >= {a.min_bin} and kh >= {a.min_bin} group by all""").fetchall()
        bin_rows += con.execute(f"""
          with b as (select family, nov,
              case when busco is null then 'NA' when busco < 80 then 'busco<80' when busco < 95 then 'busco80-95' else 'busco>=95' end qb,
              case when n50 is null then 'NA' when n50 < 50000 then 'n50<50k' when n50 < 500000 then 'n50 50k-500k' else 'n50>=500k' end nb
              from n)
          select '{t}' tier, g.grp, 'busco' metric, qb bin, family, count(*) n, avg(nov::int) novel_rate from b join g using (family) group by all
          union all
          select '{t}', g.grp, 'n50', nb, family, count(*), avg(nov::int) from b join g using (family) group by all""").fetchall()
    con.close()
    return fam_rows, bin_rows, ord_rows


FAM_COLS = ["tier", "family", "grp", "n_new", "novel_rate", "n_novel", "novel_species", "top_species_share",
            "frac_novel_private_species", "trunc_novel", "trunc_known", "mcov_novel", "mcov_known",
            "rlen_novel", "rlen_known", "top_phylum_novel"]


def fmt(v):
    return "NA" if v is None else (f"{v:.5g}" if isinstance(v, float) else str(v))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--members-dir", required=True)
    ap.add_argument("--clusters-dir", required=True)
    ap.add_argument("--tables-dir", required=True)
    ap.add_argument("--groups", required=True, help="TSV: family_key<TAB>group")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tiers", nargs="+", default=["c70", "c50"])
    ap.add_argument("--min-bin", type=int, default=50, help="min regions of a family in a quality bin")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    fam_rows, bin_rows, ord_rows = run(a)
    with open(os.path.join(a.outdir, "family_qc.tsv"), "w") as fh:
        fh.write("\t".join(FAM_COLS) + "\n")
        for r in sorted(fam_rows):
            fh.write("\t".join(fmt(v) for v in r) + "\n")
    con = duckdb.connect()
    con.execute("create table b(tier varchar, grp varchar, metric varchar, bin varchar, family varchar, n bigint, novel_rate double)")
    con.executemany("insert into b values (?,?,?,?,?,?,?)", bin_rows)
    con.execute(f"""copy (select tier, grp, metric, bin, count(*) families, round(avg(novel_rate), 4) mean_novel_rate,
                      round(median(novel_rate), 4) median_novel_rate, sum(n) regions
                    from b where n >= {a.min_bin} group by all order by tier, grp, metric, bin)
                    to '{os.path.join(a.outdir, "qc_by_quality.tsv")}' (delimiter '\t', header)""")
    con.execute("create table o(tier varchar, grp varchar, metric varchar, family varchar, n_orders bigint, d double)")
    con.executemany("insert into o values (?,?,?,?,?,?)", ord_rows)
    con.execute(f"""copy (select tier, grp, metric, count(*) families, round(avg(d), 4) mean_diff_low_minus_high,
                      round(median(d), 4) median_diff, round(avg((d > 0)::int), 3) frac_families_low_higher,
                      round(avg(n_orders), 1) mean_orders_compared
                    from o group by all order by tier, grp, metric)
                    to '{os.path.join(a.outdir, "qc_within_order.tsv")}' (delimiter '\t', header)""")
    con.close()
    print("families", len({r[1] for r in fam_rows}), "rows", len(fam_rows))
    return 0


if __name__ == "__main__":
    sys.exit(main())
