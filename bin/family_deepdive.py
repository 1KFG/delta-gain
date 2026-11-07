#!/usr/bin/env python3
"""F5 deep dive for selected families (S08, GH18, AA9, M14 ...): novelty by sub-family label (`sub` column of the members table),
by order, and taxonomic spread of novel clusters. Tier t: a region is novel if its cluster has no baseline member."""
import argparse, os
import duckdb

ap = argparse.ArgumentParser()
ap.add_argument("--members-dir", required=True); ap.add_argument("--clusters-dir", required=True)
ap.add_argument("--tables-dir", required=True); ap.add_argument("--families", nargs="+", required=True)
ap.add_argument("--tier", default="c70"); ap.add_argument("--outdir", required=True)
ap.add_argument("--min-sub", type=int, default=200)
a = ap.parse_args()
c = duckdb.connect(); c.execute("pragma threads=8")
c.execute(f"""create table sp as select ASMID asmid, any_value("ORDER") ord, any_value(CLASS) cl from read_parquet('{a.tables_dir}/species.parquet') group by 1""")
out = []
for fam in a.families:
    mt, ct = f"{a.members_dir}/{fam}.members.tsv", f"{a.clusters_dir}/{fam}.clusters.tsv"
    c.execute(f"""create or replace table j as select m.region_id, m.asmid, replace(m."class", chr(13), '') cls, m.species, m."sub" sb,
        try_cast(m.model_cov as double) mcov, try_cast(m.region_len as int) rlen, k.{a.tier} cid, sp.ord, sp.cl
      from read_csv('{mt}', delim='\t', header=true, all_varchar=true) m join read_csv('{ct}', delim='\t', header=true, all_varchar=true) k using (region_id)
      left join sp on sp.asmid = m.asmid""")
    c.execute("create or replace table cl as select cid, max((cls='baseline')::int) hb, count(distinct species) nsp from j group by 1")
    c.execute("create or replace table n as select j.*, (cl.hb = 0) nov, cl.nsp from j join cl using (cid) where cls = 'new'")
    print("==", fam, a.tier)
    print(c.execute(f"""select sb sub, count(*) new_regions, round(avg(nov::int),3) novel_rate, count(distinct case when nov then cid end) novel_clusters,
        count(distinct case when nov then species end) novel_species, round(avg(mcov),2) mcov, round(median(rlen)) med_len
        from n group by 1 having count(*) >= {a.min_sub} order by 3 desc limit 15""").fetchdf().to_string(index=False))
    print(c.execute("""select cl, ord, count(*) new_regions, round(avg(nov::int),3) novel_rate, sum(nov::int) novel_regions from n group by 1,2 having count(*) >= 300 order by novel_regions desc limit 8""").fetchdf().to_string(index=False))
    print(c.execute("""select count(distinct cid) clusters, count(distinct case when hb=0 then cid end) baseline_free_clusters,
        round(avg(case when hb=0 then nsp end),2) mean_species_in_baseline_free from cl""").fetchdf().to_string(index=False))
