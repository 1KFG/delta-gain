#!/usr/bin/env python3
"""E6: taxonomic breadth of funnel survivor clusters (novel at 50 %).

Expands each survivor cluster through the nested tiers (50 -> 70 -> 95) to its NEW proteins and counts distinct species /
genera / orders. Classes: one_species (strain population or orphan), one_genus, one_family, one_order, multi_order.
Outputs survivor_taxa.tsv (per cluster) and survivor_taxa_summary.tsv (by genome-count class, weighted by clusters and proteins)."""
import argparse
import duckdb

ap = argparse.ArgumentParser()
ap.add_argument("--funnel-dir", required=True)
ap.add_argument("--tables-dir", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--threads", type=int, default=16)
ap.add_argument("--memory-limit", default="96GB")
ap.add_argument("--tmp-dir", default=None)
a = ap.parse_args()
F = a.funnel_dir
c = duckdb.connect()
c.execute(f"pragma threads={a.threads}")
c.execute(f"pragma memory_limit='{a.memory_limit}'")
if a.tmp_dir:
    c.execute(f"pragma temp_directory='{a.tmp_dir}'")
rd = lambda p: (f"read_csv('{p}', delim='\t', quote='', escape='', strict_mode=false, header=false, "
                f"columns={{'rep':'VARCHAR','member':'VARCHAR'}})")
c.execute(f"create table s as select rep from read_csv('{F}/summary/survivors_clusters.tsv', delim='\t', header=true)")
c.execute(f"create table l1 as select t.rep srep, t.member m from {rd(F + '/tier_50_clusters.tsv')} t semi join s on s.rep = t.rep")
c.execute(f"create table l2 as select l1.srep, t.member m from {rd(F + '/tier_70_clusters.tsv')} t join l1 on l1.m = t.rep")
c.execute(f"create table l3 as select l2.srep, t.member m from {rd(F + '/tier_95_clusters.tsv')} t join l2 on l2.m = t.rep where t.member like 'NEW__%'")
c.execute(f"""create table sp as select ASMID asmid, any_value(SPECIES) species, any_value(GENUS) genus, any_value(FAMILY) fam, any_value("ORDER") ord
  from read_parquet('{a.tables_dir}/species.parquet') group by 1""")
c.execute("""create table r as select l3.srep rep, count(*) n_proteins, count(distinct split_part(l3.m,'__',2)) n_genomes,
  count(distinct sp.species) n_species, count(distinct sp.genus) n_genera, count(distinct sp.fam) n_families, count(distinct sp.ord) n_orders
  from l3 left join sp on sp.asmid = split_part(l3.m, '__', 2) group by 1""")
c.execute(f"copy r to '{a.outdir}/survivor_taxa.tsv' (delimiter '\t', header)")
c.execute(f"""copy (select case when n_genomes = 1 then '1' when n_genomes < 5 then '2-4' when n_genomes < 20 then '5-19' else '>=20' end genomes_class,
  case when n_species <= 1 then '1 one_species' when n_genera <= 1 then '2 one_genus' when n_families <= 1 then '3 one_family'
       when n_orders <= 1 then '4 one_order' else '5 multi_order' end breadth,
  count(*) clusters, sum(n_proteins) proteins from r group by all order by 1, 2) to '{a.outdir}/survivor_taxa_summary.tsv' (delimiter '\t', header)""")
print(c.execute(f"select * from read_csv('{a.outdir}/survivor_taxa_summary.tsv', delim='\t', header=true)").fetchdf().to_string(index=False))
print(c.execute("select count(*) clusters, sum(n_proteins) proteins from r").fetchall())
