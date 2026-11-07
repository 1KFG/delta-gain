#!/usr/bin/env python3
"""F6: attach Pfam / CAZy / MEROPS labels to funnel survivor clusters (novel at 50 %).

Survivor representatives are NEW__<ASMID>__<BFD protein id>; the BFD id joins directly to the BFD annotation
tables. Labels are for the cluster REPRESENTATIVE only (a cluster of >= 2 proteins is labelled by one member).
Outputs: survivor_annot.tsv (rep, n_new_proteins, n_genomes, pfam, cazy, merops flags + best Pfam acc),
survivor_annot_summary.tsv (by genome-count class), survivor_top_pfam.tsv (most frequent Pfam among survivors seen in >= 5 genomes).
"""
import argparse
import duckdb

ap = argparse.ArgumentParser()
ap.add_argument("--survivors", required=True)
ap.add_argument("--tables-dir", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--threads", type=int, default=16)
ap.add_argument("--memory-limit", default="48GB")
ap.add_argument("--tmp-dir", default=None)
a = ap.parse_args()
t = a.tables_dir
c = duckdb.connect()
c.execute(f"pragma threads={a.threads}")
c.execute(f"pragma memory_limit='{a.memory_limit}'")
if a.tmp_dir:
    c.execute(f"pragma temp_directory='{a.tmp_dir}'")
c.execute(f"""create table s as select rep, n_new_proteins, n_genomes, split_part(rep, '__', 3) pid
  from read_csv('{a.survivors}', delim='\t', header=true)""")
c.execute(f"""create table pf as select p.protein_id pid, arg_max(split_part(p.pfam_acc, '.', 1), p.domain_score) best_pfam, count(*) n_dom
  from read_parquet('{t}/pfam.parquet') p semi join s on s.pid = p.protein_id group by 1""")
c.execute(f"""create table cz as select distinct h.protein_id pid from read_parquet('{t}/cazy.cazymes_hmm.parquet') h semi join s on s.pid = h.protein_id""")
c.execute(f"""create table mp as select distinct m.protein_id pid from read_parquet('{t}/merops.parquet') m semi join s on s.pid = m.protein_id""")
c.execute("""create table r as select s.rep, s.n_new_proteins, s.n_genomes, (pf.pid is not null) has_pfam, (cz.pid is not null) has_cazy,
  (mp.pid is not null) has_merops, pf.best_pfam from s left join pf using (pid) left join cz using (pid) left join mp using (pid)""")
c.execute(f"copy r to '{a.outdir}/survivor_annot.tsv' (delimiter '\t', header)")
c.execute(f"""copy (select case when n_genomes = 1 then '1' when n_genomes < 5 then '2-4' when n_genomes < 20 then '5-19' else '>=20' end genomes_class,
  count(*) clusters, sum(n_new_proteins) proteins,
  round(avg(has_pfam::int), 4) frac_clusters_pfam, round(sum(n_new_proteins * has_pfam::int) / sum(n_new_proteins), 4) frac_proteins_pfam,
  round(avg(has_cazy::int), 4) frac_cazy, round(avg(has_merops::int), 4) frac_merops,
  round(avg((has_pfam or has_cazy or has_merops)::int), 4) frac_any, round(avg((not (has_pfam or has_cazy or has_merops))::int), 4) frac_unannotated
  from r group by 1 order by 1) to '{a.outdir}/survivor_annot_summary.tsv' (delimiter '\t', header)""")
c.execute(f"""copy (select best_pfam, count(*) clusters, sum(n_new_proteins) proteins from r where has_pfam and n_genomes >= 5
  group by 1 order by 3 desc limit 100) to '{a.outdir}/survivor_top_pfam.tsv' (delimiter '\t', header)""")
print(c.execute(f"select * from read_csv('{a.outdir}/survivor_annot_summary.tsv', delim='\t', header=true)").fetchdf().to_string(index=False))
