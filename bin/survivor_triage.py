#!/usr/bin/env python3
"""E5: contamination / transposable-element triage of funnel survivors.
(a) genome level: GC outlier genomes (|GC - order median| > --gc-dev) and low-BUSCO genomes: their share of survivor proteins and
    their novelty compared with the rest of the same stratum. (b) cluster level: share of survivor clusters/proteins whose
    representative's best Pfam is a transposon/retroelement domain, by genomes class."""
import argparse
import duckdb

ap = argparse.ArgumentParser()
ap.add_argument("--results", required=True)
ap.add_argument("--tables-dir", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--gc-dev", type=float, default=8.0)
a = ap.parse_args()
R, T = a.results, a.tables_dir
c = duckdb.connect()
c.execute(f"""create table g as select s.*, x.gc_pct gc, x.contig_count contigs, x.N50_bp n50, b.complete_pct busco
  from read_csv('{R}/delta_funnel/summary/genome_stratum.tsv', delim='\t', header=true) s
  left join read_parquet('{T}/asm_stats.parquet') x on x.ASMID = s.genome
  left join read_parquet('{T}/busco_genome.parquet') b on b.ASMID = s.genome""")
c.execute("""create table g2 as select g.*, median(gc) over (partition by ord) ord_gc, count(*) over (partition by ord) ord_n from g""")
c.execute(f"""create table g3 as select *, (abs(gc - ord_gc) > {a.gc_dev} and ord_n >= 5) gc_out, (busco < 50) low_busco from g2""")
print("== genome flags")
q = """select stratum, count(*) genomes, sum(gc_out::int) gc_outliers, sum(low_busco::int) busco_lt50,
  round(sum(novel_50 * gc_out::int) / sum(novel_50), 4) share_novel50_gc_out, round(sum(novel_50 * low_busco::int) / sum(novel_50), 4) share_novel50_lowbusco,
  round(avg(case when gc_out then frac_50 end), 3) frac50_gc_out, round(avg(case when not gc_out then frac_50 end), 3) frac50_other
  from g3 group by 1 order by 1"""
df = c.execute(q).fetchdf(); print(df.to_string(index=False)); df.to_csv(f"{a.outdir}/triage_genome_flags.tsv", sep="\t", index=False)
print(c.execute("select count(*) genomes, sum(gc_out::int) gc_out, sum(low_busco::int) low_busco, round(sum(novel_50*(gc_out or low_busco)::int)/sum(novel_50),4) share_novel50_flagged from g3").fetchdf().to_string(index=False))
c.execute(f"copy (select genome, stratum, ord, species, new_proteins, novel_50, frac_50, gc, ord_gc, busco, contigs, n50, gc_out from g3 where gc_out or frac_50 > 0.9 order by novel_50 desc limit 300) to '{a.outdir}/triage_top_genomes.tsv' (delimiter '\t', header)")
print("== TE / retroelement Pfam among survivors")
c.execute(f"""create table te as select family from read_csv('{R}/family_census/pfam_census.tsv', delim='\t', header=true, all_varchar=true)
  where regexp_matches(name, '(?i)(transpos|^RVT|retro|Tc1|^hAT|DDE|gypsy|helitron|^rve|integrase|Tnp|^gag|polyprotein|LTR|Copia|Ty3|mariner|piggy|CENP-B|Dimer_Tnp)')""")
print(c.execute("select count(*) te_pfams from te").fetchall())
c.execute(f"""create table sa as select * from read_csv('{R}/survivor_annot/survivor_annot.tsv', delim='\t', header=true)""")
df = c.execute("""select case when n_genomes = 1 then '1' when n_genomes < 5 then '2-4' when n_genomes < 20 then '5-19' else '>=20' end genomes_class,
  count(*) clusters, sum(n_new_proteins) proteins, sum(has_pfam::int) with_pfam,
  round(sum((sa.best_pfam in (select family from te))::int) / sum(has_pfam::int), 4) te_frac_of_pfam_clusters,
  round(sum(case when sa.best_pfam in (select family from te) then n_new_proteins end) / sum(n_new_proteins), 4) te_frac_of_proteins
  from sa group by 1 order by 1""").fetchdf(); print(df.to_string(index=False)); df.to_csv(f"{a.outdir}/triage_te.tsv", sep="\t", index=False)
