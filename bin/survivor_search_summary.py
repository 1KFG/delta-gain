#!/usr/bin/env python3
"""Summarise the survivor search: best hit per query representative, weighted by cluster size.
Classes: no_hit | distant (best e in (1e-10, 1e-5]) | close (<= 1e-10 and pident < 40) | homolog (<=1e-10, pident >= 40).
Outputs survivor_search_best.tsv (per query) and survivor_search_summary.tsv (counts weighted by n_new_proteins / n clusters)."""
import argparse, duckdb
ap = argparse.ArgumentParser()
ap.add_argument("--hits", required=True); ap.add_argument("--survivors", required=True)
ap.add_argument("--query-fasta", required=True); ap.add_argument("--outdir", required=True)
a = ap.parse_args()
q = []
with open(a.query_fasta) as fh:
    for l in fh:
        if l.startswith(">"):
            q.append(l[1:].split()[0])
c = duckdb.connect(); c.execute("pragma threads=16")
c.execute("create table q(rep varchar)"); c.executemany("insert into q values (?)", [(x,) for x in q])
c.execute(f"""create table surv as select rep, n_new_proteins, n_genomes from read_csv('{a.survivors}',delim='\t',header=true)""")
c.execute(f"""create table best as select query rep, arg_min(target, evalue) target, min(evalue) evalue, arg_max(pident, bits) pident, max(bits) bits, max(qcov) qcov
  from read_csv('{a.hits}',delim='\t',header=false,columns={{'query':'VARCHAR','target':'VARCHAR','pident':'DOUBLE','evalue':'DOUBLE','bits':'DOUBLE','qcov':'DOUBLE','tcov':'DOUBLE'}}) group by 1""")
c.execute("""create table r as select q.rep, s.n_new_proteins, s.n_genomes, b.target, b.evalue, b.pident, b.qcov,
  case when b.rep is null then 'no_hit' when b.evalue > 1e-10 then 'distant' when b.pident < 40 then 'close_lowid' else 'homolog' end cls
  from q join surv s using (rep) left join best b using (rep)""")
c.execute(f"copy r to '{a.outdir}/survivor_search_best.tsv' (delimiter '\t', header)")
c.execute(f"""copy (select cls, count(*) clusters, sum(n_new_proteins) proteins, round(count(*)/sum(count(*)) over(),4) frac_clusters,
   round(sum(n_new_proteins)/sum(sum(n_new_proteins)) over(),4) frac_proteins from r group by 1 order by 2 desc)
   to '{a.outdir}/survivor_search_summary.tsv' (delimiter '\t', header)""")
print(c.execute("select * from read_csv('%s/survivor_search_summary.tsv',delim='\t',header=true)" % a.outdir).fetchdf().to_string(index=False))
