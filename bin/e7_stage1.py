#!/usr/bin/env python3
"""E7 stage 1 (DuckDB): genome metadata, SSP candidates (before the cysteine filter), per-genome counts, CAZy/MEROPS family copy numbers.

Candidate = SignalP signal peptide (probability >= 0.5, best row per protein) + length <= 300 aa + no GPI-anchor (PredGPI) +
no TM helix beyond a signal anchor (TMHMM PredHel = 0, or PredHel = 1 and the helix starts <= peptide_end + 10; proteins absent
from the TMHMM table count as 0 helices). The >= 4 cysteine rule is applied in stage 2/3 once sequences are read.
"""
import argparse
import os
import sys

import duckdb


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables-dir", required=True)
    ap.add_argument("--classes", required=True)
    ap.add_argument("--merops-map", default="/srv/projects/db/MEROPS/2023-01/merops_lib.families.tab")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--memory-limit", default="120GB")
    ap.add_argument("--tmp-dir", default=None)
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    t = a.tables_dir
    o = a.outdir
    c = duckdb.connect()
    c.execute(f"pragma threads={a.threads}")
    c.execute(f"pragma memory_limit='{a.memory_limit}'")
    if a.tmp_dir:
        c.execute(f"pragma temp_directory='{a.tmp_dir}'")
    # genome metadata: one row per genome with BFD annotation prefix
    c.execute(f"""create table meta as
      select cls.ASMID asmid, sp.prefix, cls.class cl, coalesce(sp.species, cls.ASMID) species, sp."ord" ord
      from read_csv('{a.classes}', delim='\t', header=true) cls
      join (select ASMID, any_value(LOCUSTAG) prefix, any_value(SPECIES) species, any_value("ORDER") "ord"
            from read_parquet('{t}/species.parquet') group by 1) sp on sp.ASMID = cls.ASMID
      where sp.prefix is not null""")
    c.execute(f"copy meta to '{o}/genome_meta.parquet' (format parquet)")
    print("genomes with prefix", c.execute("select count(*), count(distinct prefix) from meta").fetchone(), flush=True)
    c.execute(f"""create table sp as
      select protein_id, any_value(species_prefix) prefix, max(probability) prob, arg_max(peptide_end, probability) pend
      from read_parquet('{t}/signalp.signal_peptide.parquet') where probability >= 0.5 group by protein_id""")
    print("signal peptide proteins", c.execute("select count(*) from sp").fetchone(), flush=True)
    c.execute(f"""create table tm as
      select protein_id, max(PredHel) nhel, any_value(try_cast(regexp_extract(Topology, '(\\d+)-', 1) as int)) first_start
      from read_parquet('{t}/tmhmm.parquet') group by protein_id""")
    print("tmhmm rows", c.execute("select count(*), count(*) filter (nhel = 0) from tm").fetchone(), flush=True)
    c.execute(f"""create table gpi as
      select distinct protein_id from read_parquet('{t}/predgpi.parquet') where feature = 'GPI-anchor'""")
    c.execute(f"""create table len as
      select protein_id, any_value(length) plen from read_parquet('{t}/gene_proteins.parquet') group by protein_id""")
    c.execute("""create table sp_flag as
      select sp.protein_id, sp.prefix, sp.pend, len.plen,
        (len.plen <= 300 and gpi.protein_id is null
         and (tm.protein_id is null or tm.nhel = 0 or (tm.nhel = 1 and tm.first_start <= sp.pend + 10)))::int cand
      from sp left join len using (protein_id) left join tm using (protein_id) left join gpi using (protein_id)""")
    c.execute(f"copy sp_flag to '{o}/sp_flag.parquet' (format parquet)")
    c.execute(f"copy (select * from sp_flag where cand = 1) to '{o}/cand.parquet' (format parquet)")
    print("sp / candidates", c.execute("select count(*), sum(cand) from sp_flag").fetchone(), flush=True)
    c.execute(f"""create table cnt as
      select split_part(protein_id, '_', 1) prefix, count(*) n_prot
      from read_parquet('{t}/gene_proteins.parquet') group by 1""")
    c.execute("""create table cnt2 as
      select cnt.prefix, cnt.n_prot, coalesce(s.n_sp, 0) n_sp, coalesce(s.n_cand, 0) n_cand
      from cnt left join (select prefix, count(*) n_sp, sum(cand) n_cand from sp_flag group by 1) s using (prefix)""")
    c.execute(f"copy cnt2 to '{o}/genome_counts.parquet' (format parquet)")
    # family copy numbers
    c.execute(f"""create table cz as
      select split_part(protein_id, '_', 1) prefix, 'cazy' kind, regexp_replace(HMM_id, '_[0-9]+$', '') fam, count(distinct protein_id) n
      from read_parquet('{t}/cazy.cazymes_hmm.parquet') group by all""")
    c.execute(f"""create table mmap as
      select column0 mer, split_part(column2, '.', 1) fam
      from read_csv('{a.merops_map}', delim='\t', header=false, all_varchar=true)""")
    c.execute(f"""create table mp as
      select split_part(m.protein_id, '_', 1) prefix, 'merops' kind, f.fam, count(distinct m.protein_id) n
      from read_parquet('{t}/merops.parquet') m join mmap f on f.mer = m.merops_id group by all""")
    c.execute(f"copy (select * from cz union all select * from mp) to '{o}/fam_counts.parquet' (format parquet)")
    print("family rows", c.execute("select kind, count(distinct fam), count(*) from (select * from cz union all select * from mp) group by 1").fetchall(), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
