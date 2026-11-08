#!/usr/bin/env python3
"""HGT screen proxy: taxonomic origin of the best Swiss-Prot hit of novel (survivor) proteins vs background proteins.

For each protein with a Swiss-Prot hit (e <= --evalue, query coverage >= --qcov), the best hit (max bitscore) organism is
mapped to a broad group with the NCBI taxdump (Bacteria/Archaea = prokaryote, Viruses, Fungi, Metazoa, Viridiplantae, other
Eukaryota). Compared sets: survivor representatives (by number of genomes) vs a hash-sampled background of all BFD proteins
that are not survivors. Metagenome-bin genomes (ASMID containing '_bin.') are excluded because contamination, not HGT, would
dominate there. This is a SCREEN: a prokaryote best hit is compatible with HGT, contamination, or an unsampled-clade
artefact; no phylogeny is built.
"""
import argparse
import io
import os
import re
import sys
import tarfile
import urllib.request

import duckdb

URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"


def fetch_taxa(outdir):
    path = os.path.join(outdir, "taxa_group.tsv")
    if os.path.exists(path):
        return path
    local = os.path.join(outdir, "new_taxdump.tar.gz")
    if not os.path.exists(local):
        urllib.request.urlretrieve(URL, local)
    tf = tarfile.open(local, mode="r:gz")
    rl = io.TextIOWrapper(tf.extractfile("rankedlineage.dmp"), encoding="utf-8")
    seen = {}
    gen = {}
    for line in rl:
        if not line.strip():
            continue
        f = [x.strip() for x in line.strip().rstrip("|").strip().split("\t|\t")]
        # tax_id, name, species, genus, family, order, class, phylum, kingdom, superkingdom
        if len(f) < 10:
            continue
        name, genus, kingdom, sk = f[1], f[3], f[8], f[9]
        if sk in ("Bacteria", "Archaea"):
            g = "prokaryote"
        elif sk == "Viruses":
            g = "virus"
        elif kingdom in ("Fungi", "Metazoa", "Viridiplantae"):
            g = kingdom
        elif sk == "Eukaryota":
            g = "other_euk"
        else:
            continue
        seen.setdefault(name.lower(), g)
        if genus:
            gen.setdefault(genus.lower(), g)
    with open(path, "w") as fh:
        fh.write("key\tlevel\tgrp\n")
        for k, g in seen.items():
            fh.write(f"{k}\tname\t{g}\n")
        for k, g in gen.items():
            fh.write(f"{k}\tgenus\t{g}\n")
    return path


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--survivors", required=True, help="survivor_annot.tsv (rep, n_new_proteins, n_genomes, ...)")
    ap.add_argument("--tables-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--qcov", type=float, default=0.5, help="query_cov is a fraction in the BFD table")
    ap.add_argument("--bg-mod", type=int, default=50, help="background = hash(protein_id) %% N == 0")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--memory-limit", default="90GB")
    ap.add_argument("--tmp-dir", default=None)
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    tx = fetch_taxa(a.outdir)
    t = a.tables_dir
    c = duckdb.connect()
    c.execute(f"pragma threads={a.threads}")
    c.execute(f"pragma memory_limit='{a.memory_limit}'")
    if a.tmp_dir:
        c.execute(f"pragma temp_directory='{a.tmp_dir}'")
    c.execute(f"create table tx as select * from read_csv('{tx}', delim='\t', header=true)")
    c.execute("create table tname as select key, any_value(grp) grp from tx where level = 'name' group by 1")
    c.execute("create table tgen as select key, any_value(grp) grp from tx where level = 'genus' group by 1")
    # accession -> group via the cleaned organism string
    c.execute(f"""create table acc as
      select accession, coalesce(tn.grp, tg.grp, 'unmapped') grp, any_value(organism) organism
      from (select accession, organism, lower(trim(regexp_replace(organism, '\\s*\\(.*$', ''))) nm,
                   lower(split_part(trim(regexp_replace(organism, '\\s*\\(.*$', '')), ' ', 1)) gn
            from read_parquet('{t}/swissprot_annot.parquet')) s
      left join tname tn on tn.key = s.nm left join tgen tg on tg.key = s.gn
      group by accession, tn.grp, tg.grp""")
    print(c.execute("select grp, count(*) from acc group by 1 order by 2 desc").fetchall(), flush=True)
    c.execute(f"""create table best as
      select protein_id pid, arg_max(swissprot_acc, bitscore) acc, max(bitscore) bs, arg_max(pident, bitscore) pident
      from read_parquet('{t}/swissprot.parquet')
      where try_cast(evalue as double) <= {a.evalue} and try_cast(query_cov as double) >= {a.qcov}
      group by 1""")
    c.execute(f"""create table s as select split_part(rep, '__', 2) asm, split_part(rep, '__', 3) pid, n_genomes
      from read_csv('{a.survivors}', delim='\t', header=true) where split_part(rep, '__', 2) not like '%\\_bin.%' escape '\\'""")
    c.execute("""create table setA as select case when n_genomes = 1 then '1' when n_genomes < 5 then '2-4' else '>=5' end cls,
      b.pident, coalesce(ac.grp, 'unmapped') grp, true has_hit from s join best b using (pid) left join acc ac on ac.accession = b.acc""")
    c.execute("create table nS as select case when n_genomes = 1 then '1' when n_genomes < 5 then '2-4' else '>=5' end cls, count(*) n from s group by 1")
    c.execute(f"""create table setB as select 'background' cls, b.pident, coalesce(ac.grp, 'unmapped') grp
      from best b left join acc ac on ac.accession = b.acc
      where hash(b.pid) % {a.bg_mod} = 0 and b.pid not in (select pid from s)""")
    out = os.path.join(a.outdir, "hgt_proxy.tsv")
    c.execute(f"""copy (
      with u as (select cls, pident, grp from setA union all select cls, pident, grp from setB)
      select cls, count(*) with_hit, round(avg((grp = 'prokaryote')::int), 4) frac_prokaryote,
        round(avg((grp = 'virus')::int), 4) frac_virus, round(avg((grp = 'Fungi')::int), 4) frac_fungi,
        round(avg((grp = 'Metazoa')::int), 4) frac_metazoa, round(avg((grp = 'Viridiplantae')::int), 4) frac_plant,
        round(avg((grp = 'other_euk')::int), 4) frac_other_euk, round(avg((grp = 'unmapped')::int), 4) frac_unmapped,
        round(avg((grp = 'prokaryote' and pident >= 40)::int), 4) frac_prok_pid40,
        round(avg((grp = 'prokaryote' and pident >= 60)::int), 4) frac_prok_pid60
      from u group by 1 order by 1) to '{out}' (delimiter '\t', header)""")
    c.execute(f"copy (select * from nS order by cls) to '{a.outdir}/hgt_survivor_counts.tsv' (delimiter '\t', header)")
    print(open(out).read())
    return 0


if __name__ == "__main__":
    sys.exit(main())
