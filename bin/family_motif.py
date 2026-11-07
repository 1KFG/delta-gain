#!/usr/bin/env python3
"""F5b: catalytic-motif integrity and domain architecture of novel vs known family regions (GH18 chitinase, S08 subtilase).

Novel = region's c70 cluster has no baseline member (new genomes only). Motifs are APPROXIMATE regexes for the canonical
catalytic sites; a missing motif suggests a degenerate/pseudo-enzyme, fragment, or model error - it is not proof.
  GH18: DGxDxDWE-like strand-4 motif  [LIVMF]D[GA][LIVMF]D[LIVMF][DE][WF]E
  S08 : Asp  [LIVMF][LIVMF]D[ST]G   His  H[GA][ST]H   Ser  G[TS]S[MA]A[AS]P   (all three = triad intact)
Architecture = other Pfam domains on the same protein (BFD pfam table), excluding the focal Pfam."""
import argparse
import re
import duckdb

MOT = {"cazy_family_GH18": {"cat": r"[LIVMF]D[GA][LIVMF]D[LIVMF][DE][WF]E"},
       "merops_family_S08": {"asp": r"[LIVMF][LIVMF]D[ST]G", "his": r"H[GA][ST]H", "ser": r"G[TS]S[MA]A[AS]P"}}
FOCAL = {"cazy_family_GH18": "PF00704", "merops_family_S08": "PF00082"}
ap = argparse.ArgumentParser()
ap.add_argument("--regions-dir", required=True); ap.add_argument("--clusters-dir", required=True)
ap.add_argument("--tables-dir", required=True); ap.add_argument("--outdir", required=True)
ap.add_argument("--tier", default="c70")
a = ap.parse_args()
c = duckdb.connect(); c.execute("pragma threads=8")
for fam, mot in MOT.items():
    pats = {k: re.compile(v) for k, v in mot.items()}
    seq = {}
    rid = None
    with open(f"{a.regions_dir}/{fam}.faa") as fh:
        for l in fh:
            if l.startswith(">"):
                rid = l[1:].strip(); seq[rid] = []
            else:
                seq[rid].append(l.strip())
    rows = []
    for r, parts in seq.items():
        sq = "".join(parts)
        rows.append((r,) + tuple(int(bool(p.search(sq))) for p in pats.values()))
    cols = ["region_id"] + list(pats)
    c.execute(f"create or replace table mo({', '.join(x + (' varchar' if i == 0 else ' int') for i, x in enumerate(cols))})")
    c.executemany(f"insert into mo values ({','.join('?' * len(cols))})", rows)
    c.execute(f"""create or replace table j as select m.region_id, m.protein_id, replace(m."class", chr(13), '') cls,
          try_cast(m.region_len as int) rlen, try_cast(m.model_cov as double) mcov, k.{a.tier} cid
      from read_csv('{a.regions_dir}/{fam}.members.tsv', delim='\t', header=true, all_varchar=true) m
      join read_csv('{a.clusters_dir}/{fam}.clusters.tsv', delim='\t', header=true, all_varchar=true) k using (region_id)""")
    c.execute("create or replace table cl as select cid, max((cls='baseline')::int) hb from j group by 1")
    allm = " and ".join(f"mo.{k}=1" for k in pats)
    print("==", fam, "motifs", mot)
    print(c.execute(f"""select case when j.cls='baseline' then 'baseline' when cl.hb=0 then 'new, novel (c70)' else 'new, known cluster' end grp,
        count(*) regions, round(avg(({allm})::int),3) all_motifs_intact, round(avg(rlen)) mean_len from j join cl using (cid) join mo using (region_id)
        group by 1 order by 1""").fetchdf().to_string(index=False))
    if len(pats) > 1:
        for k in pats:
            print(c.execute(f"select case when j.cls='baseline' then 'baseline' when cl.hb=0 then 'new, novel' else 'new, known' end grp, round(avg(mo.{k}::int),3) {k} from j join cl using (cid) join mo using (region_id) group by 1 order by 1").fetchdf().to_string(index=False))
    # architecture: co-domains
    c.execute(f"""create or replace table cod as select p.protein_id pid, list_sort(list(distinct split_part(p.pfam_acc, '.', 1))) doms
      from read_parquet('{a.tables_dir}/pfam.parquet') p semi join j on j.protein_id = p.protein_id
      where split_part(p.pfam_acc, '.', 1) <> '{FOCAL[fam]}' group by 1""")
    c.execute("create or replace table arch as select j.region_id, j.cls, cl.hb, (cod.pid is not null) has_co, cod.doms from j join cl using (cid) left join cod on cod.pid = j.protein_id")
    print(c.execute("""select case when cls='baseline' then 'baseline' when hb=0 then 'new, novel' else 'new, known' end grp, count(*) regions,
        round(avg(has_co::int),3) frac_with_other_pfam from arch group by 1 order by 1""").fetchdf().to_string(index=False))
    print(c.execute("""with u as (select unnest(doms) d, hb, cls from arch where doms is not null and cls='new'),
        t as (select d, sum((hb=0)::int) n_novel, sum((hb=1)::int) n_known from u group by 1)
        select d co_domain, n_novel, n_known, round(n_novel/(n_novel+n_known),3) frac_novel from t where n_novel+n_known >= 200 order by 4 desc limit 10""").fetchdf().to_string(index=False))
