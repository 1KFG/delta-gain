#!/usr/bin/env python3
"""F0 family census: per-family presence/copy number in baseline vs new genomes.

One DuckDB pass over the BFD parquet tables (Pfam, CAZy/dbCAN, MEROPS) joined to
genome classes and taxonomy.  See docs/FAMILY_DIVERSITY_PLAN.md (stage F0).

Levels written (column `level`):
  pfam               Pfam accession (version stripped), GA-gated hits as run
  cazy_family        dbCAN HMM family, subfamily suffix stripped (GH5_7 -> GH5)
  cazy_subfamily     dbCAN HMM ids that carry a subfamily suffix (GH5_7)
  cazy_efam          dbCAN-sub e-family from cazy.overview.sub_fam (GH27_e9)
  merops_family      MEROPS family (S01) from merops_lib.families.tab
  merops_subfamily   MEROPS subfamily (S01A)

Per family: genomes / proteins (members) / species, baseline vs new, copy number,
"new-taxon" species (species with no baseline genome), "gain within known species",
Fisher exact p-values (genome level and species level) and BH q-values per level.
Counts are over distinct genomes/proteins, never cluster representatives.

Only genomes with class in {baseline,new} and pep_status == matched are used.
Species level uses species.parquet SPECIES (fallback: ASMID).
"""
import argparse
import math
import os
import sys

import duckdb

LEVELS = {
    "pfam": ["pfam"],
    "cazy": ["cazy_family", "cazy_subfamily", "cazy_efam"],
    "merops": ["merops_family", "merops_subfamily"],
}
MIN_TESTED = 5  # families with fewer genomes than this get no p-value


def bh(pvals):
    """Benjamini-Hochberg q-values for a list that may contain None."""
    idx = [i for i, p in enumerate(pvals) if p is not None]
    out = [None] * len(pvals)
    m = len(idx)
    if m == 0:
        return out
    order = sorted(idx, key=lambda i: pvals[i])
    prev = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        prev = min(prev, pvals[i] * m / rank)
        out[i] = prev
    return out


def lfc(k_new, n_new, k_base, n_base):
    return math.log2(((k_new + 0.5) / (n_new + 1.0)) / ((k_base + 0.5) / (n_base + 1.0)))


def fisher_p(a, b, c, d):
    from scipy.stats import fisher_exact
    return float(fisher_exact([[a, b], [c, d]])[1])


def build_genomes(con, a):
    con.execute(f"""
      create table gen as
      select s.LOCUSTAG prefix, s.ASMID asmid,
             replace(c."class", chr(13), '') cls,
             coalesce(nullif(s.SPECIES, ''), s.ASMID) species
      from read_parquet('{a.tables_dir}/species.parquet') s
      join read_csv('{a.classes}', delim='\t', header=true, all_varchar=true) c
        on c.ASMID = s.ASMID
      where replace(c."class", chr(13), '') in ('baseline', 'new')
        and replace(c.pep_status, chr(13), '') = 'matched'
    """)
    if a.max_genomes:
        n = int(a.max_genomes)
        con.execute(f"""
          create table gen2 as
          select * from (select *, row_number() over (partition by cls order by hash(asmid)) rn from gen)
          where rn <= {n}""")
        con.execute("drop table gen")
        con.execute("create table gen as select prefix, asmid, cls, species from gen2")
    con.execute("""
      create table sp as
      select species, bool_or(cls = 'baseline') has_base, bool_or(cls = 'new') has_new
      from gen group by 1""")
    con.execute("create table gsel as select distinct prefix from gen")


def build_membership(con, a):
    con.execute("""create table fp (level varchar, fam varchar, name varchar,
                   prefix varchar, n_prot bigint, n_inst bigint, n_low bigint)""")
    t = a.tables_dir
    pre = "split_part(protein_id, '_', 1)"
    if "pfam" in a.sources:
        con.execute(f"""
          insert into fp
          select 'pfam', split_part(pfam_acc, '.', 1), any_value(pfam_id), {pre},
                 count(distinct protein_id), count(*), 0
          from read_parquet('{t}/pfam.parquet')
          where {pre} in (select prefix from gsel)
          group by split_part(pfam_acc, '.', 1), {pre}""")
    if "cazy" in a.sources:
        hmm = f"read_parquet('{t}/cazy.cazymes_hmm.parquet')"
        con.execute(f"""
          insert into fp
          select 'cazy_family', regexp_replace(HMM_id, '_[0-9]+$', ''), any_value(HMM_id), {pre},
                 count(distinct protein_id), count(*), 0
          from {hmm} where {pre} in (select prefix from gsel)
          group by regexp_replace(HMM_id, '_[0-9]+$', ''), {pre}""")
        con.execute(f"""
          insert into fp
          select 'cazy_subfamily', HMM_id, HMM_id, {pre},
                 count(distinct protein_id), count(*), 0
          from {hmm}
          where regexp_matches(HMM_id, '_[0-9]+$') and {pre} in (select prefix from gsel)
          group by HMM_id, {pre}""")
        con.execute(f"""
          insert into fp
          select 'cazy_efam', e, e, {pre},
                 count(distinct protein_id), count(*), 0
          from (select protein_id, unnest(string_split(sub_fam, '+')) e
                from read_parquet('{t}/cazy.overview.parquet')
                where sub_fam is not null and sub_fam not in ('-', ''))
          where e not in ('-', '') and {pre} in (select prefix from gsel)
          group by e, {pre}""")
    if "merops" in a.sources:
        con.execute(f"""
          create table mmap as
          select column0 mer, column1 subfam, split_part(column2, '.', 1) fam
          from read_csv('{a.merops_map}', delim='\t', header=false, all_varchar=true)""")
        con.execute(f"""
          create table mb as
          select m.protein_id, f.fam, f.subfam,
                 max(m.percent_identity) pid, min(m.evalue) ev
          from read_parquet('{t}/merops.parquet') m
          join mmap f on f.mer = m.merops_id
          where split_part(m.protein_id, '_', 1) in (select prefix from gsel)
          group by m.protein_id, f.fam, f.subfam""")
        # best hit per protein per family / subfamily (a protein can hit >1 family)
        for lvl, col in (("merops_family", "fam"), ("merops_subfamily", "subfam")):
            con.execute(f"""
              insert into fp
              select '{lvl}', {col}, {col}, split_part(protein_id, '_', 1),
                     count(distinct protein_id), count(*),
                     sum((pid < 40)::int)
              from mb group by {col}, split_part(protein_id, '_', 1)""")


def aggregate(con):
    con.execute("""
      create table j as
      select f.level, f.fam, f.name, g.asmid, g.cls, g.species, f.n_prot, f.n_inst, f.n_low
      from fp f join gen g on g.prefix = f.prefix""")
    con.execute("""
      create table agg_g as
      select level, fam, any_value(name) name,
        count(*) filter (cls = 'baseline') g_base, count(*) filter (cls = 'new') g_new,
        coalesce(sum(n_prot) filter (cls = 'baseline'), 0) prot_base,
        coalesce(sum(n_prot) filter (cls = 'new'), 0) prot_new,
        coalesce(sum(n_inst) filter (cls = 'baseline'), 0) inst_base,
        coalesce(sum(n_inst) filter (cls = 'new'), 0) inst_new,
        coalesce(sum(n_low) filter (cls = 'baseline'), 0) lowid_base,
        coalesce(sum(n_low) filter (cls = 'new'), 0) lowid_new
      from j group by level, fam""")
    con.execute("""
      create table agg_s as
      with sf as (select level, fam, species,
                    bool_or(cls = 'baseline') in_base, bool_or(cls = 'new') in_new
                  from j group by level, fam, species)
      select sf.level, sf.fam,
        count(*) filter (in_base) sp_base, count(*) filter (in_new) sp_new,
        count(*) filter (in_new and not sp.has_base) sp_newtaxon,
        count(*) filter (in_new and not in_base and sp.has_base) sp_gain_known
      from sf join sp on sp.species = sf.species group by sf.level, sf.fam""")
    return con.execute("""
      select g.*, s.sp_base, s.sp_new, s.sp_newtaxon, s.sp_gain_known
      from agg_g g join agg_s s using (level, fam) order by level, fam""").fetchall()


COLS = ["level", "family", "name", "g_base", "g_new", "prev_base", "prev_new",
        "prot_base", "prot_new", "copies_base", "copies_new", "inst_base", "inst_new",
        "lowid_base", "lowid_new", "sp_base", "sp_new", "sp_newtaxon", "sp_gain_known",
        "log2_prev_ratio_genome", "log2_prev_ratio_species", "p_genome", "p_species",
        "q_genome", "q_species", "new_only", "base_only"]


def fmt(v):
    if v is None:
        return "NA"
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--tables-dir", required=True)
    ap.add_argument("--classes", required=True, help="genome_classes.tsv (ASMID, class, pep_status)")
    ap.add_argument("--merops-map", default="/srv/projects/db/MEROPS/2023-01/merops_lib.families.tab")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--sources", nargs="+", default=["pfam", "cazy", "merops"],
                    choices=["pfam", "cazy", "merops"])
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--memory-limit", default="64GB")
    ap.add_argument("--tmp-dir", default=None)
    ap.add_argument("--max-genomes", type=int, default=0,
                    help="test mode: keep N genomes per class (hash order)")
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)

    con = duckdb.connect()
    con.execute(f"PRAGMA threads={a.threads}")
    con.execute(f"PRAGMA memory_limit='{a.memory_limit}'")
    if a.tmp_dir:
        con.execute(f"PRAGMA temp_directory='{a.tmp_dir}'")

    build_genomes(con, a)
    n_gb, n_gn = con.execute(
        "select count(*) filter (cls='baseline'), count(*) filter (cls='new') from gen").fetchone()
    n_sb, n_sn, n_snt = con.execute(
        "select count(*) filter (has_base), count(*) filter (has_new), "
        "count(*) filter (has_new and not has_base) from sp").fetchone()
    print(f"genomes baseline={n_gb} new={n_gn}; species with baseline={n_sb} with new={n_sn} "
          f"new-taxon species={n_snt}", flush=True)
    with open(os.path.join(a.outdir, "census_totals.tsv"), "w") as fh:
        fh.write("metric\tvalue\n")
        for k, v in (("genomes_baseline", n_gb), ("genomes_new", n_gn), ("species_with_baseline", n_sb),
                     ("species_with_new", n_sn), ("species_new_taxon", n_snt)):
            fh.write(f"{k}\t{v}\n")

    build_membership(con, a)
    rows = aggregate(con)
    print(f"{len(rows)} family rows", flush=True)

    by_level = {}
    for r in rows:
        (level, fam, name, gb, gn, pb, pn, ib, in_, lb, ln, sb, sn, snt, sgk) = r
        pg = ps = None
        if gb + gn >= MIN_TESTED:
            pg = fisher_p(gb, n_gb - gb, gn, n_gn - gn)
            ps = fisher_p(sb, n_sb - sb, sn, n_sn - sn)
        rec = [level, fam, name, gb, gn, gb / n_gb if n_gb else 0, gn / n_gn if n_gn else 0,
               pb, pn, pb / gb if gb else None, pn / gn if gn else None, ib, in_, lb, ln,
               sb, sn, snt, sgk,
               lfc(gn, n_gn, gb, n_gb), lfc(sn, n_sn, sb, n_sb), pg, ps, None, None,
               int(gb == 0 and gn >= 3), int(gn == 0 and gb >= 3)]
        by_level.setdefault(level, []).append(rec)
    ip, is_ = COLS.index("p_genome"), COLS.index("p_species")
    iq, iqs = COLS.index("q_genome"), COLS.index("q_species")
    for recs in by_level.values():
        qg = bh([r[ip] for r in recs])
        qs = bh([r[is_] for r in recs])
        for r, g, s in zip(recs, qg, qs):
            r[iq], r[iqs] = g, s

    for src, levels in LEVELS.items():
        recs = [r for l in levels for r in by_level.get(l, [])]
        if not recs:
            continue
        path = os.path.join(a.outdir, f"{src}_census.tsv")
        with open(path, "w") as fh:
            fh.write("\t".join(COLS) + "\n")
            for r in recs:
                fh.write("\t".join(fmt(v) for v in r) + "\n")
        print("wrote", path, len(recs), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
