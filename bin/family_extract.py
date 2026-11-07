#!/usr/bin/env python3
"""F1 family region extraction: one pass over the per-genome BFD protein files.

For each requested family writes
  <outdir>/<key>.faa            region sequences  (>protein_id|start-end)
  <outdir>/<key>.members.tsv    region_id, protein_id, asmid, class, species, start, end,
                                region_len, protein_len, score, evalue, model_cov, sub,
                                pident, efam, ec, substrate
  <outdir>/regions_summary.tsv  per-family counts

family keys:  pfam:PF00227   cazy_family:GH18   merops_family:S08   (file key uses '_' for ':')
Regions
  pfam            one region per domain instance, env_from..env_to; model_cov = hmm span / pfam_len
  cazy_family     s_start..s_end of dbCAN HMM hits; overlapping hits of one protein/family are merged
                  (best-E hit supplies sub = HMM id, plus efam/EC/substrate from cazy.overview)
  merops_family   q_start..q_end of MEROPS BLAST hits mapped to the family; overlapping hits merged
                  (best hit supplies sub = MER id and pident)
No identity/coverage filtering is applied here (Pfam GA and dbCAN/MEROPS cuts were applied upstream);
downstream stages threshold on model_cov / score.  Only genomes with class in {baseline,new}
and pep_status == matched are used (same selection as bin/family_census.py).
"""
import argparse
import csv
import gzip
import os
import sys
from collections import defaultdict
import multiprocessing as mp

import duckdb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import family_census  # noqa: E402

DEFAULT_FAMILIES = ["pfam:PF00227", "pfam:PF00091", "pfam:PF00704", "pfam:PF03443",
                    "cazy_family:GH18", "cazy_family:AA9",
                    "merops_family:S08", "merops_family:M14"]
MIN_REGION = 30
HEADER = ["region_id", "protein_id", "asmid", "class", "species", "start", "end", "region_len",
          "protein_len", "score", "evalue", "model_cov", "sub", "pident", "efam", "ec", "substrate"]


def fam_key(spec):
    return spec.replace(":", "_")


def merge_hits(rows):
    """rows: (protein_id, start, end, evalue, rest) ; merge overlapping per protein, keep best-E."""
    rows = sorted(rows, key=lambda r: (r[0], r[1], r[2]))
    out, cur = [], None
    for r in rows:
        if cur and r[0] == cur[0] and r[1] <= cur[2]:
            cur[2] = max(cur[2], r[2])
            if r[3] < cur[3]:
                cur[3], cur[4] = r[3], r[4]
        else:
            if cur:
                out.append(tuple(cur))
            cur = [r[0], r[1], r[2], r[3], r[4]]
    if cur:
        out.append(tuple(cur))
    return out


def fetch_members(con, a, specs):
    """Return {famkey: [(protein_id, start, end, meta dict)]} for requested families."""
    t = a.tables_dir
    pre = "split_part(protein_id, '_', 1)"
    by = defaultdict(list)
    pf = [s.split(":", 1)[1] for s in specs if s.startswith("pfam:")]
    if pf:
        ids = ",".join(f"'{x}'" for x in pf)
        q = con.execute(f"""
          select split_part(pfam_acc, '.', 1), protein_id, env_from, env_to, domain_score,
                 domain_i_evalue, hmm_from, hmm_to, pfam_len
          from read_parquet('{t}/pfam.parquet')
          where split_part(pfam_acc, '.', 1) in ({ids}) and {pre} in (select prefix from gsel)""")
        for fam, pid, s, e, sc, ev, hf, ht, pl in q.fetchall():
            by["pfam_" + fam].append((pid, int(s), int(e), dict(
                score=sc, evalue=ev, model_cov=round((ht - hf + 1) / pl, 3) if pl else None)))
    cf = [s.split(":", 1)[1] for s in specs if s.startswith("cazy_family:")]
    if cf:
        ids = ",".join(f"'{x}'" for x in cf)
        per = defaultdict(list)
        for fam, pid, hid, s, e, ev, cov, ec, efam, sub in con.execute(f"""
              select regexp_replace(h.HMM_id, '_[0-9]+$', ''), h.protein_id, h.HMM_id, h.s_start, h.s_end,
                     h.evalue, h.coverage, o.EC, o.sub_fam, o.substrate
              from read_parquet('{t}/cazy.cazymes_hmm.parquet') h
              left join read_parquet('{t}/cazy.overview.parquet') o on o.protein_id = h.protein_id
              where regexp_replace(h.HMM_id, '_[0-9]+$', '') in ({ids})
                and split_part(h.protein_id, '_', 1) in (select prefix from gsel)""").fetchall():
            per[fam].append((pid, int(s), int(e), float(ev), (hid, cov, ec, efam, sub)))
        for fam, rows in per.items():
            for pid, s, e, ev, (hid, cov, ec, efam, sub) in merge_hits(rows):
                by["cazy_family_" + fam].append((pid, s, e, dict(
                    evalue=ev, model_cov=cov, sub=hid, efam=efam, ec=ec, substrate=sub)))
    mf = [s.split(":", 1)[1] for s in specs if s.startswith("merops_family:")]
    if mf:
        ids = ",".join(f"'{x}'" for x in mf)
        con.execute(f"""create or replace table mmap as
          select column0 mer, column1 subfam, split_part(column2, '.', 1) fam
          from read_csv('{a.merops_map}', delim='\t', header=false, all_varchar=true)""")
        per = defaultdict(list)
        for fam, pid, mer, s, e, ev, pid_, bs in con.execute(f"""
              select f.fam, m.protein_id, m.merops_id, m.q_start, m.q_end, m.evalue, m.percent_identity, m.bitscore
              from read_parquet('{t}/merops.parquet') m join mmap f on f.mer = m.merops_id
              where f.fam in ({ids}) and split_part(m.protein_id, '_', 1) in (select prefix from gsel)"""
                                              ).fetchall():
            per[fam].append((pid, int(min(s, e)), int(max(s, e)), float(ev), (mer, pid_, bs)))
        for fam, rows in per.items():
            for pid, s, e, ev, (mer, pident, bs) in merge_hits(rows):
                by["merops_family_" + fam].append((pid, s, e, dict(
                    evalue=ev, score=bs, sub=mer, pident=pident)))
    return by


def read_fasta_wanted(path, wanted):
    opener = gzip.open if path.endswith(".gz") else open
    seqs, cur, buf = {}, None, []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None and cur in wanted:
                    seqs[cur] = "".join(buf).replace("*", "")
                tok = line[1:].split(None, 1)
                cur, buf = (tok[0] if tok else ""), []
            elif cur in wanted:
                buf.append(line.strip())
        if cur is not None and cur in wanted:
            seqs[cur] = "".join(buf).replace("*", "")
    return seqs


def work(task):
    prefix, path, rows = task  # rows: (famkey, pid, s, e, meta)
    try:
        seqs = read_fasta_wanted(path, {r[1] for r in rows})
    except OSError as ex:
        return prefix, None, str(ex)
    out, miss = [], 0
    for famkey, pid, s, e, meta in rows:
        sq = seqs.get(pid)
        if sq is None:
            miss += 1
            continue
        e = min(e, len(sq))
        if e - s + 1 < MIN_REGION:
            continue
        out.append((famkey, pid, s, e, meta, sq[s - 1:e], len(sq)))
    return prefix, out, miss


def read_manifest(path):
    paths = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            paths[r["ASMID"].strip()] = r["pep_path"].strip()
    return paths


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--tables-dir", required=True)
    ap.add_argument("--classes", required=True)
    ap.add_argument("--manifest", required=True, help="bfd_pep_manifest.tsv (ASMID, pep_path)")
    ap.add_argument("--merops-map", default="/srv/projects/db/MEROPS/2023-01/merops_lib.families.tab")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--families", nargs="+", default=DEFAULT_FAMILIES)
    ap.add_argument("--census-dir", default=None, help="add --auto-duf families from pfam_census.tsv")
    ap.add_argument("--auto-duf", type=int, default=0)
    ap.add_argument("--auto-pfam", type=int, default=0,
                    help="top N Pfam families (any name) by species-level prevalence ratio in the census "
                         "(g_new >= 100, q_species < 0.01)")
    ap.add_argument("--auto-levels", nargs="*", default=[],
                    help="add EVERY family of these census levels, e.g. cazy_family merops_family")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--memory-limit", default="32GB")
    ap.add_argument("--tmp-dir", default=None)
    ap.add_argument("--max-genomes", type=int, default=0)
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    specs = list(a.families)

    if a.auto_duf:
        # Pfam DUF families most enriched in new genomes (species-level), well supported
        rows = []
        with open(os.path.join(a.census_dir, "pfam_census.tsv")) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["name"].startswith("DUF") and r["q_species"] != "NA" and int(r["g_new"]) >= 100 \
                        and float(r["q_species"]) < 0.01:
                    rows.append((float(r["log2_prev_ratio_species"]), r["family"], r["name"]))
        rows.sort(reverse=True)
        duf = [f"pfam:{f}" for _, f, _n in rows[:a.auto_duf]]
        print("auto DUF families:", [(f, n, round(l, 2)) for l, f, n in rows[:a.auto_duf]], flush=True)
        specs += [s for s in duf if s not in specs]

    if a.auto_pfam:
        # NB: ranking by enrichment in new genomes selects the best case for 'new adds diversity'
        rows = []
        with open(os.path.join(a.census_dir, "pfam_census.tsv")) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["q_species"] != "NA" and int(r["g_new"]) >= 100 and float(r["q_species"]) < 0.01:
                    rows.append((float(r["log2_prev_ratio_species"]), r["family"], r["name"]))
        rows.sort(reverse=True)
        top = [f"pfam:{f}" for _, f, _n in rows[:a.auto_pfam]]
        print(f"auto Pfam: {len(top)} of {len(rows)} eligible; top5", rows[:5], flush=True)
        specs += [s for s in top if s not in specs]
    for lev in a.auto_levels:
        fams = []
        with open(os.path.join(a.census_dir, lev.split("_")[0] + "_census.tsv")) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["level"] == lev:
                    fams.append(f"{lev}:{r['family']}")
        print(f"auto level {lev}: {len(fams)} families", flush=True)
        specs += [s for s in fams if s not in specs]

    con = duckdb.connect()
    con.execute(f"PRAGMA threads={a.threads}")
    con.execute(f"PRAGMA memory_limit='{a.memory_limit}'")
    if a.tmp_dir:
        con.execute(f"PRAGMA temp_directory='{a.tmp_dir}'")
    family_census.build_genomes(con, a)
    gen = {p: (asm, cls, sp) for p, asm, cls, sp in
           con.execute("select prefix, asmid, cls, species from gen").fetchall()}
    by = fetch_members(con, a, specs)
    print({k: len(v) for k, v in by.items()}, "member regions before sequence extraction", flush=True)
    con.close()  # no DuckDB threads alive when worker processes are created

    paths = read_manifest(a.manifest)
    per_prefix = defaultdict(list)
    for famkey, regs in by.items():
        for pid, s, e, meta in regs:
            per_prefix[pid.split("_", 1)[0]].append((famkey, pid, s, e, meta))
    tasks, no_file = [], 0
    for prefix, rows in per_prefix.items():
        asm = gen[prefix][0]
        p = paths.get(asm)
        if not p or not os.path.exists(p):
            no_file += 1
            continue
        tasks.append((prefix, p, rows))
    print(f"{len(tasks)} genome files to read; {no_file} genomes without protein file", flush=True)

    fas, mem, cnt = {}, {}, defaultdict(lambda: defaultdict(int))
    genomes = defaultdict(lambda: defaultdict(set))
    miss_total = 0
    try:
        with mp.get_context("forkserver").Pool(a.threads) as pool:
            for prefix, out, miss in pool.imap_unordered(work, tasks, chunksize=4):
                if out is None:
                    print("READ FAIL", prefix, miss, flush=True)
                    continue
                miss_total += miss
                asm, cls, sp = gen[prefix]
                for famkey, pid, s, e, meta, seq, plen in out:
                    if famkey not in fas:
                        fas[famkey] = open(os.path.join(a.outdir, famkey + ".faa"), "w")
                        mem[famkey] = open(os.path.join(a.outdir, famkey + ".members.tsv"), "w")
                        mem[famkey].write("\t".join(HEADER) + "\n")
                    rid = f"{pid}|{s}-{e}"
                    fas[famkey].write(f">{rid}\n{seq}\n")
                    g = lambda k: ("NA" if meta.get(k) is None else str(meta[k]))  # noqa: E731
                    mem[famkey].write("\t".join(map(str, [
                        rid, pid, asm, cls, sp, s, e, e - s + 1, plen, g("score"), g("evalue"),
                        g("model_cov"), g("sub"), g("pident"), g("efam"), g("ec"), g("substrate")])) + "\n")
                    cnt[famkey]["regions"] += 1
                    cnt[famkey]["len"] += e - s + 1
                    genomes[famkey][cls].add(asm)
                    genomes[famkey]["prot"].add(pid)
    finally:
        for h in list(fas.values()) + list(mem.values()):
            h.close()

    with open(os.path.join(a.outdir, "regions_summary.tsv"), "w") as fh:
        fh.write("family\tregions\tproteins\tgenomes_baseline\tgenomes_new\tmean_len\n")
        for k in sorted(cnt):
            fh.write(f"{k}\t{cnt[k]['regions']}\t{len(genomes[k]['prot'])}\t{len(genomes[k]['baseline'])}\t"
                     f"{len(genomes[k]['new'])}\t{cnt[k]['len'] / cnt[k]['regions']:.1f}\n")
    print(f"done; ids not found in protein files: {miss_total}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
