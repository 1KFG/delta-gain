#!/usr/bin/env python3
"""Per-species pangenome summaries from the genome x cluster incidence matrix.

Genomes are grouped by samples.csv SPECIES. For each species with at least
--min-genomes genomes this reports pangenome size, core / soft-core / shell /
cloud counts, a rarefied pangenome size and core at a common genome count (so
species with different sampling compare), and a power-law fit of the
new-cluster rate.

Definitions (integer arithmetic; n = genomes in the species):
  core       present in all n genomes
  soft_core  present in >= ceil(0.95 n) genomes but not all. Empty by
             construction when n < 20, so it is left blank there.
  cloud      present in <= max(1, floor(0.15 n)) genomes (singletons are
             always cloud)
  shell      everything else
Openness: new clusters per added genome, averaged over --n-permutations random
genome orders, is fitted as kappa * N^-alpha for N = 2..n on the MEAN curve
(N = 1 adds a whole genome and is a different quantity, so it is excluded).
alpha < 1 -> open, alpha >= 1 -> closed. A species whose genomes add nothing
after the first is "closed" (fit_status no_new_after_first), not dropped.
Fewer than 3 usable points -> "unresolved". alpha is not comparable across
species with very different n; compare rarefied_* columns instead.

Species names that are placeholders ("Genus sp.", "cf.", "aff.", "uncultured",
"unclassified") are excluded unless --include-placeholder-species, because
they merge unrelated lineages into one artificially open pangenome.

Clusters are whatever clustering produced the incidence matrix (currently
MMSEQS linclust 95% id / 90% cov). Pass --exclude-asmids to drop redundant
strains before counting.
"""
import argparse
import csv
import re
import sys
import zlib
from multiprocessing import Pool

import numpy as np
import pyarrow.parquet as pq

PLACEHOLDER = re.compile(r"(\bsp\.?(?=\s|$)|\bsp\d|\bcf\.|\baff\.|uncultured|unclassified)", re.I)
FIELDS = ["species", "n_genomes", "pan_size", "core", "soft_core", "shell", "cloud",
          "singleton_clusters", "mean_genome_clusters", "core_per_mean_genome",
          "rarefied_pan_size", "rarefied_pan_size_sd", "rarefied_core",
          "fit_status", "alpha", "alpha_n_points", "alpha_r2", "kappa", "openness"]


def load_pairs(path):
    """Dictionary-encoded read: int32 codes instead of 2 Python strings per row."""
    t = pq.read_table(path, columns=["asmid", "cluster_id"],
                      read_dictionary=["asmid", "cluster_id"]).unify_dictionaries()
    a = t.column("asmid").combine_chunks()
    c = t.column("cluster_id").combine_chunks()
    return (a.indices.to_numpy(zero_copy_only=False).astype(np.int32),
            c.indices.to_numpy(zero_copy_only=False).astype(np.int32),
            a.dictionary.to_pylist(), len(c.dictionary))


def load_species(samples_csv):
    species_of, seen = {}, set()
    with open(samples_csv) as fh:
        for row in csv.DictReader(fh):
            if row["ASMID"] in seen:
                sys.exit(f"duplicate ASMID in {samples_csv}: {row['ASMID']}")
            seen.add(row["ASMID"])
            species_of[row["ASMID"]] = (row.get("SPECIES") or "").strip()
    return species_of


def load_exclude(path):
    if not path:
        return set()
    with open(path) as fh:
        return {line.strip().split("\t")[0] for line in fh if line.strip()}


def group_species(asm_names, species_of, exclude, min_genomes, include_placeholder):
    by, unknown, placeholder = {}, 0, set()
    for code, a in enumerate(asm_names):
        if a in exclude:
            continue
        sp = species_of.get(a, "")
        if not sp:
            unknown += 1
        elif not include_placeholder and PLACEHOLDER.search(sp):
            placeholder.add(sp)
        else:
            by.setdefault(sp, []).append(code)
    return {sp: c for sp, c in by.items() if len(c) >= min_genomes}, unknown, placeholder


def class_counts(counts, n):
    soft_min = -(-95 * n // 100)          # ceil(0.95 n)
    cloud_max = max(1, (15 * n) // 100)   # floor(0.15 n), singletons always cloud
    core = int((counts == n).sum())
    cloud = int((counts <= cloud_max).sum())
    soft = int(((counts >= soft_min) & (counts < n) & (counts > cloud_max)).sum())
    return core, soft, cloud


def fit_mean_curve(mean_new):
    """log-log OLS of the mean new-cluster count at N = 2..n."""
    n = len(mean_new)
    pos = np.arange(2, n + 1, dtype=float)
    y = mean_new[1:]
    ok = y > 0
    if ok.sum() == 0:
        return {"fit_status": "no_new_after_first", "openness": "closed"}
    if ok.sum() < 3:
        return {"fit_status": "too_few_points", "openness": "unresolved",
                "alpha_n_points": int(ok.sum())}
    lx, ly = np.log(pos[ok]), np.log(y[ok])
    slope, icpt = np.polyfit(lx, ly, 1)
    pred = icpt + slope * lx
    ss_tot = ((ly - ly.mean()) ** 2).sum()
    r2 = 1 - ((ly - pred) ** 2).sum() / ss_tot if ss_tot > 0 else 1.0
    alpha = -slope
    return {"fit_status": "ok", "alpha": alpha, "alpha_n_points": int(ok.sum()),
            "alpha_r2": r2, "kappa": float(np.exp(icpt)),
            "openness": "open" if alpha < 1 else "closed"}


def summarize(args):
    name, g_idx, c_raw, n, rarefy_n, n_perm, seed = args
    uc, c_idx = np.unique(c_raw, return_inverse=True)
    m = len(uc)
    key = np.unique(g_idx.astype(np.int64) * m + c_idx)       # drop duplicate pairs
    g_idx = (key // m).astype(np.int32)
    c_idx = (key % m).astype(np.int32)
    order = np.argsort(c_idx, kind="stable")
    g_sorted = g_idx[order]
    starts = np.flatnonzero(np.r_[True, c_idx[order][1:] != c_idx[order][:-1]])
    counts = np.bincount(c_idx, minlength=m)

    rng = np.random.default_rng([seed, zlib.crc32(name.encode())])
    new = np.zeros((n_perm, n))
    rare_core = np.zeros(n_perm)
    for p in range(n_perm):
        rank = rng.permutation(n)                              # rank[g] = position of genome g
        first = np.minimum.reduceat(rank[g_sorted], starts)    # first position of each cluster
        new[p] = np.bincount(first, minlength=n)
        sel = rank[g_idx] < rarefy_n
        rare_core[p] = (np.bincount(c_idx[sel], minlength=m) == rarefy_n).sum()
    cum = new.cumsum(axis=1)

    core, soft, cloud = class_counts(counts, n)
    mean_clusters = len(g_idx) / n
    row = {f: "" for f in FIELDS}
    row.update({
        "species": name, "n_genomes": n, "pan_size": m, "core": core,
        "soft_core": soft if n >= 20 else "", "cloud": cloud,
        "shell": m - core - cloud - (soft if n >= 20 else 0),
        "singleton_clusters": int((counts == 1).sum()),
        "mean_genome_clusters": mean_clusters, "core_per_mean_genome": core / mean_clusters,
        "rarefied_pan_size": float(cum[:, rarefy_n - 1].mean()),
        "rarefied_pan_size_sd": float(cum[:, rarefy_n - 1].std()),
        "rarefied_core": float(rare_core.mean()),
    })
    row.update(fit_mean_curve(new.mean(axis=0)))
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--incidence-matrix", required=True)
    ap.add_argument("--samples-csv", required=True)
    ap.add_argument("--exclude-asmids", default="", help="file of ASMIDs to drop (first column)")
    ap.add_argument("--min-genomes", type=int, default=10)
    ap.add_argument("--rarefy-n", type=int, default=10)
    ap.add_argument("--n-permutations", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--include-placeholder-species", action="store_true")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    if args.rarefy_n > args.min_genomes:
        ap.error("--rarefy-n must be <= --min-genomes")

    a_codes, c_codes, asm_names, _ = load_pairs(args.incidence_matrix)
    groups, unknown, placeholder = group_species(
        asm_names, load_species(args.samples_csv), load_exclude(args.exclude_asmids),
        args.min_genomes, args.include_placeholder_species)
    if not groups:
        sys.exit("no species met --min-genomes; nothing to write")

    sp_names = sorted(groups)
    sp_of_asm = np.full(len(asm_names), -1, dtype=np.int32)
    local_g = np.zeros(len(asm_names), dtype=np.int32)
    for i, sp in enumerate(sp_names):
        codes = np.array(groups[sp], dtype=np.int32)
        sp_of_asm[codes] = i
        local_g[codes] = np.arange(len(codes), dtype=np.int32)
    s = sp_of_asm[a_codes]
    keep = s >= 0
    a_codes, c_codes, s = a_codes[keep], c_codes[keep], s[keep]
    order = np.argsort(s, kind="stable")
    a_codes, c_codes, s = a_codes[order], c_codes[order], s[order]
    bounds = np.searchsorted(s, np.arange(len(sp_names) + 1))
    jobs = [(sp, local_g[a_codes[bounds[i]:bounds[i + 1]]], c_codes[bounds[i]:bounds[i + 1]],
             len(groups[sp]), args.rarefy_n, args.n_permutations, args.seed)
            for i, sp in enumerate(sp_names)]
    if args.threads > 1:
        with Pool(args.threads) as pool:
            rows = pool.map(summarize, jobs, chunksize=1)
    else:
        rows = [summarize(j) for j in jobs]

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"species written: {len(rows)}; genomes with no species: {unknown}; "
          f"placeholder species skipped: {len(placeholder)}", file=sys.stderr)


if __name__ == "__main__":
    main()
