#!/usr/bin/env python3
"""Streaming, bounded-memory version of bin_novelty_hits.py.

Memory is O(queries in the largest chunk), not O(all queries).

INPUT: --manifest TSV, one row per chunk IN DESIRED OUTPUT ORDER. Header:
    chunk <TAB> queries <TAB> <target1> <TAB> <target2> ...
  'queries' is the chunk FASTA that was given to DIAMOND; each target column
  is that chunk's hit TSV (outfmt 6 + stitle) against that target.

INVARIANT relied on: every query id belongs to exactly one chunk, and hits
for a query appear only in that chunk's hit file for each target. It is
CHECKED, not assumed: every qseqid read from a hit file must be an id of the
chunk FASTA named on the same manifest row (dict lookup, per-chunk set).
Chunk FASTAs have unique ids (checked) so a query can't be in two chunks'
hit files without raising. Output order = manifest order x chunk FASTA order;
--queries (optional, the full FASTA) is streamed in lockstep to prove that
equals the global FASTA order and that no query was dropped or added.

Per query and target: each hit row is classified; the best row wins by
(class strong>weak>no_hit, lower evalue, higher bitscore, first seen).
"""
import argparse
import math
import os
import sys
import tempfile

NAMES = ["no_hit", "weak_hit", "strong_hit"]


def die(msg):
    raise SystemExit(f"bin_novelty_hits: ERROR: {msg}")


def iter_fasta_ids(path):
    with open(path, encoding="utf-8", errors="replace") as fh:
        for n, line in enumerate(fh, 1):
            if line.startswith(">"):
                parts = line[1:].split(None, 1)
                if not parts:
                    die(f"{path}:{n}: empty FASTA id")
                yield parts[0]


def read_manifest(path):
    with open(path) as fh:
        rows = [l.rstrip("\r\n").split("\t") for l in fh
                if l.strip() and not l.startswith("#")]
    if not rows:
        die(f"{path}: empty manifest")
    header = rows[0]
    if header[:2] != ["chunk", "queries"] or len(header) < 3:
        die("manifest header must be: chunk, queries, <target>...")
    targets = header[2:]
    if len(set(targets)) != len(targets):
        die("duplicate target column names in manifest")
    seen, out = set(), []
    for r in rows[1:]:
        if len(r) != len(header):
            die(f"manifest row has {len(r)} columns, expected {len(header)}: {r[0]}")
        if r[0] in seen:
            die(f"duplicate chunk id in manifest: {r[0]}")
        seen.add(r[0])
        out.append((r[0], r[1], r[2:]))
    if not out:
        die("manifest has no chunks")
    return targets, out


def classify(ev, qcov, scov, ev_strong, ev_weak, cov_strong):
    if ev <= ev_strong and qcov >= cov_strong and scov >= cov_strong:
        return 2
    return 1 if ev <= ev_weak else 0


def bin_file(path, idx, thr, mode):
    """Return ({pos: (key, ev_str, qcov_str, scov_str)}, n_rows)."""
    best, nrows = {}, 0
    with open(path, encoding="utf-8", errors="replace") as fh:
        for n, line in enumerate(fh, 1):
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            f = line.split("\t", 14)
            where = f"{path}:{n}"
            if len(f) < 14:
                die(f"{where}: {len(f)} fields, need >=14 (outfmt 6 + stitle); truncated file?")
            pos = idx.get(f[0])
            if pos is None:
                die(f"{where}: qseqid {f[0]!r} is not in this chunk's query FASTA "
                    "(invariant violated: wrong chunk/FASTA pairing or query in two chunks)")
            try:
                ev, bits, qc, sc = float(f[10]), float(f[11]), float(f[12]), float(f[13])
            except ValueError:
                die(f"{where}: non-numeric evalue/bitscore/qcov/scov: {f[10:14]}")
            if not (all(map(math.isfinite, (ev, bits, qc, sc))) and ev >= 0
                    and 0 <= qc <= 100 and 0 <= sc <= 100):
                die(f"{where}: out-of-range evalue/coverage: {f[10:14]}")
            nrows += 1
            cls = classify(ev, qc, sc, *thr)
            # best-evalue (default): the lowest-E-value row alone decides the
            # class, first row wins ties -- identical to the original script.
            # best-class: every row is classified and the best class wins
            # (ties: lower E-value, then higher bitscore, then first seen).
            key = (cls, -ev, bits) if mode == "best-class" else (-ev,)
            cur = best.get(pos)
            if cur is None or key > cur[0]:
                best[pos] = (key, cls, f[10], f[12], f[13])   # raw DIAMOND strings
    return best, nrows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--queries", default="", help="full query FASTA; if given, verify chunk order == FASTA order")
    ap.add_argument("--evalue-strong", type=float, default=1e-5)
    ap.add_argument("--evalue-weak", type=float, default=1e-3)
    ap.add_argument("--cov-strong", type=float, default=70.0)
    ap.add_argument("--classify", choices=["best-evalue", "best-class"], default="best-evalue",
                    help="best-evalue: classify only the lowest-E-value row per query/target "
                         "(original behaviour). best-class: classify every row, keep the best class.")
    ap.add_argument("--fail-on-empty", action="store_true",
                    help="error if any (chunk,target) hit file has zero rows")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if a.evalue_weak < a.evalue_strong:
        die("--evalue-weak must be >= --evalue-strong")
    thr = (a.evalue_strong, a.evalue_weak, a.cov_strong)

    targets, chunks = read_manifest(a.manifest)
    glob = iter_fasta_ids(a.queries) if a.queries else None
    header = ["query_id"]
    for t in targets:
        header += [f"{t}_status", f"{t}_evalue", f"{t}_qcov", f"{t}_scov"]
    header.append("overall_status")

    counts = [0, 0, 0]
    total = 0
    empties = []
    d = os.path.dirname(os.path.abspath(a.out))
    fd, tmp = tempfile.mkstemp(dir=d, prefix=".novelty.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as out:
            out.write("\t".join(header) + "\n")
            for cid, fasta, paths in chunks:
                ids = list(iter_fasta_ids(fasta))
                idx = {q: i for i, q in enumerate(ids)}
                if len(idx) != len(ids):
                    dup = next(q for i, q in enumerate(ids) if idx[q] != i)
                    die(f"chunk {cid}: duplicate query id {dup!r} in {fasta}")
                if glob is not None:
                    for q in ids:
                        g = next(glob, None)
                        if g != q:
                            die(f"chunk {cid}: query {q!r} != global FASTA record {g!r} "
                                "(chunk order/content differs from --queries)")
                per_t = []
                for t, p in zip(targets, paths):
                    best, nrows = bin_file(p, idx, thr, a.classify)
                    if nrows == 0:
                        empties.append((cid, t))
                        if a.fail_on_empty:
                            die(f"chunk {cid} target {t}: hit file {p} has no rows")
                    per_t.append(best)
                for i, q in enumerate(ids):
                    row, overall = [q], 0
                    for best in per_t:
                        h = best.get(i)
                        if h is None:
                            row += ["no_hit", "", "", ""]
                        else:
                            r = h[1]
                            row += [NAMES[r], h[2], h[3], h[4]]
                            overall = max(overall, r)
                    row.append(NAMES[overall])
                    counts[overall] += 1
                    out.write("\t".join(row) + "\n")
                total += len(ids)
            if glob is not None and next(glob, None) is not None:
                die("--queries has records beyond the last chunk")
        os.replace(tmp, a.out)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise
    for cid, t in empties:
        print(f"WARNING: chunk {cid} target {t}: zero hit rows", file=sys.stderr)
    print(f"bin_novelty_hits: {total} queries in {len(chunks)} chunks -> "
          f"strong_hit={counts[2]} weak_hit={counts[1]} no_hit={counts[0]} "
          f"(targets: {', '.join(targets)})", file=sys.stderr)


if __name__ == "__main__":
    main()
