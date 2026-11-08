#!/usr/bin/env python3
"""Pick no-hit profile-cluster representatives for the structure pilot (large, multi-member, 100-400 aa)."""
import random
import sys
from collections import Counter

R = "/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/results/nohit_profile"
OUT = sys.argv[1]
N = int(sys.argv[2])
MODE = sys.argv[3] if len(sys.argv) > 3 else "nohit"
size = Counter()
with open(f"{R}/profile_clusters.tsv") as fh:
    for l in fh:
        rep, mem = l.rstrip("\n").split("\t")[:2]
        size[rep] += 1
hit = set()
with open(f"{R}/profile_hits.tsv") as fh:
    for l in fh:
        hit.add(l.split("\t", 1)[0])
cand = {r: n for r, n in size.items() if n >= 20 and ((r not in hit) if MODE == "nohit" else (r in hit))}
print("profile clusters", len(size), "no-hit with >=20 members", len(cand), flush=True)
rnd = random.Random(7)
want = set(rnd.sample(sorted(cand), min(len(cand), N * 4)))
out = []
name = None
seq = []
def flush():
    if name in want:
        s = "".join(seq)
        if 100 <= len(s) <= 400:
            out.append((name, s))
with open(f"{R}/members.faa") as fh:
    for l in fh:
        if l.startswith(">"):
            flush()
            name = l[1:].strip().split()[0]
            seq = []
        else:
            seq.append(l.strip())
flush()
rnd.shuffle(out)
out = out[:N]
with open(OUT, "w") as fh:
    for n, s in out:
        fh.write(f">{n} members={cand[n]}\n{s}\n")
print("selected", len(out))
