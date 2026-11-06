"""Deterministic checks for bin/family_null.py (baseline-vs-baseline null)."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
import family_null as fn  # noqa: E402


def scenario(shared_species=False):
    # order o1: baseline b1..b3, new n1..n3
    sp = ["A", "A", "B", "A", "C", "C"] if shared_species else ["s1", "s2", "s3", "s4", "s5", "s6"]
    g_rows = [(f"b{i + 1}", "baseline", "o1", sp[i]) for i in range(3)] + \
             [(f"n{i + 1}", "new", "o1", sp[3 + i]) for i in range(3)]
    # clusters: X=0 shared by all baseline genomes, U_i private to baseline genome i, V_* private to new genomes
    g_r, cl = [], []
    for i in range(3):
        g_r += [i, i]
        cl += [0, 1 + i]
    for j in range(3):
        g_r += [3 + j, 3 + j]
        cl += [4 + 2 * j, 5 + 2 * j]
    return g_rows, {"fam": (["c90"], np.array(g_r), {"c90": np.array(cl)})}


def test_new_beats_baseline_null():
    g_rows, fam = scenario()
    rows = fn.run_null(g_rows, fam, [0.34], 20, 1)
    d = {r[2]: dict(zip(fn.COLS, r)) for r in rows}
    for variant in ("genome", "species"):
        r = d[variant]
        assert (r["genomes_H"], r["genomes_N"], r["orders"]) == (1.0, 1.0, 1)
        # held-out baseline genome: X is covered by the others, its private cluster is not -> 0.5
        assert r["stat_null"] == 0.5 and r["null_lo"] == r["null_hi"] == 0.5
        # new genome: both regions fall in clusters absent from the reference -> 1.0
        assert r["stat_new"] == 1.0
        assert r["diff_mean"] == 0.5 and r["diff_lo"] == r["diff_hi"] == 0.5 and r["p_diff_le0"] == 0.0
        assert r["known_new"] is None  # no species shared with the reference


def test_known_species_stratum():
    g_rows, fam = scenario(shared_species=True)  # n1 shares species A with the baseline; C only new
    rows = fn.run_null(g_rows, fam, [0.34], 60, 1)
    r = {x[2]: dict(zip(fn.COLS, x)) for x in rows}["genome"]
    # when a drawn new genome is n1 (species A, present in the reference unless both A baseline genomes
    # are held out - impossible with 1 held out), its regions count in the known stratum and are all novel
    assert r["known_new"] in (None, 1.0)


def test_no_difference_when_new_equals_baseline():
    g_rows, _ = scenario()
    g_r, cl = [], []
    for i in range(3):  # new genomes mirror the baseline pattern: shared X + a private cluster
        g_r += [i, i, 3 + i, 3 + i]
        cl += [0, 1 + i, 0, 4 + i]
    fam = {"fam": (["c90"], np.array(g_r), {"c90": np.array(cl)})}
    r = {x[2]: dict(zip(fn.COLS, x)) for x in fn.run_null(g_rows, fam, [0.34], 20, 1)}["genome"]
    assert r["stat_new"] == r["stat_null"] == 0.5 and r["diff_mean"] == 0.0 and r["p_diff_le0"] == 1.0


def test_species_collapse_one_genome_per_species():
    rng = np.random.default_rng(0)
    g_sp = np.array([0, 0, 0, 1, 1, 2])
    idx = np.arange(6)
    for _ in range(20):
        pick = fn.one_per_species(rng, idx, g_sp)
        assert sorted(g_sp[pick]) == [0, 1, 2]
