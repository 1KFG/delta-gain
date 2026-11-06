"""Deterministic checks for bin/family_null_order.py (leave-one-order-out null)."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
import family_null_order as fo  # noqa: E402


def scenario():
    # o1: baseline b1..b3 ; o2: baseline b4..b6 ; o3: new-only n1,n2 ; o1 also has new n3
    g_rows = [(f"b{i + 1}", "baseline", "o1" if i < 3 else "o2", f"s{i}") for i in range(6)]
    g_rows += [("n1", "new", "o3", "t1"), ("n2", "new", "o3", "t2"), ("n3", "new", "o1", "t3")]
    g_r, cl = [], []
    for i in range(6):  # each baseline genome: shared cluster 0 + an order-specific cluster (1 for o1, 2 for o2)
        g_r += [i, i]
        cl += [0, 1 if i < 3 else 2]
    g_r += [6, 6, 7, 7]  # new-only order: both regions in private clusters
    cl += [3, 4, 5, 6]
    g_r += [8, 8]  # new genome of o1: shared + o1-specific -> o1-specific is novel when o1 is held out
    cl += [0, 1]
    return g_rows, {"fam": (["c90"], np.array(g_r), {"c90": np.array(cl)})}


def test_order_null():
    g_rows, fam = scenario()
    rows = fo.run_order_null(g_rows, fam, K=3, n_draws=10, seed=1, min_base=3, min_new=2)
    d = {(r[2], r[3], r[4]): dict(zip(fo.COLS, r)) for r in rows}
    for v in ("genome", "species"):
        # o1 held out: cluster 0 covered by o2, cluster 1 absent -> half novel
        assert d[(v, "o1", "base_out")]["stat_mean"] == 0.5
        assert d[(v, "o2", "base_out")]["stat_mean"] == 0.5
        # new-only order: everything novel (nb == 0, reference is all baseline)
        r = d[(v, "NEW_ONLY_POOL", "new_only_pool")]
        assert r["stat_mean"] == 1.0 and r["n_base"] == 0 and r["n_new"] == 2 and r["k"] == 2
        # o3 has 2 new genomes (>= min_new) so it is also reported on its own, reference = whole baseline
        assert d[(v, "o3", "new")]["stat_mean"] == 1.0
        # o1 has only one new genome -> below min_new, not reported
        assert (v, "o1", "new") not in d
