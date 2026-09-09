from lib_accumulation import run_permutation, run_all_permutations

ASMID_CLUSTERS = {
    "A": {"c1", "c2"},
    "B": {"c2", "c3"},   # c2 shared with A -> not new when B follows A
    "C": {"c4"},
}
EXTERNAL_STATUS = {"c1": "strong_hit", "c2": "no_hit", "c3": "no_hit", "c4": "weak_hit"}

def test_single_walk_marginal_and_cumulative_counts():
    records = run_permutation(["A", "B", "C"], ASMID_CLUSTERS, EXTERNAL_STATUS)
    by_asmid = {r["asmid"]: r for r in records}
    assert by_asmid["A"]["new_internal_count"] == 2       # c1, c2 both new
    assert by_asmid["A"]["internal_cumulative"] == 2
    assert by_asmid["A"]["new_external_count"] == 1        # only c2 is no_hit
    assert by_asmid["B"]["new_internal_count"] == 1         # only c3 is new (c2 seen)
    assert by_asmid["B"]["internal_cumulative"] == 3
    assert by_asmid["B"]["new_external_count"] == 1          # c3 is no_hit and new
    assert by_asmid["C"]["new_internal_count"] == 1
    assert by_asmid["C"]["internal_cumulative"] == 4
    assert by_asmid["C"]["new_external_count"] == 0           # c4 is weak_hit, not no_hit

def test_run_all_permutations_is_reproducible_with_seed():
    walks_1 = run_all_permutations(["A", "B", "C"], ASMID_CLUSTERS, EXTERNAL_STATUS,
                                    n_permutations=5, seed=42)
    walks_2 = run_all_permutations(["A", "B", "C"], ASMID_CLUSTERS, EXTERNAL_STATUS,
                                    n_permutations=5, seed=42)
    orders_1 = [[r["asmid"] for r in w] for w in walks_1]
    orders_2 = [[r["asmid"] for r in w] for w in walks_2]
    assert orders_1 == orders_2

def test_every_permutation_visits_every_genome_exactly_once():
    walks = run_all_permutations(["A", "B", "C"], ASMID_CLUSTERS, EXTERNAL_STATUS,
                                  n_permutations=10, seed=1)
    for w in walks:
        assert sorted(r["asmid"] for r in w) == ["A", "B", "C"]
