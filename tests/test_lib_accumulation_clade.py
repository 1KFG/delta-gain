from lib_accumulation import run_all_permutations, aggregate_clade_contribution

ASMID_CLUSTERS = {"A": {"c1", "c2"}, "B": {"c2", "c3"}, "C": {"c4"}}
EXTERNAL_STATUS = {"c1": "strong_hit", "c2": "no_hit", "c3": "no_hit", "c4": "weak_hit"}
CLADE_OF = {"A": ("ORDER", "Pleosporales"), "B": ("ORDER", "Pleosporales"),
            "C": ("CLASS", "Agaricomycetes")}

def test_clade_contribution_sums_to_total_clusters():
    walks = run_all_permutations(["A", "B", "C"], ASMID_CLUSTERS, EXTERNAL_STATUS,
                                  n_permutations=200, seed=7)
    rows = aggregate_clade_contribution(walks, CLADE_OF)
    by_label = {r["clade_label"]: r for r in rows}
    assert set(by_label) == {"Pleosporales", "Agaricomycetes"}
    assert by_label["Pleosporales"]["n_genomes"] == 2
    assert by_label["Agaricomycetes"]["n_genomes"] == 1
    # Total internal clusters = 4 (c1,c2,c3,c4); every genome's expected
    # marginal contribution under random order sums exactly to that total
    # (Shapley efficiency property) regardless of insertion order.
    total_internal = sum(r["total_marginal_internal"] for r in rows)
    assert abs(total_internal - 4.0) < 1e-9

def test_mean_marginal_is_total_divided_by_n_genomes():
    walks = run_all_permutations(["A", "B", "C"], ASMID_CLUSTERS, EXTERNAL_STATUS,
                                  n_permutations=50, seed=3)
    rows = aggregate_clade_contribution(walks, CLADE_OF)
    for r in rows:
        assert abs(r["mean_marginal_internal"] - r["total_marginal_internal"] / r["n_genomes"]) < 1e-9
