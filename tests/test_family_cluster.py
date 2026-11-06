"""Hand-checked novelty arithmetic for bin/family_cluster.py (no MMseqs2 needed)."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
import family_cluster as fc  # noqa: E402


def test_compose():
    chains = fc.compose({"a": "r0", "b": "r0", "c": "r2"}, [{"r0": "R", "r2": "R"}, {"R": "S"}])
    assert chains == {"a": ["r0", "R", "S"], "b": ["r0", "R", "S"], "c": ["r2", "R", "S"]}
    # a rep absent from a later step stays itself
    assert fc.compose({"x": "q"}, [{}]) == {"x": ["q", "q"]}


def test_summarize(tmp_path):
    mem = tmp_path / "m.tsv"
    mem.write_text("region_id\tclass\tspecies\n"
                   "r1\tbaseline\tspA\nr2\tbaseline\tspA\nr3\tnew\tspA\n"
                   "r4\tnew\tspB\nr5\tnew\tspB\nr6\tnew\tspC\n")
    asg = tmp_path / "a.tsv"
    asg.write_text("region_id\tc90\tc70\n"
                   "r1\tX\tX\nr2\tX\tX\nr3\tX\tX\nr4\tY\tY\nr5\tY\tY\nr6\tZ\tY\n")
    rows = fc.summarize("fam", str(mem), str(asg), ["c90", "c70"], {"spA": True, "spB": False, "spC": False})
    d = [dict(zip(fc.SUMMARY_COLS, r)) for r in rows]
    t90, t70 = d
    assert (t90["regions_base"], t90["regions_new"], t90["clusters"]) == (2, 4, 3)
    assert (t90["clusters_with_base"], t90["clusters_new_only"], t90["clusters_base_only"]) == (1, 2, 0)
    assert t90["new_regions_in_new_only"] == 3 and t90["frac_new_novel"] == 0.75
    assert (t90["new_known_species"], t90["novel_known_species"], t90["frac_novel_known_species"]) == (1, 0, 0.0)
    assert (t90["new_new_taxon"], t90["novel_new_taxon"], t90["frac_novel_new_taxon"]) == (3, 3, 1.0)
    assert (t90["pairs_new_species"], t90["pairs_novel"]) == (3, 2)
    assert (t90["pairs_known_species"], t90["pairs_novel_known_species"]) == (1, 0)
    assert (t70["clusters"], t70["clusters_new_only"], t70["new_regions_in_new_only"]) == (2, 1, 3)
    assert (t70["pairs_new_species"], t70["pairs_novel"]) == (3, 2)
