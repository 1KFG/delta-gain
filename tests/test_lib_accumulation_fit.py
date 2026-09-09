import random
from lib_accumulation import fit_power_law_per_permutation

def _synthetic_walk(n, kappa, alpha, noise_seed):
    """Build a walk whose new_internal_count(N) exactly follows
    kappa * N^-alpha (rounded, with tiny multiplicative noise), to check the
    fit recovers alpha within a reasonable tolerance."""
    rng = random.Random(noise_seed)
    records = []
    for position in range(1, n + 1):
        expected = kappa * position ** (-alpha)
        noisy = max(0, round(expected * rng.uniform(0.9, 1.1)))
        records.append({"position": position, "asmid": f"g{position}",
                         "new_internal_count": noisy, "new_external_count": noisy})
    return records

def test_recovers_known_alpha_within_tolerance():
    walks = [_synthetic_walk(n=200, kappa=500, alpha=0.6, noise_seed=s) for s in range(30)]
    result = fit_power_law_per_permutation(walks, count_key="new_internal_count")
    assert abs(result["alpha_mean"] - 0.6) < 0.15
    assert result["alpha_ci_low"] < result["alpha_mean"] < result["alpha_ci_high"]
    assert len(result["alphas"]) == 30

def test_closed_pan_proteome_gives_alpha_above_one():
    walks = [_synthetic_walk(n=200, kappa=500, alpha=1.5, noise_seed=s) for s in range(30)]
    result = fit_power_law_per_permutation(walks, count_key="new_internal_count")
    assert result["alpha_mean"] > 1.0

def test_zero_count_points_are_excluded_not_fatal():
    walks = [_synthetic_walk(n=500, kappa=50, alpha=1.2, noise_seed=s) for s in range(10)]
    # alpha=1.2 over N=500 guarantees some exact-zero new_internal_count points
    # at large N -- must not raise (log(0)) and must still return a result.
    result = fit_power_law_per_permutation(walks, count_key="new_internal_count")
    assert result["alpha_mean"] > 0
