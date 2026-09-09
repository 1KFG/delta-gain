import pandas as pd
from plot_accumulation_plotly import build_figures

def test_build_figures_returns_all_three():
    summary_df = pd.DataFrame({
        "N": [1, 2, 3],
        "internal_mean": [2.0, 3.0, 4.0], "internal_p2_5": [2.0, 3.0, 4.0], "internal_p97_5": [2.0, 3.0, 4.0],
        "external_mean": [1.0, 2.0, 2.0], "external_p2_5": [1.0, 2.0, 2.0], "external_p97_5": [1.0, 2.0, 2.0],
    })
    fit_df = pd.DataFrame({
        "series": ["internal", "external"],
        "alpha_mean": [0.6, 0.4], "alpha_ci_low": [0.4, 0.2], "alpha_ci_high": [0.8, 0.6],
        "kappa_mean": [3.0, 2.0], "r_squared_mean": [0.9, 0.85],
    })
    clade_df = pd.DataFrame({
        "clade_label": ["Pleosporales", "Agaricomycetes"],
        "mean_marginal_internal": [1.5, 1.0],
    })
    figures = build_figures(summary_df, fit_df, clade_df, "test-v1-abc1234")
    assert set(figures) == {"internal_curve", "external_curve", "clade_contribution"}
    for fig in figures.values():
        assert len(fig.data) > 0
        # Design doc: the run's version tag goes in the figure's own title,
        # not just in the output filename.
        assert "test-v1-abc1234" in fig.layout.title.text
    # Design doc "Known caveat": stated on the two accumulation curves,
    # where a saturation claim would be read.
    for name in ("internal_curve", "external_curve"):
        assert "not a calibrated estimate" in figures[name].layout.title.text
    assert "not a calibrated estimate" not in figures["clade_contribution"].layout.title.text
