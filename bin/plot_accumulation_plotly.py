#!/usr/bin/env python3
"""Interactive figures 1-3 (design doc: Outputs/plots) -- only the three
figures worth exploring interactively; the diagnostic/validation figures
stay ggplot-only (bin/plot_accumulation_ggplot.R)."""
import argparse

import pandas as pd
import plotly.graph_objects as go


# Design doc, "Known caveat": must be stated explicitly wherever any result
# from this analysis is reported. Shown on figures 1 and 2 specifically --
# the accumulation curves, where a saturation claim would be read.
CAVEAT = ("Describes discovery within the annotated subset only, "
          "not a calibrated estimate of total fungal pan-proteome size.")


def _curve_figure(df, mean_col, lo_col, hi_col, title, version_tag):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=df["N"], y=df[hi_col], mode="lines",
                              line=dict(width=0), showlegend=False))
    fig.add_trace(go.Scatter(x=df["N"], y=df[lo_col], mode="lines", fill="tonexty",
                              line=dict(width=0), name="95% CI"))
    fig.add_trace(go.Scatter(x=df["N"], y=df[mean_col], mode="lines", name="mean"))
    fig.update_layout(
        # Design doc, Outputs/plots: the version tag goes in the figure's own
        # title, not only in the output filename, so a figure pulled out of
        # context still identifies which run produced it.
        title=f"{title} (version: {version_tag})<br><sub>{CAVEAT}</sub>",
        xaxis_title="Genomes added (N)",
        yaxis_title="Cumulative distinct clusters",
    )
    return fig


def build_figures(summary_df, fit_df, clade_df, version_tag):
    internal_curve = _curve_figure(summary_df, "internal_mean", "internal_p2_5",
                                    "internal_p97_5", "Internal accumulation curve",
                                    version_tag)
    external_curve = _curve_figure(summary_df, "external_mean", "external_p2_5",
                                    "external_p97_5", "External-novelty accumulation curve",
                                    version_tag)
    clade_sorted = clade_df.sort_values("mean_marginal_internal", ascending=True)
    clade_contribution = go.Figure(go.Bar(
        x=clade_sorted["mean_marginal_internal"], y=clade_sorted["clade_label"], orientation="h"
    ))
    clade_contribution.update_layout(
        title=f"Per-clade mean marginal contribution (version: {version_tag})",
        xaxis_title="Mean new clusters per genome added",
    )
    return {"internal_curve": internal_curve, "external_curve": external_curve,
            "clade_contribution": clade_contribution}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", required=True)
    ap.add_argument("--fit", required=True)
    ap.add_argument("--clade", required=True)
    ap.add_argument("--version-tag", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    summary_df = pd.read_csv(args.summary, sep="\t")
    fit_df = pd.read_csv(args.fit, sep="\t")
    clade_df = pd.read_csv(args.clade, sep="\t")
    figures = build_figures(summary_df, fit_df, clade_df, args.version_tag)
    for name, fig in figures.items():
        fig.write_html(f"{args.outdir}/{name}_{args.version_tag}.html")
