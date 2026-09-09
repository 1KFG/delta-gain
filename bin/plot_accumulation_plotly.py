#!/usr/bin/env python3
"""Interactive figures 1-3 (design doc: Outputs/plots) -- only the three
figures worth exploring interactively; the diagnostic/validation figures
stay ggplot-only (bin/plot_accumulation_ggplot.R)."""
import argparse

import pandas as pd
import plotly.graph_objects as go


def _curve_figure(df, mean_col, lo_col, hi_col, title):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=df["N"], y=df[hi_col], mode="lines",
                              line=dict(width=0), showlegend=False))
    fig.add_trace(go.Scatter(x=df["N"], y=df[lo_col], mode="lines", fill="tonexty",
                              line=dict(width=0), name="95% CI"))
    fig.add_trace(go.Scatter(x=df["N"], y=df[mean_col], mode="lines", name="mean"))
    fig.update_layout(title=title, xaxis_title="Genomes added (N)",
                       yaxis_title="Cumulative distinct clusters")
    return fig


def build_figures(summary_df, fit_df, clade_df):
    internal_curve = _curve_figure(summary_df, "internal_mean", "internal_p2_5",
                                    "internal_p97_5", "Internal accumulation curve")
    external_curve = _curve_figure(summary_df, "external_mean", "external_p2_5",
                                    "external_p97_5", "External-novelty accumulation curve")
    clade_sorted = clade_df.sort_values("mean_marginal_internal", ascending=True)
    clade_contribution = go.Figure(go.Bar(
        x=clade_sorted["mean_marginal_internal"], y=clade_sorted["clade_label"], orientation="h"
    ))
    clade_contribution.update_layout(title="Per-clade mean marginal contribution",
                                      xaxis_title="Mean new clusters per genome added")
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
    figures = build_figures(summary_df, fit_df, clade_df)
    for name, fig in figures.items():
        fig.write_html(f"{args.outdir}/{name}_{args.version_tag}.html")
