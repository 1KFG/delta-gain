#!/usr/bin/env Rscript
# Static, publication-styled figures for the protein-content accumulation
# curve (design doc: Outputs/plots). Primary plotting implementation --
# covers the full 10-figure set; Plotly (bin/plot_accumulation_plotly.py)
# duplicates only figures 1-3 for interactive exploration.

suppressMessages({
  library(ggplot2)
  library(arrow)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
opt <- list(
  summary_tsv = args[which(args == "--summary") + 1],
  fit_tsv     = args[which(args == "--fit") + 1],
  clade_tsv   = args[which(args == "--clade") + 1],
  version_tag = args[which(args == "--version-tag") + 1],
  outdir      = args[which(args == "--outdir") + 1]
)

summary_df <- read_tsv_arrow(opt$summary_tsv)
fit_df     <- read_tsv_arrow(opt$fit_tsv)
clade_df   <- read_tsv_arrow(opt$clade_tsv)

caption <- paste("version:", opt$version_tag)

# Figure 1/2: internal + external accumulation curves with power-law overlay.
plot_accumulation_curve <- function(df, mean_col, lo_col, hi_col, fit_row, title) {
  alpha <- fit_row$alpha_mean
  kappa <- fit_row$kappa_mean
  # Integrated rate-form model as a predicted cumulative trajectory overlay.
  predicted_cumulative <- cumsum(kappa * df$N ^ (-alpha))
  ggplot(df, aes(x = N)) +
    geom_ribbon(aes(ymin = .data[[lo_col]], ymax = .data[[hi_col]]), alpha = 0.2) +
    geom_line(aes(y = .data[[mean_col]]), linewidth = 0.8) +
    geom_line(aes(y = predicted_cumulative), linetype = "dashed", color = "firebrick") +
    labs(title = title,
         subtitle = sprintf("alpha=%.2f [%.2f, %.2f], kappa=%.1f", alpha, fit_row$alpha_ci_low, fit_row$alpha_ci_high, kappa),
         x = "Genomes added (N)", y = "Cumulative distinct clusters", caption = caption) +
    theme_minimal()
}

fig1 <- plot_accumulation_curve(summary_df, "internal_mean", "internal_p2_5", "internal_p97_5",
                                 fit_df[fit_df$series == "internal", ], "Internal accumulation curve")
fig2 <- plot_accumulation_curve(summary_df, "external_mean", "external_p2_5", "external_p97_5",
                                 fit_df[fit_df$series == "external", ], "External-novelty accumulation curve")

# Figure 3: per-clade marginal contribution, sorted descending, roll-up
# labels shown verbatim from clade_contribution.tsv (already carries the
# "ORDER (other families)" convention from Task 5/7).
clade_df_sorted <- clade_df[order(-clade_df$mean_marginal_internal), ]
fig3 <- ggplot(clade_df_sorted, aes(x = reorder(clade_label, mean_marginal_internal), y = mean_marginal_internal)) +
  geom_col() + coord_flip() +
  labs(title = "Per-clade mean marginal contribution (internal)",
       x = NULL, y = "Mean new clusters per genome added", caption = caption) +
  theme_minimal()

pdf(file.path(opt$outdir, sprintf("accumulation_figures_%s.pdf", opt$version_tag)), width = 8, height = 6)
print(fig1); print(fig2); print(fig3)
dev.off()
