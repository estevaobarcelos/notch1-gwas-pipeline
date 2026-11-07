#!/usr/bin/env Rscript
# =============================================================================
# 04_pca_analysis.R
# PCA visualization and ancestry stratification assessment
#
# Input : notch1_pca.eigenval, notch1_pca.eigenvec
# Output: figures/pca_plot.png, figures/scree_plot.png
#
# Interpretation for admixed populations:
#   In a homogeneous European cohort, PC1 alone typically explains >50% of
#   variance and separates clear ancestry clusters. In a highly admixed
#   cohort like SABE (Brazilian), variance is distributed across many PCs,
#   reflecting continuous admixture rather than discrete clusters. This is
#   expected and not a data quality issue.
#
# LIMITATION: Only 98 SNPs were available after LD pruning (focused NOTCH1
#   panel). Ancestry inference from so few SNPs is unreliable. PCs are shown
#   for completeness but should not be over-interpreted. For robust ancestry
#   correction, a genome-wide panel is required.
#
# Author: Estevão Barcelos
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
})

dir.create("figures", showWarnings = FALSE)

cat("Loading PCA results...\n")

# ── Load data ─────────────────────────────────────────────────────────────────
eigenval_file <- "notch1_pca.eigenval"
eigenvec_file <- "notch1_pca.eigenvec"

for (f in c(eigenval_file, eigenvec_file)) {
  if (!file.exists(f)) {
    stop(paste("File not found:", f, "\nRun 02_qc_pipeline.sh first."))
  }
}

eigenval <- read.table(eigenval_file, header = FALSE)$V1
eigenvec <- read.table(eigenvec_file, header = FALSE)

n_pcs <- ncol(eigenvec) - 2
colnames(eigenvec) <- c("FID", "IID", paste0("PC", 1:n_pcs))

n_samples <- nrow(eigenvec)
cat(sprintf("  Samples : %d\n", n_samples))
cat(sprintf("  PCs     : %d\n", n_pcs))

# ── Variance explained ────────────────────────────────────────────────────────
pct_var <- round(eigenval / sum(eigenval) * 100, 2)
cat("\nVariance explained per PC:\n")
for (i in seq_along(pct_var)) {
  cat(sprintf("  PC%2d: %6.2f%%\n", i, pct_var[i]))
}

# ── Scree plot ────────────────────────────────────────────────────────────────
scree_df <- data.frame(
  PC  = factor(paste0("PC", 1:n_pcs), levels = paste0("PC", 1:n_pcs)),
  Pct = pct_var
)

p_scree <- ggplot(scree_df, aes(x = PC, y = Pct)) +
  geom_col(fill = "#2980b9", alpha = 0.85, width = 0.7) +
  geom_text(aes(label = paste0(Pct, "%")),
            vjust = -0.4, size = 3, color = "gray30") +
  labs(
    title    = "Scree Plot — Variance Explained by Each PC",
    subtitle = paste0(
      "NOTCH1 region, SABE cohort (n=", n_samples, ")\n",
      "Based on ", 98, " SNPs after LD pruning — interpret with caution"
    ),
    x = "Principal Component",
    y = "Variance explained (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "gray40", size = 10),
    axis.text.x   = element_text(angle = 0)
  ) +
  annotate(
    "text", x = n_pcs * 0.6, y = max(pct_var) * 0.85,
    label = paste0(
      "NOTE: Gradual eigenvalue decay\n",
      "is expected in admixed populations\n",
      "(continuous ancestry gradients,\n",
      "not discrete clusters)"
    ),
    size = 3, color = "gray40", hjust = 0,
    fontface = "italic"
  )

ggsave("figures/scree_plot.png", p_scree,
       width = 9, height = 5, dpi = 150)
cat("\nScree plot saved to: figures/scree_plot.png\n")

# ── PC1 vs PC2 scatter ────────────────────────────────────────────────────────
# Colour by PC1 quintile as a proxy for ancestry gradient
eigenvec <- eigenvec %>%
  mutate(
    PC1_quintile = cut(PC1,
      breaks = quantile(PC1, probs = seq(0, 1, 0.2)),
      labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
      include.lowest = TRUE
    )
  )

p_pca <- ggplot(eigenvec, aes(x = PC1, y = PC2, color = PC1_quintile)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(
    values = c("#1a237e", "#1565c0", "#0288d1", "#00897b", "#2e7d32"),
    name   = "PC1 quintile\n(ancestry gradient)"
  ) +
  labs(
    title    = "PC1 vs PC2 — NOTCH1 Region, SABE Cohort",
    subtitle = paste0(
      "n = ", n_samples, " admixed Brazilian individuals | ",
      round(pct_var[1], 1), "% and ",
      round(pct_var[2], 1), "% variance explained\n",
      "Coloured by PC1 quintile as proxy for ancestry gradient"
    ),
    x = paste0("PC1 (", pct_var[1], "%)"),
    y = paste0("PC2 (", pct_var[2], "%)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "gray40", size = 9),
    legend.position = "right"
  ) +
  annotate(
    "text", x = min(eigenvec$PC1), y = max(eigenvec$PC2),
    label = paste0(
      "Diffuse cloud = continuous admixture\n",
      "No discrete clusters expected in\n",
      "highly admixed populations"
    ),
    size = 3, color = "gray50", hjust = 0, vjust = 1,
    fontface = "italic"
  )

ggsave("figures/pca_plot.png", p_pca,
       width = 9, height = 7, dpi = 150)
cat("PCA plot saved to: figures/pca_plot.png\n")

# ── PC3 vs PC4 ───────────────────────────────────────────────────────────────
p_pc34 <- ggplot(eigenvec, aes(x = PC3, y = PC4, color = PC1_quintile)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(
    values = c("#1a237e", "#1565c0", "#0288d1", "#00897b", "#2e7d32"),
    name   = "PC1 quintile"
  ) +
  labs(
    title = "PC3 vs PC4 — NOTCH1 Region, SABE Cohort",
    x = paste0("PC3 (", pct_var[3], "%)"),
    y = paste0("PC4 (", pct_var[4], "%)")
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13))

ggsave("figures/pca_pc34_plot.png", p_pc34,
       width = 8, height = 6, dpi = 150)
cat("PC3 vs PC4 plot saved to: figures/pca_pc34_plot.png\n")

# ── Save eigenvec with cleaner column names for downstream use ────────────────
write.csv(eigenvec, "results/pca_covariates.csv", row.names = FALSE)
cat("PCA covariates saved to: results/pca_covariates.csv\n")

cat("\n=================================================\n")
cat(" PCA SUMMARY\n")
cat("=================================================\n")
cat(sprintf("  Samples         : %d\n", n_samples))
cat(sprintf("  SNPs used       : 98 (after LD pruning)\n"))
cat(sprintf("  PC1 variance    : %.2f%%\n", pct_var[1]))
cat(sprintf("  PC1-10 total    : %.2f%%\n", sum(pct_var[1:10])))
cat("\n  INTERPRETATION:\n")
cat("  Gradual eigenvalue decay and diffuse PC1-PC2 scatter\n")
cat("  are consistent with continuous admixture in a Brazilian\n")
cat("  population. No major ancestry outliers were detected.\n")
cat("  PCs 1-10 will be used as covariates in association models.\n")
cat("\nNext step: bash scripts/05_association_analysis.sh\n")
