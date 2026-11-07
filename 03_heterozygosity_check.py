#!/usr/bin/env python3
"""
03_heterozygosity_check.py
==========================
Detect and remove heterozygosity outliers from PLINK .het output.

Methodology:
    - Compute observed heterozygosity: (N_NM - O_HOM) / N_NM
    - Calculate mean and SD across all samples
    - Flag samples outside ±3 SD as outliers
    - Generate a diagnostic plot
    - Write outlier list for PLINK --remove

Why heterozygosity matters:
    - Very LOW het (→ 0): possible sample swap, inbreeding, or complete
      genotyping failure for that individual
    - Very HIGH het (> mean + 3 SD): possible sample contamination
      (two DNA sources mixed), or mapping artefacts
    - In admixed populations, the MEAN heterozygosity is higher than in
      homogeneous populations — this is expected and not a QC failure.
      The ±3 SD criterion accounts for this automatically.

Input : notch1_het.het  (from PLINK --het)
Output: outliers_het.txt, figures/heterozygosity_plot.png

Author: Estevão Barcelos
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import sys
import os

# ── Configuration ─────────────────────────────────────────────────────────────
HET_FILE   = "notch1_het.het"
OUT_LIST   = "outliers_het.txt"
OUT_FIGURE = "figures/heterozygosity_plot.png"
SD_THRESH  = 3.0

os.makedirs("figures", exist_ok=True)

# ── Load data ─────────────────────────────────────────────────────────────────
print(f"Loading {HET_FILE} ...")
try:
    df = pd.read_csv(HET_FILE, sep=r"\s+")
except FileNotFoundError:
    print(f"ERROR: {HET_FILE} not found. Run 02_qc_pipeline.sh first.")
    sys.exit(1)

# PLINK .het columns: FID IID O(HOM) E(HOM) N(NM) F
# Observed heterozygosity = (N_NM - O_HOM) / N_NM
df["HET"] = (df["N(NM)"] - df["O(HOM)"]) / df["N(NM)"]

# ── Statistics ────────────────────────────────────────────────────────────────
mean_het = df["HET"].mean()
sd_het   = df["HET"].std()
lo_thresh = mean_het - SD_THRESH * sd_het
hi_thresh = mean_het + SD_THRESH * sd_het

print(f"\nHeterozygosity statistics:")
print(f"  N samples : {len(df)}")
print(f"  Mean      : {mean_het:.4f}")
print(f"  SD        : {sd_het:.4f}")
print(f"  Lower {SD_THRESH} SD : {lo_thresh:.4f}")
print(f"  Upper {SD_THRESH} SD : {hi_thresh:.4f}")

# ── Flag outliers ─────────────────────────────────────────────────────────────
df["OUTLIER"] = (df["HET"] < lo_thresh) | (df["HET"] > hi_thresh)
n_outliers = df["OUTLIER"].sum()

print(f"\nOutliers (outside ±{SD_THRESH} SD): {n_outliers}")
if n_outliers > 0:
    print(df[df["OUTLIER"]][["FID", "IID", "HET"]].to_string(index=False))

# ── Write outlier list ────────────────────────────────────────────────────────
outliers = df[df["OUTLIER"]][["FID", "IID"]]
outliers.to_csv(OUT_LIST, sep=" ", header=False, index=False)
print(f"\nOutlier list written to: {OUT_LIST}")
if n_outliers == 0:
    print("  (empty — no samples to remove)")

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.suptitle(
    "Sample Heterozygosity QC — NOTCH1 Region, SABE Cohort (n=1,171)",
    fontsize=13, fontweight="bold", y=1.02
)

colors = ["#e74c3c" if o else "#3498db" for o in df["OUTLIER"]]

# Left: scatter plot (sample index vs het)
ax1 = axes[0]
ax1.scatter(range(len(df)), df["HET"], c=colors, s=8, alpha=0.6, linewidths=0)
ax1.axhline(mean_het,  color="black",  lw=1.5, linestyle="--", label=f"Mean ({mean_het:.3f})")
ax1.axhline(hi_thresh, color="#e74c3c", lw=1.5, linestyle=":",  label=f"+{SD_THRESH} SD ({hi_thresh:.3f})")
ax1.axhline(lo_thresh, color="#e74c3c", lw=1.5, linestyle=":",  label=f"−{SD_THRESH} SD ({lo_thresh:.3f})")
ax1.set_xlabel("Sample index", fontsize=11)
ax1.set_ylabel("Observed heterozygosity", fontsize=11)
ax1.set_title("Heterozygosity per sample", fontsize=11)
ax1.legend(fontsize=9)

normal_patch  = mpatches.Patch(color="#3498db", label=f"Normal (n={len(df)-n_outliers})")
outlier_patch = mpatches.Patch(color="#e74c3c", label=f"Outlier (n={n_outliers})")
ax1.legend(handles=[normal_patch, outlier_patch], fontsize=9, loc="upper right")

# Right: histogram
ax2 = axes[1]
ax2.hist(df[~df["OUTLIER"]]["HET"], bins=50, color="#3498db", alpha=0.7, label="Normal")
if n_outliers > 0:
    ax2.hist(df[df["OUTLIER"]]["HET"],  bins=20, color="#e74c3c", alpha=0.8, label="Outlier")
ax2.axvline(mean_het,  color="black",  lw=1.5, linestyle="--")
ax2.axvline(hi_thresh, color="#e74c3c", lw=1.5, linestyle=":")
ax2.axvline(lo_thresh, color="#e74c3c", lw=1.5, linestyle=":")
ax2.set_xlabel("Observed heterozygosity", fontsize=11)
ax2.set_ylabel("Number of samples", fontsize=11)
ax2.set_title("Distribution of heterozygosity", fontsize=11)
ax2.legend(fontsize=9)

# Annotation box
note = (
    f"NOTE: Higher mean heterozygosity ({mean_het:.3f}) compared to\n"
    f"homogeneous European cohorts (~0.14) is EXPECTED in admixed\n"
    f"Brazilian populations. This is not a QC failure."
)
ax2.text(
    0.02, 0.97, note,
    transform=ax2.transAxes,
    fontsize=8, verticalalignment="top",
    bbox=dict(boxstyle="round,pad=0.4", facecolor="lightyellow", alpha=0.8)
)

plt.tight_layout()
plt.savefig(OUT_FIGURE, dpi=150, bbox_inches="tight")
print(f"Figure saved to: {OUT_FIGURE}")

# ── Summary ───────────────────────────────────────────────────────────────────
print(f"\n{'='*50}")
print(f" HETEROZYGOSITY QC SUMMARY")
print(f"{'='*50}")
print(f"  Total samples  : {len(df)}")
print(f"  Outliers flagged: {n_outliers}")
print(f"  Samples retained: {len(df) - n_outliers}")
print(f"\nNext step: bash scripts/02_qc_pipeline.sh already ran LD pruning and PCA.")
print(f"           If outliers were found, re-run PCA after removing them.")
print(f"           Run: Rscript scripts/04_pca_analysis.R")
