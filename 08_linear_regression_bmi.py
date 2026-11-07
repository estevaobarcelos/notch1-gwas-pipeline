#!/usr/bin/env python3
"""
08_linear_regression_bmi.py
============================
Linear regression: BMI-2010 (continuous) ~ NOTCH1 SNP genotypes.

This complements the logistic regression (case/control) by treating
BMI as a continuous outcome. It answers: does the risk genotype (TT)
predict a higher BMI value, beyond just being overweight vs not?

Models run for each SNP (additive dosage 0/1/2):
  1. Unadjusted: BMI ~ genotype
  2. Adjusted:   BMI ~ genotype + age + sex + EUR + AFR + NAM
                       + HDL + TG + hsCRP + HbA1c + hypertension

For rs9411207 specifically, also runs:
  3. Genotype-level: BMI by TT/TC/CC with ANOVA and post-hoc tests
  4. BMI distribution plot by genotype (violin + boxplot)
  5. Trend test: does BMI increase with T allele dosage?

Inputs:
  sample_full_metadata.csv          — BMI-2010 and ID per individual
  sample_covariates_model2.txt      — covariates
  notch1_phenotyped.{bed,bim,fam}  — genotype data

Outputs:
  results/linear_regression_all_snps.csv    — β, SE, p for all SNPs
  results/linear_regression_rs9411207.csv   — full rs9411207 results
  figures/bmi_by_genotype_rs9411207.png     — BMI distribution by genotype
  figures/linear_regression_volcano.png     — effect size vs p-value

Dependencies: pandas, numpy, scipy, statsmodels, matplotlib, seaborn

Author: Estevão Barcelos
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings("ignore")

try:
    import seaborn as sns
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False

os.makedirs("results", exist_ok=True)
os.makedirs("figures", exist_ok=True)

TARGET_SNP = "rs9411207"
BIM_FILE   = "notch1_phenotyped.bim"
FAM_FILE   = "notch1_phenotyped.fam"
BED_FILE   = "notch1_phenotyped.bed"
META_FILE  = "sample_full_metadata.csv"
COV_FILE   = "sample_covariates_model2.txt"


# ── Load PLINK binary ─────────────────────────────────────────────────────────
def load_plink_bed(bed, bim, fam):
    bim_df = pd.read_csv(bim, sep="\t", header=None,
                         names=["chr","snp","cm","pos","a1","a2"])
    fam_df = pd.read_csv(fam, sep=r"\s+", header=None,
                         names=["fid","iid","pid","mid","sex","pheno"])
    n_snps    = len(bim_df)
    n_samples = len(fam_df)
    n_bytes   = (n_samples + 3) // 4

    with open(bed, "rb") as f:
        if f.read(3) != b"\x6c\x1b\x01":
            raise ValueError("Not a valid PLINK .bed file")
        raw = np.frombuffer(f.read(), dtype=np.uint8)

    geno = np.full((n_snps, n_samples), np.nan, dtype=np.float32)
    for snp_idx in range(n_snps):
        byte_block = raw[snp_idx * n_bytes: (snp_idx + 1) * n_bytes]
        for byte_idx, byte in enumerate(byte_block):
            for bit_pair in range(4):
                sample_idx = byte_idx * 4 + bit_pair
                if sample_idx >= n_samples:
                    break
                code = (byte >> (bit_pair * 2)) & 0b11
                if   code == 0b00: geno[snp_idx, sample_idx] = 0
                elif code == 0b10: geno[snp_idx, sample_idx] = 1
                elif code == 0b11: geno[snp_idx, sample_idx] = 2
    return bim_df, fam_df, geno


# ── Linear regression ─────────────────────────────────────────────────────────
def linear_model(y, X_design, snp_name="SNP"):
    """OLS regression. Returns β, SE, t, p for the SNP term (column 1)."""
    mask = ~np.isnan(y) & ~np.any(np.isnan(X_design), axis=1)
    y_v = y[mask]
    X_v = X_design[mask]
    if mask.sum() < 20 or X_v[:, 1].std() == 0:
        return dict(snp=snp_name, n=mask.sum(), beta=np.nan,
                    se=np.nan, t=np.nan, p=np.nan, r2=np.nan)
    try:
        model = sm.OLS(y_v, X_v).fit()
        return dict(
            snp  = snp_name,
            n    = int(mask.sum()),
            beta = model.params[1],
            se   = model.bse[1],
            t    = model.tvalues[1],
            p    = model.pvalues[1],
            r2   = model.rsquared,
            r2_adj = model.rsquared_adj,
        )
    except Exception:
        return dict(snp=snp_name, n=mask.sum(), beta=np.nan,
                    se=np.nan, t=np.nan, p=np.nan, r2=np.nan)


def main():
    print("=" * 60)
    print(" Linear Regression: BMI-2010 ~ NOTCH1 SNP Genotypes")
    print("=" * 60)

    # ── Load data ─────────────────────────────────────────────────
    for f in [BED_FILE, BIM_FILE, FAM_FILE, META_FILE]:
        if not os.path.exists(f):
            print(f"ERROR: {f} not found.")
            sys.exit(1)

    print("\nLoading data...")
    bim, fam, geno = load_plink_bed(BED_FILE, BIM_FILE, FAM_FILE)
    n_snps, n_samples = geno.shape
    print(f"  {n_snps} SNPs × {n_samples} samples")

    # BMI from metadata
    meta = pd.read_csv(META_FILE)
    # Align BMI to fam order by IID
    # fam IID is "176501664_176501664"; meta FID is "176501664_176501664"
    bmi_map = dict(zip(meta["FID"].astype(str), meta["BMI_2010"]))
    bmi_vec = np.array([
        bmi_map.get(str(fam.loc[i, "fid"]), np.nan)
        for i in range(len(fam))
    ], dtype=float)

    n_bmi = (~np.isnan(bmi_vec)).sum()
    bmi_valid = bmi_vec[~np.isnan(bmi_vec)]
    print(f"  BMI-2010 available: {n_bmi} individuals")
    print(f"  BMI range: {bmi_valid.min():.1f} – {bmi_valid.max():.1f}")
    print(f"  BMI mean ± SD: {bmi_valid.mean():.1f} ± {bmi_valid.std():.1f}")

    # Covariates
    covar_exists = os.path.exists(COV_FILE)
    if covar_exists:
        cov_df = pd.read_csv(COV_FILE, sep=r"\s+")
        cov_map = {str(r["FID"]): r for _, r in cov_df.iterrows()}
        cov_cols = ["AGE","SEX","EUR","AFR","NAM","HDL","TG","hsCRP","HbA1c","Hypertension"]
        cov_cols = [c for c in cov_cols if c in cov_df.columns]
        print(f"  Covariates loaded: {cov_cols}")
    else:
        print("  WARNING: Covariate file not found — running unadjusted only")

    # ── Run linear regression for all SNPs ───────────────────────
    print(f"\nRunning linear regression for {n_snps} SNPs...")
    results_unadj = []
    results_adj   = []

    for i, row in bim.iterrows():
        g = geno[i].astype(float)
        snp = row["snp"]

        # Unadjusted: BMI ~ intercept + genotype
        X_u = sm.add_constant(g)
        res_u = linear_model(bmi_vec, X_u, snp)
        res_u["pos"] = row["pos"]
        res_u["a1"]  = row["a1"]
        res_u["a2"]  = row["a2"]
        results_unadj.append(res_u)

        # Adjusted: BMI ~ intercept + genotype + covariates
        if covar_exists:
            cov_matrix = []
            for j in range(n_samples):
                fid_j = str(fam.loc[j, "fid"])
                if fid_j in cov_map:
                    r = cov_map[fid_j]
                    row_cov = [r.get(c, np.nan) for c in cov_cols]
                else:
                    row_cov = [np.nan] * len(cov_cols)
                cov_matrix.append(row_cov)
            cov_arr = np.array(cov_matrix, dtype=float)
            X_a = np.column_stack([np.ones(n_samples), g, cov_arr])
            res_a = linear_model(bmi_vec, X_a, snp)
            res_a["pos"] = row["pos"]
            res_a["a1"]  = row["a1"]
            res_a["a2"]  = row["a2"]
            results_adj.append(res_a)

    df_unadj = pd.DataFrame(results_unadj).sort_values("p")
    df_adj   = pd.DataFrame(results_adj).sort_values("p") if results_adj else pd.DataFrame()

    # BH correction on linear results
    for df, label in [(df_unadj, "unadj"), (df_adj, "adj")]:
        if df.empty: continue
        p_valid = df["p"].fillna(1.0).values
        n = len(p_valid)
        idx_sorted = np.argsort(p_valid)
        q = np.zeros(n); cummin = np.inf
        for k in range(n - 1, -1, -1):
            q_k = p_valid[idx_sorted[k]] * n / (k + 1)
            cummin = min(cummin, q_k)
            q[idx_sorted[k]] = cummin
        df["q_BH"] = np.minimum(q, 1.0)
        df["BH_sig"] = df["q_BH"] < 0.05

    # Save
    df_unadj.to_csv("results/linear_regression_all_snps_unadj.csv", index=False)
    if not df_adj.empty:
        df_adj.to_csv("results/linear_regression_all_snps_adj.csv", index=False)
    print(f"  Linear regression results saved.")

    # ── Focus on rs9411207 ────────────────────────────────────────
    print(f"\n── rs9411207 — BMI linear regression ──")
    target_idx = bim[bim["snp"] == TARGET_SNP].index
    if len(target_idx) == 0:
        print(f"  WARNING: {TARGET_SNP} not found in dataset.")
    else:
        idx = target_idx[0]
        g_target = geno[idx].astype(float)

        # Unadjusted
        X_u = sm.add_constant(g_target)
        mask_u = ~np.isnan(bmi_vec) & ~np.isnan(g_target)
        model_u = sm.OLS(bmi_vec[mask_u], X_u[mask_u]).fit()

        print(f"  Unadjusted (BMI ~ genotype dosage):")
        print(f"    β = {model_u.params[1]:.3f} kg/m²  per T allele")
        print(f"    SE = {model_u.bse[1]:.3f}")
        print(f"    t  = {model_u.tvalues[1]:.3f}")
        print(f"    p  = {model_u.pvalues[1]:.5f}")
        print(f"    R² = {model_u.rsquared:.4f}")

        # Adjusted
        if covar_exists:
            cov_matrix = []
            for j in range(n_samples):
                fid_j = str(fam.loc[j, "fid"])
                if fid_j in cov_map:
                    r = cov_map[fid_j]
                    row_cov = [r.get(c, np.nan) for c in cov_cols]
                else:
                    row_cov = [np.nan] * len(cov_cols)
                cov_matrix.append(row_cov)
            cov_arr = np.array(cov_matrix, dtype=float)
            X_a = np.column_stack([np.ones(n_samples), g_target, cov_arr])
            mask_a = ~np.isnan(bmi_vec) & ~np.any(np.isnan(X_a), axis=1)
            model_a = sm.OLS(bmi_vec[mask_a], X_a[mask_a]).fit()

            print(f"\n  Adjusted (BMI ~ genotype + covariates):")
            print(f"    β = {model_a.params[1]:.3f} kg/m²  per T allele")
            print(f"    SE = {model_a.bse[1]:.3f}")
            print(f"    t  = {model_a.tvalues[1]:.3f}")
            print(f"    p  = {model_a.pvalues[1]:.5f}")
            print(f"    R² = {model_a.rsquared:.4f}")

        # Genotype-level BMI stats (CC / TC / TT)
        genotype_labels = {0.0: "CC", 1.0: "TC", 2.0: "TT"}
        bmi_by_gt = {}
        for gt_code, gt_label in genotype_labels.items():
            mask = (g_target == gt_code) & ~np.isnan(bmi_vec)
            bmi_by_gt[gt_label] = bmi_vec[mask]

        print(f"\n  BMI by genotype:")
        print(f"    {'Genotype':<8} {'N':>5} {'Mean':>8} {'SD':>7} {'Median':>8}")
        print(f"    {'─'*42}")
        for gt in ["CC","TC","TT"]:
            vals = bmi_by_gt[gt]
            print(f"    {gt:<8} {len(vals):>5} {vals.mean():>8.2f} "
                  f"{vals.std():>7.2f} {np.median(vals):>8.2f}")

        # ANOVA
        f_stat, p_anova = stats.f_oneway(
            bmi_by_gt["CC"], bmi_by_gt["TC"], bmi_by_gt["TT"])
        print(f"\n  One-way ANOVA: F={f_stat:.3f}, p={p_anova:.5f}")

        # Pairwise t-tests (Bonferroni-corrected within 3 pairs)
        pairs = [("CC","TC"), ("CC","TT"), ("TC","TT")]
        print(f"  Pairwise t-tests (Bonferroni × 3):")
        for g1, g2 in pairs:
            t, p = stats.ttest_ind(bmi_by_gt[g1], bmi_by_gt[g2])
            sig = "★" if p * 3 < 0.05 else ""
            print(f"    {g1} vs {g2}: t={t:.3f}, p={p:.5f} "
                  f"(adj p={p*3:.5f}) {sig}")

        # Trend test (Jonckheere-Terpstra via Spearman correlation)
        gt_vec  = g_target[~np.isnan(bmi_vec) & ~np.isnan(g_target)]
        bmi_vec2 = bmi_vec[~np.isnan(bmi_vec) & ~np.isnan(g_target)]
        rho, p_trend = stats.spearmanr(gt_vec, bmi_vec2)
        print(f"\n  Allele dosage trend (Spearman):")
        print(f"    ρ = {rho:.4f},  p = {p_trend:.5f}")
        trend_dir = "↑ BMI increases with T allele count" \
                    if rho > 0 else "↓ BMI decreases with T allele count"
        print(f"    {trend_dir}")

        # Save rs9411207 summary
        rs_summary = pd.DataFrame([
            {"Model": "Unadjusted", "N": int(mask_u.sum()),
             "beta_per_T_allele": model_u.params[1], "SE": model_u.bse[1],
             "t": model_u.tvalues[1], "p": model_u.pvalues[1],
             "R2": model_u.rsquared},
        ])
        if covar_exists:
            rs_summary = pd.concat([rs_summary, pd.DataFrame([{
                "Model": "Adjusted (Model 2 covariates)", "N": int(mask_a.sum()),
                "beta_per_T_allele": model_a.params[1], "SE": model_a.bse[1],
                "t": model_a.tvalues[1], "p": model_a.pvalues[1],
                "R2": model_a.rsquared,
            }])], ignore_index=True)
        rs_summary.to_csv("results/linear_regression_rs9411207.csv", index=False)

        # ── Figure: BMI by genotype ───────────────────────────────
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))

        colors_gt = {"CC": "#3498db", "TC": "#95a5a6", "TT": "#e74c3c"}

        # Panel 1: Violin + boxplot
        ax1 = axes[0]
        positions = [0, 1, 2]
        gt_order  = ["CC", "TC", "TT"]
        data_plot = [bmi_by_gt[gt] for gt in gt_order]

        if HAS_SEABORN:
            df_plot = pd.DataFrame({
                "Genotype": np.concatenate([[gt]*len(bmi_by_gt[gt]) for gt in gt_order]),
                "BMI":      np.concatenate([bmi_by_gt[gt] for gt in gt_order])
            })
            sns.violinplot(data=df_plot, x="Genotype", y="BMI",
                           palette=colors_gt, ax=ax1, inner=None, alpha=0.5)
            sns.boxplot(data=df_plot, x="Genotype", y="BMI",
                        palette=colors_gt, ax=ax1, width=0.3,
                        fliersize=2, linewidth=1.2)
        else:
            vp = ax1.violinplot(data_plot, positions=positions,
                                showmedians=False, showextrema=False)
            for pc, gt in zip(vp["bodies"], gt_order):
                pc.set_facecolor(colors_gt[gt]); pc.set_alpha(0.5)
            bp = ax1.boxplot(data_plot, positions=positions,
                             widths=0.3, patch_artist=True, medianprops=dict(lw=2))
            for patch, gt in zip(bp["boxes"], gt_order):
                patch.set_facecolor(colors_gt[gt]); patch.set_alpha(0.7)
            ax1.set_xticks(positions)
            ax1.set_xticklabels([f"{gt}\n(n={len(bmi_by_gt[gt])})" for gt in gt_order])

        ax1.set_ylabel("BMI-2010 (kg/m²)", fontsize=11)
        ax1.set_title(f"BMI by rs9411207 genotype\n"
                      f"ANOVA: F={f_stat:.2f}, p={p_anova:.4f}",
                      fontsize=10, fontweight="bold")

        # Significance brackets
        y_max = max(np.percentile(bmi_by_gt[gt], 95) for gt in gt_order)
        y_br  = y_max + 1.5
        _, p_cc_tt = stats.ttest_ind(bmi_by_gt["CC"], bmi_by_gt["TT"])
        if p_cc_tt < 0.05:
            ax1.plot([0, 2], [y_br, y_br], "k-", lw=1.2)
            ax1.plot([0, 0], [y_max + 0.5, y_br], "k-", lw=1.2)
            ax1.plot([2, 2], [y_max + 0.5, y_br], "k-", lw=1.2)
            stars = "***" if p_cc_tt < 0.001 else "**" if p_cc_tt < 0.01 else "*"
            ax1.text(1, y_br + 0.3, stars, ha="center", fontsize=12)

        # Panel 2: Mean BMI ± SE by genotype
        ax2 = axes[1]
        means = [bmi_by_gt[gt].mean() for gt in gt_order]
        sems  = [bmi_by_gt[gt].std() / np.sqrt(len(bmi_by_gt[gt])) for gt in gt_order]
        ns    = [len(bmi_by_gt[gt]) for gt in gt_order]
        bar_colors = [colors_gt[gt] for gt in gt_order]
        bars = ax2.bar(gt_order, means, yerr=sems, color=bar_colors,
                       alpha=0.8, capsize=5, error_kw=dict(lw=1.5))
        for bar, mean, n in zip(bars, means, ns):
            ax2.text(bar.get_x() + bar.get_width()/2,
                     mean + sems[means.index(mean)] + 0.3,
                     f"{mean:.1f}", ha="center", va="bottom", fontsize=9)
        ax2.set_ylabel("Mean BMI ± SE (kg/m²)", fontsize=11)
        ax2.set_title(f"Mean BMI by genotype\n"
                      f"Trend: ρ={rho:.3f}, p={p_trend:.4f}",
                      fontsize=10, fontweight="bold")
        ax2.set_ylim(min(means) - 3, max(means) + 4)

        # Arrow showing trend
        if rho > 0:
            ax2.annotate("", xy=(2, means[2]+0.5),
                         xytext=(0, means[0]+0.5),
                         arrowprops=dict(arrowstyle="->", color="#e74c3c", lw=1.5))
            ax2.text(1, max(means)+2.5, "↑ T allele → higher BMI",
                     ha="center", fontsize=9, color="#e74c3c", style="italic")

        # Panel 3: Linear regression plot (scatter + line)
        ax3 = axes[2]
        jitter_x = g_target + np.random.normal(0, 0.07, len(g_target))
        mask_plot = ~np.isnan(bmi_vec) & ~np.isnan(g_target)
        scatter_colors = [colors_gt.get(genotype_labels.get(g, ""), "#bdc3c7")
                          for g in g_target[mask_plot]]
        ax3.scatter(jitter_x[mask_plot], bmi_vec[mask_plot],
                    c=scatter_colors, s=8, alpha=0.3)

        # Regression line
        x_line = np.array([0, 1, 2])
        y_line = model_u.params[0] + model_u.params[1] * x_line
        ax3.plot(x_line, y_line, color="#2c3e50", lw=2.5,
                 label=f"β={model_u.params[1]:.3f} kg/m² per T allele\np={model_u.pvalues[1]:.4f}")

        ax3.set_xticks([0, 1, 2])
        ax3.set_xticklabels(["CC", "TC", "TT"], fontsize=11)
        ax3.set_xlabel("rs9411207 genotype", fontsize=11)
        ax3.set_ylabel("BMI-2010 (kg/m²)", fontsize=11)
        ax3.set_title("Linear regression\nBMI ~ genotype dosage (unadjusted)",
                      fontsize=10, fontweight="bold")
        ax3.legend(fontsize=9, loc="upper left")

        plt.suptitle(
            "rs9411207 (NOTCH1) — BMI-2010 Linear Regression\n"
            "SABE Cohort, São Paulo, Brazil",
            fontsize=11, fontweight="bold", y=1.02
        )
        plt.tight_layout()
        plt.savefig("figures/bmi_by_genotype_rs9411207.png",
                    dpi=150, bbox_inches="tight")
        print("\nFigure saved: figures/bmi_by_genotype_rs9411207.png")

    # ── Figure 2: Linear regression volcano (all SNPs) ───────────
    if not df_unadj.empty:
        fig2, ax = plt.subplots(figsize=(10, 5))
        df_plot = df_unadj.dropna(subset=["beta","p"])
        colors = ["#e74c3c" if s == TARGET_SNP
                  else "#e67e22" if q < 0.05
                  else "#2980b9" if p < 0.05
                  else "#bdc3c7"
                  for s, q, p in zip(df_plot["snp"],
                                      df_plot.get("q_BH", [1]*len(df_plot)),
                                      df_plot["p"])]
        ax.scatter(df_plot["beta"], -np.log10(df_plot["p"]),
                   c=colors, s=20, alpha=0.8)
        ax.axhline(-np.log10(0.05), color="gray", lw=1, linestyle="--",
                   label="p = 0.05")
        ax.axvline(0, color="black", lw=0.8, linestyle="-", alpha=0.3)

        if TARGET_SNP in df_plot["snp"].values:
            t = df_plot[df_plot["snp"] == TARGET_SNP].iloc[0]
            ax.annotate(TARGET_SNP,
                        xy=(t["beta"], -np.log10(t["p"])),
                        xytext=(8, 8), textcoords="offset points",
                        fontsize=8, color="#e74c3c", fontweight="bold",
                        arrowprops=dict(arrowstyle="->", color="#e74c3c", lw=1))

        ax.set_xlabel("β coefficient (BMI kg/m² per T allele)", fontsize=11)
        ax.set_ylabel("−log₁₀(p)", fontsize=11)
        ax.set_title("Linear regression — BMI ~ genotype\nAll NOTCH1 region SNPs (unadjusted)",
                     fontsize=11, fontweight="bold")
        handles = [
            mpatches.Patch(color="#e74c3c", label=TARGET_SNP),
            mpatches.Patch(color="#e67e22", label="BH q < 0.05"),
            mpatches.Patch(color="#2980b9", label="p < 0.05"),
            mpatches.Patch(color="#bdc3c7", label="Not significant"),
        ]
        ax.legend(handles=handles, fontsize=9)
        plt.tight_layout()
        plt.savefig("figures/linear_regression_volcano.png",
                    dpi=150, bbox_inches="tight")
        print("Volcano plot saved: figures/linear_regression_volcano.png")

    # ── Top linear association results ────────────────────────────
    print("\n" + "="*60)
    print(" TOP 10 SNPs — BMI linear regression (unadjusted)")
    print("="*60)
    cols_p = ["snp","pos","a1","beta","se","p","q_BH"]
    cols_p = [c for c in cols_p if c in df_unadj.columns]
    print(df_unadj[cols_p].head(10).to_string(index=False))

    print("\nStep 08 complete.")


if __name__ == "__main__":
    main()
