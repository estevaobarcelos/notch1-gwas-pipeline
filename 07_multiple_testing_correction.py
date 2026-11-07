#!/usr/bin/env python3
"""
07_multiple_testing_correction.py
==================================
Multiple testing correction for NOTCH1 SNP association results.

Two methods applied:
  1. Benjamini-Hochberg (BH) FDR — controls false discovery rate.
     More appropriate than Bonferroni when tests are correlated (LD).
     Threshold: FDR q < 0.05 and q < 0.10.

  2. Permutation-based correction — gold standard for regional studies.
     Shuffles phenotype labels N times, recomputes association for each
     SNP, derives empirical p-values that automatically account for LD
     structure between SNPs.
     Default: 10,000 permutations.

Why not Bonferroni alone:
  Bonferroni assumes all tests are independent. In a 50-kb region,
  SNPs are correlated (LD), so Bonferroni overcorrects — it is too
  conservative. For 567 SNPs: threshold = 0.05/567 = 8.8e-5, which
  rs9411207 (p=0.0002) does NOT survive. BH and permutation correctly
  account for the correlation structure.

Inputs:
  notch1_model2.assoc.logistic   — Model 2 association results (PLINK)
  notch1_phenotyped.{bed,bim,fam} — genotype data for permutation

Outputs:
  results/multiple_testing_results.csv   — all SNPs with BH q-values
  results/permutation_results.csv        — empirical p-values per SNP
  figures/multiple_testing_plot.png      — BH + permutation summary figure

Dependencies: pandas, numpy, scipy, matplotlib

Author: Estevão Barcelos
Reference: Benjamini & Hochberg (1995), J R Stat Soc B
           Errera et al. Scientific Reports 2024
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings("ignore")

os.makedirs("results", exist_ok=True)
os.makedirs("figures", exist_ok=True)

# ── Configuration ─────────────────────────────────────────────────────────────
ASSOC_FILE  = "notch1_model2.assoc.logistic"
BIM_FILE    = "notch1_phenotyped.bim"
FAM_FILE    = "notch1_phenotyped.fam"
BED_FILE    = "notch1_phenotyped.bed"
N_PERMS     = 10000
FDR_THRESH  = 0.05
TARGET_SNP  = "rs9411207"
SEED        = 42
np.random.seed(SEED)


# ── Load association results ──────────────────────────────────────────────────
def load_assoc(path):
    if not os.path.exists(path):
        print(f"ERROR: {path} not found.")
        print("Run bash scripts/05_association_analysis.sh first.")
        sys.exit(1)
    df = pd.read_csv(path, sep=r"\s+")
    df = df[df["TEST"] == "ADD"].copy()
    df = df.dropna(subset=["P"])
    df = df[df["P"] > 0]
    return df.reset_index(drop=True)


# ── BH FDR correction ─────────────────────────────────────────────────────────
def bh_correction(p_values):
    """
    Benjamini-Hochberg procedure.
    Returns q-values (adjusted p-values) controlling FDR at given level.
    """
    n = len(p_values)
    p_sorted_idx = np.argsort(p_values)
    p_sorted     = p_values[p_sorted_idx]

    q_values = np.zeros(n)
    cummin   = np.inf
    for i in range(n - 1, -1, -1):
        q = p_sorted[i] * n / (i + 1)
        cummin = min(cummin, q)
        q_values[p_sorted_idx[i]] = cummin

    return np.minimum(q_values, 1.0)


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
        magic = f.read(3)
        if magic != b"\x6c\x1b\x01":
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


# ── Logistic regression test statistic ───────────────────────────────────────
def logistic_stat(geno_vec, pheno_vec):
    """
    Fast logistic regression z-statistic for permutation testing.
    Uses Wald test on the allele dosage coefficient.
    Returns |z| (absolute value for two-sided test).
    """
    mask = ~np.isnan(geno_vec) & ~np.isnan(pheno_vec) & (pheno_vec >= 0)
    g = geno_vec[mask]
    y = pheno_vec[mask]
    if g.std() == 0 or y.sum() < 5 or (1 - y).sum() < 5:
        return 0.0

    # Newton-Raphson logistic regression (1 predictor + intercept)
    X = np.column_stack([np.ones(len(g)), g])
    beta = np.zeros(2)
    for _ in range(20):
        eta  = X @ beta
        mu   = 1 / (1 + np.exp(-eta))
        W    = mu * (1 - mu)
        XtWX = X.T @ (W[:, None] * X)
        try:
            XtWX_inv = np.linalg.inv(XtWX + 1e-10 * np.eye(2))
        except np.linalg.LinAlgError:
            return 0.0
        grad = X.T @ (y - mu)
        beta_new = beta + XtWX_inv @ grad
        if np.max(np.abs(beta_new - beta)) < 1e-6:
            beta = beta_new
            break
        beta = beta_new

    se = np.sqrt(np.diag(XtWX_inv)[1])
    return abs(beta[1] / se) if se > 0 else 0.0


# ── Permutation test ──────────────────────────────────────────────────────────
def permutation_test(geno_matrix, pheno_vec, snp_names, assoc_df,
                     n_perms=10000):
    """
    Permutation-based multiple testing correction.
    For each permutation: shuffle phenotype, compute max statistic across
    all SNPs → build null distribution of max statistics.
    Empirical p-value = proportion of permutations where max stat >= observed.
    """
    # Only test SNPs that have association results
    snp_set  = set(assoc_df["SNP"].values)
    snp_idx  = [i for i, s in enumerate(snp_names) if s in snp_set]
    snp_sel  = [snp_names[i] for i in snp_idx]
    geno_sel = geno_matrix[snp_idx, :]

    # Observed statistics
    obs_stats = {}
    for i, snp in enumerate(snp_sel):
        obs_stats[snp] = logistic_stat(geno_sel[i], pheno_vec)

    # Build null distribution (max statistic per permutation)
    null_max = np.zeros(n_perms)
    pheno_valid = pheno_vec[~np.isnan(pheno_vec)]
    pheno_idx   = np.where(~np.isnan(pheno_vec))[0]

    print(f"  Running {n_perms:,} permutations ({len(snp_sel)} SNPs)...")
    for perm in range(n_perms):
        if (perm + 1) % 1000 == 0:
            print(f"    {perm+1:,} / {n_perms:,} done")
        perm_pheno          = pheno_vec.copy()
        shuffled            = pheno_valid.copy()
        np.random.shuffle(shuffled)
        perm_pheno[pheno_idx] = shuffled

        max_stat = 0.0
        for i in range(len(snp_sel)):
            stat = logistic_stat(geno_sel[i], perm_pheno)
            if stat > max_stat:
                max_stat = stat
        null_max[perm] = max_stat

    # Empirical p-values
    perm_results = []
    for snp in snp_sel:
        obs  = obs_stats[snp]
        p_emp = (null_max >= obs).mean()
        perm_results.append({"SNP": snp, "stat_obs": obs,
                              "p_empirical": max(p_emp, 1 / n_perms)})

    return pd.DataFrame(perm_results), null_max


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════
def main():
    print("=" * 60)
    print(" Multiple Testing Correction — NOTCH1 SNPs")
    print(f" Method 1: Benjamini-Hochberg FDR")
    print(f" Method 2: Permutation test ({N_PERMS:,} permutations)")
    print("=" * 60)

    # ── Load association results ──────────────────────────────────
    print("\nLoading association results...")
    assoc = load_assoc(ASSOC_FILE)
    n_snps = len(assoc)
    print(f"  {n_snps} SNPs with association results (Model 2, ADD test)")

    # ── BH correction ─────────────────────────────────────────────
    print("\nApplying Benjamini-Hochberg FDR correction...")
    p_vals   = assoc["P"].values.astype(float)
    q_vals   = bh_correction(p_vals)
    assoc["q_BH"]          = q_vals
    assoc["BH_sig_05"]     = q_vals < 0.05
    assoc["BH_sig_10"]     = q_vals < 0.10
    assoc["bonferroni_p"]  = np.minimum(p_vals * n_snps, 1.0)
    assoc["bonferroni_sig"]= assoc["bonferroni_p"] < 0.05

    n_bh05  = assoc["BH_sig_05"].sum()
    n_bh10  = assoc["BH_sig_10"].sum()
    n_bonf  = assoc["bonferroni_sig"].sum()

    print(f"  Bonferroni threshold  : p < {0.05/n_snps:.2e}  → {n_bonf} significant")
    print(f"  BH FDR q < 0.05       : {n_bh05} significant")
    print(f"  BH FDR q < 0.10       : {n_bh10} significant")

    # Target SNP
    if TARGET_SNP in assoc["SNP"].values:
        t = assoc[assoc["SNP"] == TARGET_SNP].iloc[0]
        print(f"\n  {TARGET_SNP}:")
        print(f"    Nominal p        = {t['P']:.6f}")
        print(f"    Bonferroni p     = {t['bonferroni_p']:.6f}  "
              f"({'sig' if t['bonferroni_sig'] else 'NOT sig'})")
        print(f"    BH q-value       = {t['q_BH']:.6f}  "
              f"({'sig q<0.05' if t['BH_sig_05'] else 'q<0.10' if t['BH_sig_10'] else 'NOT sig'})")

    # ── Permutation test ──────────────────────────────────────────
    print("\nLoading genotype data for permutation test...")
    for f in [BED_FILE, BIM_FILE, FAM_FILE]:
        if not os.path.exists(f):
            print(f"  WARNING: {f} not found — skipping permutation test.")
            print("  Copy PLINK files to project directory and re-run.")
            perm_df  = pd.DataFrame()
            null_max = np.array([])
            break
    else:
        bim, fam, geno = load_plink_bed(BED_FILE, BIM_FILE, FAM_FILE)
        pheno_raw = fam["pheno"].values.astype(float)
        pheno_vec = np.where(pheno_raw == 2, 1,
                    np.where(pheno_raw == 1, 0, np.nan))

        perm_df, null_max = permutation_test(
            geno, pheno_vec, bim["snp"].values, assoc, N_PERMS)

        # Merge permutation results
        assoc = assoc.merge(perm_df[["SNP","p_empirical"]], on="SNP", how="left")
        assoc["perm_sig_05"] = assoc["p_empirical"] < 0.05

        n_perm = (assoc["perm_sig_05"] == True).sum()
        print(f"\n  Permutation p < 0.05: {n_perm} significant SNPs")

        if TARGET_SNP in assoc["SNP"].values:
            t = assoc[assoc["SNP"] == TARGET_SNP].iloc[0]
            print(f"\n  {TARGET_SNP}:")
            print(f"    Empirical p      = {t.get('p_empirical', 'N/A'):.6f}  "
                  f"({'sig' if t.get('perm_sig_05', False) else 'NOT sig'})")

    # ── Save results ──────────────────────────────────────────────
    out_cols = ["CHR","SNP","BP","A1","NMISS","OR","L95","U95","P",
                "bonferroni_p","bonferroni_sig","q_BH","BH_sig_05","BH_sig_10"]
    if "p_empirical" in assoc.columns:
        out_cols += ["p_empirical","perm_sig_05"]

    assoc_out = assoc[[c for c in out_cols if c in assoc.columns]]
    assoc_out = assoc_out.sort_values("P")
    assoc_out.to_csv("results/multiple_testing_results.csv", index=False)
    print(f"\nFull results saved: results/multiple_testing_results.csv")

    # ── Figure ────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Panel 1: Nominal p vs BH q — volcano-style
    ax1 = axes[0]
    colors = []
    for _, row in assoc.iterrows():
        if row["SNP"] == TARGET_SNP:
            colors.append("#e74c3c")
        elif row.get("BH_sig_05", False):
            colors.append("#e67e22")
        elif row.get("BH_sig_10", False):
            colors.append("#f1c40f")
        else:
            colors.append("#bdc3c7")

    ax1.scatter(-np.log10(assoc["P"]),
                -np.log10(assoc["q_BH"]),
                c=colors, s=25, alpha=0.8)
    ax1.axhline(-np.log10(0.05), color="#e67e22", lw=1.2,
                linestyle="--", label="q = 0.05")
    ax1.axhline(-np.log10(0.10), color="#f1c40f", lw=1.2,
                linestyle=":", label="q = 0.10")
    ax1.axvline(-np.log10(0.05), color="gray", lw=0.8,
                linestyle="--", alpha=0.5)

    if TARGET_SNP in assoc["SNP"].values:
        t = assoc[assoc["SNP"] == TARGET_SNP].iloc[0]
        ax1.annotate(TARGET_SNP,
                     xy=(-np.log10(t["P"]), -np.log10(t["q_BH"])),
                     xytext=(5, 5), textcoords="offset points",
                     fontsize=7, color="#e74c3c", fontweight="bold")

    ax1.set_xlabel("−log₁₀(nominal p)", fontsize=10)
    ax1.set_ylabel("−log₁₀(BH q-value)", fontsize=10)
    ax1.set_title("Nominal p vs BH q-value\n(diagonal = no correction needed)",
                  fontsize=9, fontweight="bold")
    ax1.legend(fontsize=8)
    handles = [
        mpatches.Patch(color="#e74c3c", label=TARGET_SNP),
        mpatches.Patch(color="#e67e22", label="BH q < 0.05"),
        mpatches.Patch(color="#f1c40f", label="BH q < 0.10"),
        mpatches.Patch(color="#bdc3c7", label="Not significant"),
    ]
    ax1.legend(handles=handles, fontsize=7)

    # Panel 2: BH q-value distribution (sorted)
    ax2 = axes[1]
    q_sorted = np.sort(assoc["q_BH"].values)
    rank     = np.arange(1, len(q_sorted) + 1)
    ax2.plot(rank, q_sorted, color="#2980b9", lw=2, label="BH q-values")
    ax2.axhline(0.05, color="#e67e22", lw=1.5,
                linestyle="--", label="q = 0.05")
    ax2.axhline(0.10, color="#f1c40f", lw=1.5,
                linestyle=":", label="q = 0.10")
    ax2.fill_between(rank, 0, q_sorted,
                     where=(q_sorted < 0.10), alpha=0.15, color="#2980b9")
    ax2.set_xlabel("SNP rank (by p-value)", fontsize=10)
    ax2.set_ylabel("BH q-value", fontsize=10)
    ax2.set_title(f"BH q-value distribution\n{n_snps} SNPs tested (Model 2)",
                  fontsize=9, fontweight="bold")
    ax2.legend(fontsize=8)
    ax2.set_ylim(0, 1.05)

    # Panel 3: Permutation null distribution
    ax3 = axes[2]
    if len(null_max) > 0:
        ax3.hist(null_max, bins=60, color="#3498db", alpha=0.75,
                 edgecolor="white", lw=0.4, label="Null max |z| per permutation")

        if TARGET_SNP in assoc["SNP"].values and "stat_obs" in perm_df.columns:
            obs_stat_target = perm_df[perm_df["SNP"] == TARGET_SNP]["stat_obs"].values
            if len(obs_stat_target) > 0:
                ax3.axvline(obs_stat_target[0], color="#e74c3c",
                            lw=2, label=f"{TARGET_SNP} |z|={obs_stat_target[0]:.2f}")

        pct_95 = np.percentile(null_max, 95)
        ax3.axvline(pct_95, color="#e67e22", lw=1.5, linestyle="--",
                    label=f"95th percentile = {pct_95:.2f}")
        ax3.set_xlabel("Max |z| statistic across all SNPs", fontsize=10)
        ax3.set_ylabel("Frequency", fontsize=10)
        ax3.set_title(f"Permutation null distribution\n"
                      f"{N_PERMS:,} permutations",
                      fontsize=9, fontweight="bold")
        ax3.legend(fontsize=8)
    else:
        ax3.text(0.5, 0.5, "Permutation test\nnot run\n(PLINK files missing)",
                 ha="center", va="center", transform=ax3.transAxes,
                 fontsize=12, color="gray")
        ax3.set_title("Permutation null distribution", fontsize=9)

    plt.suptitle(
        "Multiple Testing Correction — NOTCH1 Region SNPs\n"
        "SABE Cohort | Model 2 (age, sex, ancestry, HDL, TG, hsCRP, HbA1c, hypertension)",
        fontsize=10, fontweight="bold", y=1.02
    )
    plt.tight_layout()
    plt.savefig("figures/multiple_testing_plot.png",
                dpi=150, bbox_inches="tight")
    print("Figure saved: figures/multiple_testing_plot.png")

    # ── Print summary table ───────────────────────────────────────
    print("\n" + "="*65)
    print(" MULTIPLE TESTING SUMMARY — Top SNPs")
    print("="*65)
    cols_print = ["SNP", "P", "bonferroni_p", "q_BH"]
    if "p_empirical" in assoc.columns:
        cols_print.append("p_empirical")
    top = assoc.sort_values("P").head(15)
    print(top[cols_print].to_string(index=False))

    print("\n" + "="*65)
    print(f" SIGNIFICANCE SUMMARY (n={n_snps} SNPs tested)")
    print("="*65)
    print(f"  Method              Threshold     Significant SNPs")
    print(f"  ─────────────────────────────────────────────────")
    print(f"  Nominal             p < 0.05      {(assoc['P']<0.05).sum()}")
    print(f"  Bonferroni          p < {0.05/n_snps:.2e}  {n_bonf}")
    print(f"  BH FDR              q < 0.05      {n_bh05}")
    print(f"  BH FDR              q < 0.10      {n_bh10}")
    if "p_empirical" in assoc.columns:
        print(f"  Permutation         p < 0.05      {n_perm}")

    print(f"\n  rs9411207 survives:")
    t = assoc[assoc["SNP"] == TARGET_SNP].iloc[0] if TARGET_SNP in assoc["SNP"].values else None
    if t is not None:
        print(f"    Bonferroni : {'✓' if t['bonferroni_sig'] else '✗'}")
        print(f"    BH q<0.05  : {'✓' if t['BH_sig_05'] else '✗'}")
        print(f"    BH q<0.10  : {'✓' if t['BH_sig_10'] else '✗'}")
        if "p_empirical" in assoc.columns:
            print(f"    Permutation: {'✓' if t.get('perm_sig_05', False) else '✗'}")

    print("\nStep 07 complete.")
    print("Next: python3 scripts/08_linear_regression_bmi.py")


if __name__ == "__main__":
    main()
