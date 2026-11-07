#!/usr/bin/env python3
"""
06_results_summary.py
=====================
Generate a comprehensive QC summary table from PLINK log files.
Saves results/qc_summary_table.csv for the README and manuscript methods.

Author: Estevão Barcelos
"""

import os
import re
import pandas as pd

os.makedirs("results", exist_ok=True)

# ── Define pipeline steps and expected files ──────────────────────────────────
steps = [
    {
        "step": "01. Raw VCF input",
        "bim":  "notch1_plink.bim",
        "fam":  "notch1_plink.fam",
        "filter": "—",
        "note": "WGS, VQSR-filtered, biallelic SNPs"
    },
    {
        "step": "02a. Call rate (SNP)",
        "bim":  "notch1_cr_snp.bim",
        "fam":  "notch1_cr_snp.fam",
        "filter": "geno > 5% missing",
        "note": "Removes low-quality genotyped sites"
    },
    {
        "step": "02b. Call rate (sample)",
        "bim":  "notch1_cr_ind.bim",
        "fam":  "notch1_cr_ind.fam",
        "filter": "mind > 2% missing",
        "note": "Removes low-quality samples"
    },
    {
        "step": "02c. MAF filter",
        "bim":  "notch1_maf.bim",
        "fam":  "notch1_maf.fam",
        "filter": "MAF < 1%",
        "note": "Rare variants excluded from common-variant analysis"
    },
    {
        "step": "02d. HWE filter",
        "bim":  "notch1_hwe.bim",
        "fam":  "notch1_hwe.fam",
        "filter": "HWE p < 1×10⁻⁶",
        "note": "Applied to full sample (controls only when pheno available)"
    },
    {
        "step": "03. Heterozygosity QC",
        "bim":  "notch1_final_3sd.bim",
        "fam":  "notch1_final_3sd.fam",
        "filter": "±3 SD",
        "note": "Outlier removal; ±3 SD is the statistical standard"
    },
]

# ── Build summary table ───────────────────────────────────────────────────────
rows = []
prev_variants = None
prev_samples  = None

for s in steps:
    n_variants = None
    n_samples  = None

    if os.path.exists(s["bim"]):
        with open(s["bim"]) as f:
            n_variants = sum(1 for _ in f)
    if os.path.exists(s["fam"]):
        with open(s["fam"]) as f:
            n_samples = sum(1 for _ in f)

    removed_v = (prev_variants - n_variants) if (prev_variants and n_variants) else None
    removed_s = (prev_samples  - n_samples)  if (prev_samples  and n_samples)  else None

    rows.append({
        "Step":               s["step"],
        "Filter":             s["filter"],
        "Variants":           n_variants if n_variants else "N/A",
        "Samples":            n_samples  if n_samples  else "N/A",
        "Variants removed":   removed_v  if removed_v is not None else "—",
        "Samples removed":    removed_s  if removed_s is not None else "—",
        "Note":               s["note"],
    })

    if n_variants: prev_variants = n_variants
    if n_samples:  prev_samples  = n_samples

df = pd.DataFrame(rows)

# ── Save CSV ──────────────────────────────────────────────────────────────────
df.to_csv("results/qc_summary_table.csv", index=False)
print("QC summary table saved to: results/qc_summary_table.csv")

# ── Print to terminal ─────────────────────────────────────────────────────────
print("\n" + "="*80)
print(" QC PIPELINE SUMMARY — NOTCH1 REGION, SABE COHORT")
print("="*80)
print(df.to_string(index=False))
print("="*80)

# ── Genotyping rate from final log ────────────────────────────────────────────
log_file = "notch1_final_3sd.log"
if os.path.exists(log_file):
    with open(log_file) as f:
        for line in f:
            if "genotyping rate" in line.lower():
                print(f"\nFinal {line.strip()}")
                break

print("\nSummary complete. See results/qc_summary_table.csv")
