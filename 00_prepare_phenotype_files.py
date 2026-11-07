#!/usr/bin/env python3
"""
00_prepare_phenotype_files.py
=============================
Prepare all PLINK-ready input files from the three SABE cohort source files:
  - Keys.xlsx     : HLI ↔ WGS ID correspondence
  - Pheno.xlsx    : Full phenotypic and clinical data (189 variables)
  - Ancestries.xlsx: Global ancestry proportions (EUR, AFR, NAM, EAS)

Case/Control definition (from published study abstract):
  CASE    = Excessive weight → BMI ≥ 25 (overweight + obese) → PLINK code 2
  CONTROL = Normal weight   → BMI < 25                       → PLINK code 1
  MISSING = No BMI available                                  → PLINK code -9

BMI source:
  BMI-calc-2010 ONLY (most complete wave, n=1,024 with valid data)
  Extreme outliers excluded: BMI < 10 or BMI > 100
  Individuals without 2010 BMI → PHENOTYPE = -9 (missing, excluded by PLINK)

Outputs (PLINK-ready format):
  sample_pheno.txt       — FID IID PHENOTYPE
  sample_sex.txt         — FID IID SEX (1=M, 2=F)
  sample_covariates.txt  — FID IID AGE SEX EUR AFR NAM EAS
  sample_full_metadata.csv — full merged dataset for reference

Author: Estevão Barcelos
Cohort: SABE (Saúde, Bem-Estar e Envelhecimento), São Paulo, Brazil
"""

import os
import sys
import numpy as np
import pandas as pd
import openpyxl

# ── Input files ───────────────────────────────────────────────────────────────
KEYS_FILE  = "Keys.xlsx"
PHENO_FILE = "Pheno.xlsx"
ANC_FILE   = "Ancestries.xlsx"

# ── Output files ──────────────────────────────────────────────────────────────
OUT_PHENO  = "sample_pheno.txt"
OUT_SEX    = "sample_sex.txt"
OUT_COV    = "sample_covariates.txt"
OUT_META   = "sample_full_metadata.csv"

# ── BMI cutoff ────────────────────────────────────────────────────────────────
BMI_CASE_THRESHOLD = 25.0  # BMI >= 25 → excessive weight (case)


def load_pheno(path):
    """Load Pheno.xlsx with deduplicated column names."""
    wb = openpyxl.load_workbook(path, read_only=True)
    rows = list(wb.active.iter_rows(values_only=True))
    seen = {}
    cols = []
    for i, c in enumerate(rows[0]):
        c = str(c) if c is not None else f"col_{i}"
        seen[c] = seen.get(c, 0) + 1
        cols.append(c if seen[c] == 1 else f"{c}_{seen[c]}")
    return pd.DataFrame(rows[1:], columns=cols)


def load_keys(path):
    """Load Keys.xlsx: HLI ↔ WGS ID mapping."""
    wb = openpyxl.load_workbook(path, read_only=True)
    rows = list(wb.active.iter_rows(values_only=True))
    return pd.DataFrame(rows[1:], columns=["HLI", "WGS_ID"])


def load_ancestries(path):
    """Load Ancestries.xlsx: EUR, AFR, NAM, EAS proportions per individual."""
    wb = openpyxl.load_workbook(path, read_only=True)
    rows = list(wb.active.iter_rows(values_only=True))
    df = pd.DataFrame(rows[1:], columns=["ID", "_drop", "EUR", "AFR", "NAM", "EAS"])
    return df.drop(columns=["_drop"])


def parse_bmi(val):
    """Parse BMI value, returning NaN for missing or physiologically impossible values."""
    try:
        f = float(val)
        return f if 10.0 <= f <= 100.0 else np.nan
    except (TypeError, ValueError):
        return np.nan


def main():
    # Check inputs
    for f in [KEYS_FILE, PHENO_FILE, ANC_FILE]:
        if not os.path.exists(f):
            print(f"ERROR: {f} not found in current directory.")
            sys.exit(1)

    print("=" * 55)
    print(" Preparing PLINK input files — SABE Cohort")
    print("=" * 55)

    # ── Load ──────────────────────────────────────────────────────────────────
    print("\nLoading source files...")
    pheno = load_pheno(PHENO_FILE)
    keys  = load_keys(KEYS_FILE)
    anc   = load_ancestries(ANC_FILE)

    pheno["HLI-internal"] = pheno["HLI-internal"].astype(str)
    keys["HLI"]           = keys["HLI"].astype(str)
    anc["ID"]             = anc["ID"].astype(str)
    print(f"  Pheno      : {len(pheno)} individuals, {len(pheno.columns)} variables")
    print(f"  Keys       : {len(keys)} ID mappings")
    print(f"  Ancestries : {len(anc)} individuals")

    # ── BMI 2010 ONLY → Phenotype ─────────────────────────────────────────────
    print("\nDefining case/control status from BMI-2010 only...")
    pheno["BMI_2010"] = pheno["BMI-calc-2010"].apply(parse_bmi)

    # PLINK codes: 2=case, 1=control, -9=missing (no 2010 BMI → excluded)
    pheno["PHENOTYPE"] = np.where(
        pheno["BMI_2010"].isna(), -9,
        np.where(pheno["BMI_2010"] >= BMI_CASE_THRESHOLD, 2, 1)
    )
    pheno["SEX_code"] = pheno["Gender"].map({"M": 1, "F": 2}).fillna(0).astype(int)
    pheno["FID"] = pheno["HLI-internal"]
    pheno["IID"] = pheno["HLI-internal"]

    n_cases    = (pheno["PHENOTYPE"] == 2).sum()
    n_controls = (pheno["PHENOTYPE"] == 1).sum()
    n_missing  = (pheno["PHENOTYPE"] == -9).sum()
    print(f"  With BMI-2010 available  : {n_cases + n_controls}")
    print(f"  Missing BMI-2010 (→ -9)  : {n_missing}")
    print(f"  Cases (BMI ≥ {BMI_CASE_THRESHOLD})       : {n_cases}")
    print(f"  Controls (BMI < {BMI_CASE_THRESHOLD})    : {n_controls}")

    # ── Merge ─────────────────────────────────────────────────────────────────
    df = pheno.merge(keys, left_on="HLI-internal", right_on="HLI", how="inner")
    df = df.merge(anc,   left_on="HLI-internal", right_on="ID",  how="left")
    for col in ["EUR", "AFR", "NAM", "EAS"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    print(f"\nAfter merging all sources: {len(df)} individuals")

    # ── Write output files ────────────────────────────────────────────────────
    print("\nWriting PLINK input files...")

    # 1. Phenotype
    df[["FID", "IID", "PHENOTYPE"]].to_csv(
        OUT_PHENO, sep=" ", index=False, header=False)
    print(f"  ✓ {OUT_PHENO}")

    # 2. Sex
    df[["FID", "IID", "SEX_code"]].to_csv(
        OUT_SEX, sep=" ", index=False, header=False)
    print(f"  ✓ {OUT_SEX}")

    # 3. Covariates: Age, Sex, Ancestry (mean-impute missing ancestry)
    cov = df[["FID", "IID", "Age-at-blood-collection",
              "SEX_code", "EUR", "AFR", "NAM", "EAS"]].copy()
    cov.columns = ["FID", "IID", "AGE", "SEX", "EUR", "AFR", "NAM", "EAS"]
    cov["AGE"] = pd.to_numeric(cov["AGE"], errors="coerce").round(2)
    for col in ["EUR", "AFR", "NAM", "EAS"]:
        cov[col] = cov[col].fillna(cov[col].median())
    cov.to_csv(OUT_COV, sep=" ", index=False)
    print(f"  ✓ {OUT_COV}  (covariates: AGE SEX EUR AFR NAM EAS)")

    # 4b. BMI quantitative phenotype file for PLINK --linear
    # Format: FID IID BMI_value  (-9 for missing)
    bmi_quant = df[["FID", "IID", "BMI_2010"]].copy()
    bmi_quant["BMI_2010"] = bmi_quant["BMI_2010"].fillna(-9)
    bmi_quant.to_csv("sample_bmi_quantitative.txt",
                     sep=" ", index=False, header=False)
    print(f"  ✓ sample_bmi_quantitative.txt  (BMI-2010 as continuous phenotype for --linear)")

    # 4. Full metadata for reference
    meta = df[["FID", "IID", "WGS_ID", "Gender", "SEX_code",
               "Age-at-blood-collection", "BMI_2010", "PHENOTYPE",
               "EUR", "AFR", "NAM", "EAS"]].copy()
    meta.columns = ["FID", "IID", "WGS_ID", "Sex", "Sex_PLINK", "Age",
                    "BMI_2010", "Phenotype_PLINK",
                    "EUR", "AFR", "NAM", "EAS"]
    meta.to_csv(OUT_META, index=False)
    print(f"  ✓ {OUT_META}  (full merged dataset)")

    # ── Cohort summary ────────────────────────────────────────────────────────
    print("\n" + "=" * 55)
    print(" COHORT SUMMARY")
    print("=" * 55)
    n_valid = n_cases + n_controls
    print(f"  Total in dataset         : {len(df)}")
    print(f"  With BMI-2010            : {n_valid} ({n_valid/len(df)*100:.1f}%)")
    print(f"  Missing BMI-2010 (→ -9)  : {n_missing}")
    print(f"  Cases (BMI ≥ 25)         : {n_cases} ({n_cases/n_valid*100:.1f}% of valid)")
    print(f"    Overweight (25–29.9)   : {((df['BMI_2010']>=25)&(df['BMI_2010']<30)).sum()}")
    print(f"    Obese (≥ 30)           : {(df['BMI_2010']>=30).sum()}")
    print(f"  Controls (BMI < 25)      : {n_controls} ({n_controls/n_valid*100:.1f}% of valid)")

    for code, label in [(1, "Controls"), (2, "Cases")]:
        sub = df[df["PHENOTYPE"] == code]
        bmi = sub["BMI_2010"]
        age = pd.to_numeric(sub["Age-at-blood-collection"], errors="coerce")
        print(f"\n  {label} (n={len(sub)}):")
        print(f"    BMI  : {bmi.mean():.1f} ± {bmi.std():.1f}  [{bmi.min():.1f}–{bmi.max():.1f}]")
        print(f"    Age  : {age.mean():.1f} ± {age.std():.1f} years")
        print(f"    Sex  : {(sub['Gender']=='M').sum()}M ({(sub['Gender']=='M').mean()*100:.0f}%) / "
              f"{(sub['Gender']=='F').sum()}F ({(sub['Gender']=='F').mean()*100:.0f}%)")
        print(f"    EUR  : {sub['EUR'].mean():.3f}  "
              f"AFR: {sub['AFR'].mean():.3f}  "
              f"NAM: {sub['NAM'].mean():.3f}")

    print("\nDone. Files ready for PLINK pipeline.")
    print("Next: bash scripts/01_vcf_to_plink.sh")


if __name__ == "__main__":
    main()
