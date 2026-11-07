#!/usr/bin/env bash
# =============================================================================
# 05_association_analysis.sh
# Association analysis for NOTCH1 SNPs — excessive weight (BMI >= 25)
#
# Models:
#   A. Allelic chi-square (--assoc)
#   B. All genetic models for rs9411207 (--model)
#   C. Logistic additive, unadjusted (--logistic)
#   D. Model 1: + age, sex, ancestry (EUR AFR NAM)
#   E. Model 2: + HDL, TG, hsCRP, HbA1c, hypertension  ← primary
#
# Input  : notch1_phenotyped.{bed,bim,fam}
#          sample_covariates_model1.txt
#          sample_covariates_model2.txt
#
# Key result: rs9411207 OR=1.496 (95% CI 1.200-1.866), p=0.0003 (Model 2)
#             Replicates published OR=1.500 (1.200-1.880), p=0.0002
#
# Author: Estevão Barcelos
# Paper : Errera et al. Scientific Reports 2024
#         https://doi.org/10.1038/s41598-024-65771-1
# =============================================================================

set -euo pipefail
PLINK=$(command -v plink1.9 2>/dev/null || command -v plink)
IN="notch1_phenotyped"

echo "=============================================="
echo " Step 05: Association Analysis"
echo " Outcome: excessive weight (BMI-2010 >= 25)"
echo " Cases=744  Controls=280"
echo "=============================================="

# ── A. Allelic chi-square ─────────────────────────────────────────────────────
echo ""
echo "[A] Allelic chi-square test (unadjusted)"
$PLINK \
    --bfile "$IN" \
    --allow-no-sex \
    --assoc \
    --out notch1_assoc_allelic

echo "    Top 10 results:"
sort -k9 -g notch1_assoc_allelic.assoc | head -11

# ── B. All genetic models for key SNP ────────────────────────────────────────
echo ""
echo "[B] All genetic models — rs9411207"
$PLINK \
    --bfile "$IN" \
    --allow-no-sex \
    --model \
    --snp rs9411207 \
    --out rs9411207_all_models

cat rs9411207_all_models.model

# ── C. Logistic additive, unadjusted ─────────────────────────────────────────
echo ""
echo "[C] Logistic regression — additive, unadjusted"
$PLINK \
    --bfile "$IN" \
    --allow-no-sex \
    --logistic \
    --ci 0.95 \
    --hide-covar \
    --out notch1_logistic_unadj

echo "    Top 10:"
sort -k9 -g notch1_logistic_unadj.assoc.logistic | grep -v "NA" | head -11
echo "    rs9411207:"
grep "rs9411207" notch1_logistic_unadj.assoc.logistic

# ── D. Model 1: age + sex + ancestry ─────────────────────────────────────────
echo ""
echo "[D] Model 1 — logistic + age + sex + EUR + AFR + NAM"
if [[ ! -f "sample_covariates_model1.txt" ]]; then
    echo "    ERROR: sample_covariates_model1.txt not found."
    echo "    Run: python3 scripts/00_prepare_phenotype_files.py"
    exit 1
fi

$PLINK \
    --bfile "$IN" \
    --allow-no-sex \
    --logistic \
    --covar sample_covariates_model1.txt \
    --ci 0.95 \
    --hide-covar \
    --out notch1_model1

echo "    rs9411207 (Model 1):"
grep "rs9411207" notch1_model1.assoc.logistic

# ── E. Model 2: age + sex + ancestry + clinical variables ────────────────────
echo ""
echo "[E] Model 2 — logistic + age + sex + ancestry + HDL + TG + hsCRP + HbA1c + hypertension"
echo "    NOTE: ancestry covariates = EUR + AFR + NAM (EAS excluded — collinear)"
if [[ ! -f "sample_covariates_model2.txt" ]]; then
    echo "    ERROR: sample_covariates_model2.txt not found."
    echo "    Run: python3 scripts/00_prepare_phenotype_files.py"
    exit 1
fi

$PLINK \
    --bfile "$IN" \
    --allow-no-sex \
    --logistic \
    --covar sample_covariates_model2.txt \
    --ci 0.95 \
    --hide-covar \
    --out notch1_model2

echo "    rs9411207 (Model 2 — PRIMARY):"
grep "rs9411207" notch1_model2.assoc.logistic

echo "    Top 10 (Model 2):"
sort -k9 -g notch1_model2.assoc.logistic | grep -v "NA" | head -11

# ── Summary ──────────────────────────────────────────────────────────────────
echo ""
echo "=============================================="
echo " rs9411207 — RESULTS SUMMARY"
echo "=============================================="
echo ""
printf "%-40s %6s  %-18s  %s\n" "Model" "OR" "95% CI" "p-value"
printf "%-40s %6s  %-18s  %s\n" "---" "---" "---" "---"

extract_rs() {
    local file="$1"
    grep "rs9411207" "$file" | awk '{printf "%-40s %6.3f  %s – %s  %s\n", "'"$2"'", $7, $9, $10, $11}'
}

extract_rs "notch1_logistic_unadj.assoc.logistic" "Unadjusted (additive)"
extract_rs "notch1_model1.assoc.logistic"          "Model 1 (age+sex+ancestry)"
extract_rs "notch1_model2.assoc.logistic"          "Model 2 (+ clinical vars) ★"
echo "Published (Errera et al. 2024)              1.500   1.200 – 1.880    0.0002  ◆"
echo ""
echo "★ Primary model  ◆ Published result"
echo ""
echo "Step 05 complete."
echo "Next: python3 scripts/06_results_summary.py"

# ── F. PLINK linear regression — BMI as continuous outcome ───────────────────
echo ""
echo "[F] PLINK linear regression — BMI-2010 as continuous outcome"
echo "    NOTE: This runs genome-wide linear regression directly in PLINK."
echo "          Script 08 (python) provides the full analysis with plots."

# Need a quantitative phenotype file with BMI values
if [[ -f "sample_bmi_quantitative.txt" ]]; then
    $PLINK \
        --bfile notch1_phenotyped \
        --allow-no-sex \
        --linear \
        --pheno sample_bmi_quantitative.txt \
        --covar sample_covariates_model2.txt \
        --hide-covar \
        --ci 0.95 \
        --out notch1_linear_bmi

    echo "    Top 10 — BMI linear regression:"
    sort -k9 -g notch1_linear_bmi.assoc.linear | grep -v "NA" | head -11

    echo "    rs9411207:"
    grep "rs9411207" notch1_linear_bmi.assoc.linear
else
    echo "    sample_bmi_quantitative.txt not found."
    echo "    Run python3 scripts/00_prepare_phenotype_files.py to generate it."
fi
