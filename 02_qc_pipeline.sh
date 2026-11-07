#!/usr/bin/env bash
# =============================================================================
# 02_qc_pipeline.sh
# Full quality control pipeline for NOTCH1 region WGS data
#
# Steps:
#   2a. SNP call rate filter  (--geno 0.05)
#   2b. Sample call rate filter (--mind 0.02)
#   2c. MAF filter            (--maf 0.01)
#   2d. HWE filter            (--hwe 1e-6, controls only when phenotype available)
#   2e. LD pruning            (for PCA input)
#   2f. PCA                   (10 principal components)
#
# IMPORTANT NOTE ON HWE IN ADMIXED POPULATIONS:
#   The SABE cohort is highly admixed (European + African + Amerindian ancestry).
#   HWE can be violated in admixed samples due to population stratification
#   (Wahlund effect), independently of genotyping quality. The HWE filter is
#   therefore applied with a conservative threshold (1e-6) and treated as a
#   supplementary filter, with call rate and differential missingness as the
#   primary QC metrics.
#
# Author: Estevão Barcelos
# Reference: Marees et al. 2018 (PMID: 29484742)
# =============================================================================

set -euo pipefail

PLINK=$(command -v plink1.9 2>/dev/null || command -v plink)
IN="notch1_plink"

echo "=============================================="
echo " Step 02: Quality Control Pipeline"
echo "=============================================="

# ── 2a. SNP call rate ─────────────────────────────────────────────────────────
echo ""
echo "[2a] SNP call rate filter (geno > 5% missing excluded)"
$PLINK \
    --bfile "$IN" \
    --geno 0.05 \
    --make-bed \
    --out notch1_cr_snp

REMOVED_GENO=$(( $(wc -l < ${IN}.bim) - $(wc -l < notch1_cr_snp.bim) ))
echo "     Variants removed: $REMOVED_GENO"

# ── 2b. Sample call rate ──────────────────────────────────────────────────────
echo ""
echo "[2b] Sample call rate filter (mind > 2% missing excluded)"
$PLINK \
    --bfile notch1_cr_snp \
    --mind 0.02 \
    --make-bed \
    --out notch1_cr_ind

REMOVED_MIND=$(( $(wc -l < notch1_cr_snp.fam) - $(wc -l < notch1_cr_ind.fam) ))
echo "     Samples removed: $REMOVED_MIND"

# ── 2c. MAF filter ────────────────────────────────────────────────────────────
echo ""
echo "[2c] Minor allele frequency filter (MAF < 1% excluded)"
echo "     NOTE: High removal expected — NOTCH1 region is rich in rare variants."
echo "           Rare variants (MAF < 1%) should be analyzed separately with"
echo "           burden tests (SKAT, SKAT-O) rather than single-variant tests."
$PLINK \
    --bfile notch1_cr_ind \
    --maf 0.01 \
    --make-bed \
    --out notch1_maf

REMOVED_MAF=$(( $(wc -l < notch1_cr_ind.bim) - $(wc -l < notch1_maf.bim) ))
echo "     Variants removed: $REMOVED_MAF"

# ── 2d. HWE filter ────────────────────────────────────────────────────────────
echo ""
echo "[2d] Hardy-Weinberg equilibrium filter (p < 1e-6)"
echo "     NOTE: Applied to full sample (no phenotype file provided)."
echo "           When phenotype is available, re-run with --filter-controls"
echo "           to apply HWE only to controls. See docs/methods.md."
$PLINK \
    --bfile notch1_maf \
    --hwe 1e-6 \
    --make-bed \
    --out notch1_hwe

REMOVED_HWE=$(( $(wc -l < notch1_maf.bim) - $(wc -l < notch1_hwe.bim) ))
echo "     Variants removed: $REMOVED_HWE"

# ── 2e. Heterozygosity report (for outlier detection in script 03) ────────────
echo ""
echo "[2e] Generating heterozygosity report"
$PLINK \
    --bfile notch1_hwe \
    --het \
    --out notch1_het
echo "     Report written to notch1_het.het"

# ── 2f. LD pruning for PCA ───────────────────────────────────────────────────
echo ""
echo "[2f] LD pruning for PCA input (window=50, step=10, r2=0.1)"
$PLINK \
    --bfile notch1_hwe \
    --indep-pairwise 50 10 0.1 \
    --out notch1_pruned

PRUNED_IN=$(wc -l < notch1_pruned.prune.in)
PRUNED_OUT=$(wc -l < notch1_pruned.prune.out)
echo "     SNPs retained : $PRUNED_IN"
echo "     SNPs pruned   : $PRUNED_OUT"
echo ""
echo "     WARNING: Only $PRUNED_IN SNPs available for PCA after pruning."
echo "     This is insufficient for robust ancestry inference (typically >100k)."
echo "     PCA results should be treated as exploratory only."

# ── 2g. PCA ──────────────────────────────────────────────────────────────────
echo ""
echo "[2g] Running PCA (10 components)"
$PLINK \
    --bfile notch1_hwe \
    --extract notch1_pruned.prune.in \
    --pca 10 \
    --out notch1_pca
echo "     Results: notch1_pca.eigenval + notch1_pca.eigenvec"

# ── Summary ──────────────────────────────────────────────────────────────────
echo ""
echo "=============================================="
echo " QC SUMMARY"
echo "=============================================="
printf "%-35s %10s %10s\n" "Step" "Variants" "Samples"
printf "%-35s %10s %10s\n" "---" "---" "---"
printf "%-35s %10s %10s\n" "Raw input"           "$(wc -l < ${IN}.bim)"           "$(wc -l < ${IN}.fam)"
printf "%-35s %10s %10s\n" "Call rate SNP (>95%)" "$(wc -l < notch1_cr_snp.bim)"  "$(wc -l < notch1_cr_snp.fam)"
printf "%-35s %10s %10s\n" "Call rate sample (>98%)" "$(wc -l < notch1_cr_ind.bim)" "$(wc -l < notch1_cr_ind.fam)"
printf "%-35s %10s %10s\n" "MAF > 1%"            "$(wc -l < notch1_maf.bim)"     "$(wc -l < notch1_maf.fam)"
printf "%-35s %10s %10s\n" "HWE p < 1e-6"        "$(wc -l < notch1_hwe.bim)"     "$(wc -l < notch1_hwe.fam)"
echo ""
echo "Step 02 complete. Run script 03 for heterozygosity outlier detection."
