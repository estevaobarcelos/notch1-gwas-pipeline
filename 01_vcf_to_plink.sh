#!/usr/bin/env bash
# =============================================================================
# 01_vcf_to_plink.sh
# Convert annotated VCF to PLINK binary format (.bed/.bim/.fam)
#
# Input : notch1_original.vcf  (WGS, VQSR-filtered, biallelic, hg38)
# Output: notch1_plink.{bed,bim,fam,log}
#
# Cohort: SABE (n=1171, admixed Brazilian, chr9 NOTCH1 ±50kb region)
# Author: Estevão Barcelos
# =============================================================================

set -euo pipefail

VCF="notch1_original.vcf"
OUT="notch1_plink"

echo "=============================================="
echo " Step 01: VCF → PLINK conversion"
echo "=============================================="
echo "Input VCF : $VCF"
echo "Output    : $OUT.*"
echo ""

# Check input exists
if [[ ! -f "$VCF" ]]; then
    echo "ERROR: $VCF not found. Place your VCF in the project root."
    exit 1
fi

# Check PLINK is available
if ! command -v plink1.9 &>/dev/null && ! command -v plink &>/dev/null; then
    echo "ERROR: PLINK not found. Install with: sudo apt install plink1.9"
    exit 1
fi

PLINK=$(command -v plink1.9 2>/dev/null || command -v plink)

$PLINK \
    --vcf "$VCF" \
    --double-id \
    --vcf-half-call missing \
    --make-bed \
    --out "$OUT"

echo ""
echo "Output files:"
ls -lh ${OUT}.{bed,bim,fam,log} 2>/dev/null

echo ""
echo "Variant and sample counts:"
echo "  Variants : $(wc -l < ${OUT}.bim)"
echo "  Samples  : $(wc -l < ${OUT}.fam)"
echo ""
echo "Step 01 complete."
