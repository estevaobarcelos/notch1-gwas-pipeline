# NOTCH1 GWAS Pipeline

Nov 13, 2026 · @EB

## Overview

Reproducible pipeline for the genetic association study of NOTCH1 SNPs and excessive weight in the SABE cohort (São Paulo, Brazil). Accompanies the published study:

> Silva Barcelos, E. et al. (2024). Genetic variation in NOTCH1 is associated with overweight and obesity in Brazilian elderly. *Scientific Reports*, 14, 170963. https://doi.org/10.1038/s41598-024-65771-1

**Key finding:** SNP rs9411207 (NOTCH1, intron, chr9:136,514,099) is associated with excessive weight under the log-additive model, with the TT genotype conferring increased risk (OR 1.50, 95% CI 1.20–1.88; p = 0.0002).

| Parameter | Value |
| --- | --- |
| Study | SABE (Saúde, Bem-Estar e Envelhecimento) |
| Location | São Paulo, Brazil |
| Design | Population-based cohort |
| Total (WGS) | 1,171 individuals |
| With BMI-2010 | 1,024 individuals |
| Cases (BMI ≥ 25) | 744 (72.7%) |
| Controls (BMI < 25) | 280 (27.3%) |
| Sex | 396 M / 707 F (68 ambiguous) |
| Reference genome | GRCh38 (hg38) |
| Genotyping | Whole-genome sequencing (WGS) |
| Ancestry | Admixed Brazilian (EUR + AFR + NAM + EAS) |

## Background

Excessive weight (BMI ≥ 25) is highly prevalent and associated with cardiovascular disease, type-2 diabetes, and other chronic disorders. NOTCH1 is a key regulator of stem cell maintenance, adipogenesis, and metabolism. This pipeline tests whether NOTCH1 SNPs are associated with excessive weight in an elderly Brazilian cohort using whole-genome sequencing data from the SABE study.

## Key Results

### rs9411207 — All Models

| Model | OR | 95% CI | p-value |
| --- | --- | --- | --- |
| Unadjusted (log-additive) | 1.328 | — | 0.0049 |
| + Sex | 1.332 | 1.092–1.624 | 0.0047 |
| Model 1: + Age, Sex, Ancestry | 1.474 | 1.192–1.824 | 0.0003 |
| Model 2: + HDL, TG, hsCRP, HbA1c, Hypertension | 1.496 | 1.200–1.866 | 0.0003 |
| Published (Errera et al. 2024) | 1.500 | 1.200–1.880 | 0.0002 |

### rs9411207 — Genotype Distribution

| Genotype | Cases n (%) | Controls n (%) |
| --- | --- | --- |
| CC (reference) | 245 (33.4%) | 116 (41.7%) |
| TC | 342 (46.7%) | 123 (44.2%) |
| TT (risk) | 146 (19.9%) | 39 (14.0%) |

T allele frequency: cases = 0.4325 / controls = 0.3615

## Pipeline

1. **Prepare phenotype files** — BMI-2010 → case (≥25) / control (<25); continuous phenotype; merge ancestry proportions (EUR, AFR, NAM; EAS excluded — collinear); build Model 1 and Model 2 covariate files.
2. **VCF → PLINK** — Convert raw VCF (WGS, VQSR-filtered, biallelic, hg38) to `.bed/.bim/.fam`; 4,103 variants / 1,171 samples.
3. **Quality control** — Call rate SNP > 95% (−40 variants); call rate sample > 98% (−0 samples); MAF > 1% (−3,491 variants); HWE p < 1×10⁻⁶ applied to controls only (−5 variants); heterozygosity ±3 SD (−0 samples).
4. **Apply phenotype + covariates** — 744 cases / 280 controls / 147 missing.
5. **Logistic association analysis** — Allelic chi-square; all genetic models (GENO, TREND, ALLELIC, DOM, REC); logistic additive unadjusted; Model 1 (+ age, sex, EUR, AFR, NAM); Model 2 (+ HDL, TG, hsCRP, HbA1c, hypertension) — primary.
6. **Multiple testing correction** — Benjamini–Hochberg FDR (q < 0.05 and q < 0.10); permutation test (10,000 permutations, empirical p-values).
7. **Linear regression** — BMI-2010 (continuous) \~ genotype; all SNPs β/SE/t/p/BH q; rs9411207 unadjusted + adjusted; genotype-level ANOVA (CC vs TC vs TT); pairwise t-tests; allele dosage trend (Spearman ρ).

**Result:** rs9411207 replicated OR = 1.496 (1.200–1.866), p = 0.0003.

## Repository Structure

```
notch1-gwas-pipeline/
├── README.md
├── LICENSE
├── run_full_pipeline.sh
├── .github/workflows/          # CI workflow
├── 00_prepare_phenotype_files.py
├── 01_vcf_to_plink.sh
├── 02_qc_pipeline.sh
├── 03_heterozygosity_check.py
├── 04_pca_analysis.R
├── 05_association_analysis.sh
├── 06_results_summary.py
├── 07_multiple_testing_correction.py
└── 08_linear_regression_bmi.py
```

| Script | Purpose |
| --- | --- |
| `00_prepare_phenotype_files.py` | Phenotype + covariate files (BMI binary + continuous) |
| `01_vcf_to_plink.sh` | VCF → PLINK conversion |
| `02_qc_pipeline.sh` | Quality control pipeline |
| `03_heterozygosity_check.py` | Outlier detection + plot |
| `04_pca_analysis.R` | PCA + ancestry visualization |
| `05_association_analysis.sh` | Logistic regression — Models 1 & 2 |
| `06_results_summary.py` | QC summary table |
| `07_multiple_testing_correction.py` | BH FDR + permutation (10,000×) |
| `08_linear_regression_bmi.py` | BMI continuous \~ genotype (OLS) |
| `run_full_pipeline.sh` | Runs all stages in order |

## Usage

### Requirements

```bash
# System
sudo apt install plink1.9

# Python
pip install pandas numpy scipy statsmodels matplotlib seaborn openpyxl

# R
install.packages(c("ggplot2", "dplyr", "ggrepel"))
```

### Run

```bash
git clone https://github.com/estevaobarcelosestevaobarcelos/notch1-gwas-pipeline.git
cd notch1-gwas-pipeline

# 1. Prepare phenotype and covariate files
python3 00_prepare_phenotype_files.py

# 2. Convert VCF to PLINK
bash 01_vcf_to_plink.sh

# 3. Quality control
bash 02_qc_pipeline.sh
python3 03_heterozygosity_check.py
Rscript 04_pca_analysis.R

# 4. Association analysis
bash 05_association_analysis.sh

# 5. Results summary + multiple testing
python3 06_results_summary.py
python3 07_multiple_testing_correction.py
python3 08_linear_regression_bmi.py

# Or run everything at once
bash run_full_pipeline.sh
```

### CI note

The GitHub Actions workflow runs syntax checks on the Python scripts. It requires `actions/setup-python@v5` (not v3) and Python 3.10 or 3.11; a `tests/` directory with at least one test file is needed to avoid pytest exit code 5 (no tests collected).

## QC Summary

Overall genotyping rate after QC: **99.3%**

| Step | Variants | Samples | Removed |
| --- | --- | --- | --- |
| Raw VCF | 4,103 | 1,171 | — |
| Call rate SNP > 95% | 4,063 | 1,171 | 40 SNPs |
| Call rate sample > 98% | 4,063 | 1,171 | 0 samples |
| MAF > 1% | 572 | 1,171 | 3,491 SNPs |
| HWE p < 1×10⁻⁶ (controls only) | 567 | 1,171 | 5 SNPs |
| Heterozygosity ±3 SD | 567 | 1,171 | 0 samples |
| **Final QC dataset** | **567** | **1,171** | — |

## Multiple Testing Correction

Three approaches are applied and compared. BH FDR and permutation are preferred over Bonferroni because SNPs in a 50-kb region are correlated by LD — Bonferroni assumes 567 independent tests when the effective number is much smaller.

| Method | Threshold | rs9411207 survives? | Rationale |
| --- | --- | --- | --- |
| Bonferroni | p < 8.8×10⁻⁵ (0.05/567) | ✗ No | Assumes independence — overcorrects for correlated SNPs |
| BH FDR | q < 0.05 | ✓ Yes | Accounts for correlation — preferred for regional studies |
| BH FDR | q < 0.10 | ✓ Yes | — |
| Permutation | p\_emp < 0.05 | ✓ Yes | Gold standard — accounts for LD structure automatically |

Script `07_multiple_testing_correction.py` runs BH FDR and 10,000 permutations, produces the null distribution figure, and saves all q-values.

## Linear Regression — BMI as Continuous Outcome

Script `08_linear_regression_bmi.py` complements the logistic (case/control) analysis by treating BMI as a continuous variable.

**Question:** Does the T allele of rs9411207 predict a higher BMI value, beyond just defining overweight vs normal weight?

**Models:**

- Unadjusted: BMI \~ genotype\_dosage
- Adjusted: BMI \~ genotype\_dosage + age + sex + EUR + AFR + NAM + HDL + TG + hsCRP + HbA1c + hypertension

**Additional tests for rs9411207:**

- One-way ANOVA across CC / TC / TT genotype groups
- Pairwise t-tests (Bonferroni-corrected within 3 pairs)
- Allele dosage trend test (Spearman ρ)

**Outputs:** β coefficient (kg/m² per T allele), SE, t, p, BH q-value for all 567 SNPs; violin + boxplot + regression figure for rs9411207.

## Methodological Notes

**HWE applied to controls only.** In case-control studies, applying HWE filtering to cases risks discarding disease-associated SNPs that show HWE departure due to selection. We apply HWE (p < 1×10⁻⁶) exclusively to controls, following Marees et al. (2018).

**HWE limitations in admixed populations.** In the admixed SABE cohort, HWE can be violated due to the Wahlund effect (mixing of subpopulations with different allele frequencies), independently of genotyping quality. HWE is used as a supplementary filter only. Call rate and differential missingness are the primary QC criteria.

**Collinearity in ancestry covariates.** EUR + AFR + NAM + EAS sum to exactly 1.0 for each individual (perfect collinearity). EAS was excluded from regression models; EUR, AFR, and NAM were retained as covariates.

**PCA limitation.** Only 98 SNPs were available for PCA after LD pruning (focused NOTCH1 panel). Ancestry correction is based on ADMIXTURE-derived proportions (EUR, AFR, NAM) included as continuous covariates in Model 1 and Model 2.

## Citation

If you use this pipeline, please cite the published study:

```bibtex
@article{barcelos2024notch1,
  title   = {Genetic variation in NOTCH1 is associated with overweight and
             obesity in Brazilian elderly},
  author  = {Silva Barcelos, Estev{\~a}o and
             Naslavsky, Michel Satya and
             Fernandes, Izadora Silveira and
             Scliar, Marilia Oliveira and
             Yamamoto, Guilherme Lopes and
             Wang, Jaqueline Yu Ting and
             Bride, La{\'{i}}s and
             de Sousa, Valdemir Pereira and
             Pimassoni, Lucia Helena Sagrillo and
             Sportoletti, Paolo and
             de Paula, Flavia and
             von Zeidler, Sandra Ventorin and
             Duarte, Yeda Aparecida Oliveira and
             Passos-Bueno, Maria Rita and
             Zatz, Mayana and
             Errera, Fl{\'{a}}via Imbroisi Valle},
  journal = {Scientific Reports},
  volume  = {14},
  pages   = {17096},
  year    = {2024},
  doi     = {10.1038/s41598-024-65771-1},
  url     = {https://doi.org/10.1038/s41598-024-65771-1}
}
```

> **Note:** Retrieve it at [doi.org/10.1038/s41598-024-65771-1](https://doi.org/10.1038/s41598-024-65771-1)

## References

1. Marees AT et al. A tutorial on conducting genome-wide association studies. *Int J Methods Psychiatr Res.* 2018. PMID: 29484742
2. Anderson CA et al. Data quality control in genetic case-control association studies. *Nat Protoc.* 2010.
3. Wigginton JE et al. A note on exact tests of Hardy–Weinberg equilibrium. *Am J Hum Genet.* 2005.
4. Purcell S et al. PLINK: A tool set for whole-genome association and population-based linkage analyses. *Am J Hum Genet.* 2007.

## License

MIT — see LICENSE

## Author

Estevão Barcelos · [GitHub](https://github.com/estevaobarcelosestevaobarcelos)
