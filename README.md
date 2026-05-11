# LRS Diagnostic Yield Meta-analysis — Analysis Scripts

R scripts and figure code accompanying the systematic review and IPD meta-analysis of long-read sequencing (LRS) diagnostic yield in short-read-negative rare disease patients.

## Contents

| File | Purpose |
|---|---|
| `Step1_Overall_Rescue_Rates_submission.R` | Crude rescue rates, Wilson CIs, internal consistency checks (Supplementary Table S4) |
| `Step2_Prior_Testing_v2_submission.R` | Rescue by prior testing modality, Fisher's exact, Cochran-Armitage trend test (Supplementary Table S8). Sourced by the Figure S2 script |
| `Step3_HPO_generate_corrected_excel_submission.R` | HPO ontology DAG traversal, bottom-up category aggregation, Fisher's exact + BH FDR per category, within-neurology sensitivity analysis. Produces `patient_hpo_categories.xlsx` used by Figure 2b |
| `Step4_Rescue_Mechanisms_submission.R` | Rescue mechanism distribution and cross-tabulations (underlies Figure 3) |
| `Step6_Variant_Characteristics_submission.R` | Variant sizes, types, genomic regions, ACMG classifications in rescued cases |
| `Step7_Meta_Analysis_submission.R` | Random-effects meta-analysis (Freeman-Tukey arcsine, REML), leave-one-out, Egger's test, subgroup analyses (Tables S5, S6, S7, S9) |
| `Step8_excluded_studies_submission.R` | Descriptive analysis of Tier 1 studies excluded from the primary meta-analysis (Tables S10, S11) |
| `Figure2a_ForestPlot_CORRECTED.R` | Forest plot for Group A primary meta-analysis (Figure 2a). Values are taken from Step 7 output |
| `Figure2b_ForestPlot_CORRECTED.R` | HPO phenotype forest plot (Figure 2b). Reads `patient_hpo_categories.xlsx` from Step 3 |

## Input data

All scripts read from `Dataset_S1.xlsx` (supplied with the manuscript). Step 3 additionally requires `hp.json` (HPO ontology, OBO Graph JSON format), downloadable from https://hpo.jax.org/data/ontology.

Place `Dataset_S1.xlsx` and `hp.json` in the same directory as the scripts before running.

## Run order

Steps are independent except for the following dependencies:

- Step 3 must be run before Figure 2b (Figure 2b reads `patient_hpo_categories.xlsx` produced by Step 3)
- Step 7 should be run before Figure 2a if you want to verify the hardcoded summary values against your own pipeline output

Otherwise the steps can be run in any order. A natural sequence is Step 1 → Step 2 → Step 3 → Step 4 → Step 6 → Step 7 → Step 8 → Figure 2a → Figure 2b.

## R version and packages

Developed on R 4.4. Required packages (the scripts will install any missing ones automatically):

- `readxl`, `dplyr`, `tidyr` (used throughout)
- `jsonlite`, `openxlsx` (Step 3)
- `meta`, `metafor` (Step 7)
- `DescTools` (Step 2; optional, used only for the Cochran-Armitage Z-statistic cross-check)
- `forestploter`, `grid` (Figures 2a, 2b)

## Output

Each Step script writes a plain-text log to a file named `Step{N}_..._OUTPUT.txt` in the working directory. The Figure scripts write PDF and PNG (Figure 2a) or PDF (Figure 2b) to the script directory.

## Notes

- Dataset_S1 is read from the sheet `Tier 1 - Group A` (11 studies, N=538) except for Step 8, which reads `Tier 1 - Group B+C`.
- The first two rows of each sheet are a group banner and the actual column headers; the scripts skip the banner row.
- Random seeds are not used; all analyses are deterministic.
