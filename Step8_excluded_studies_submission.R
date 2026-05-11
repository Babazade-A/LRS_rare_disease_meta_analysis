###############################################################################
#  STEP 8 — EXCLUDED STUDIES: DESCRIPTIVE ANALYSES (R VERIFICATION SCRIPT)
#  Systematic Review & IPD Meta-Analysis: LRS Diagnostic Yield in Rare Disease
#  All Group B + Group C studies excluded from the primary Group A analysis
#  Output: step8_excluded_studies_output.txt
###############################################################################

# ── Libraries ----------------------------------------------------------------
library(readxl)
library(dplyr)
library(tidyr)

# ── Wilson CI helper ----------------------------------------------------------
wilson_ci <- function(x, n, digits = 1) {
  if (n == 0) return("0/0 (NA)")
  p <- x / n
  z <- qnorm(0.975)
  denom <- 1 + z^2 / n
  mid   <- (p + z^2 / (2 * n)) / denom
  half  <- (z / denom) * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  lo <- max(0, mid - half) * 100
  hi <- min(1, mid + half) * 100
  pct <- p * 100
  sprintf("%d/%d (%.1f%%; 95%% CI %.1f–%.1f%%)", x, n, pct, lo, hi)
}

# ── Read data -----------------------------------------------------------------
dat <- read_excel("Dataset_S1.xlsx", sheet = "Tier 1 - Group B+C", skip = 1)

# Standardise column names
names(dat)[names(dat) == "Previous_SR_Result (or simultaneous)"] <- "SR_Result"

# 8 mechanism columns
mech_names <- c(
  "Reanalysis"   = grep("^Reanalysis_finding", names(dat), value = TRUE),
  "SV_Detection" = grep("^SV_Detection",       names(dat), value = TRUE),
  "Dark_Region"  = grep("^Dark_Region",        names(dat), value = TRUE),
  "Repeat_Exp"   = grep("^Repeat_Expansion",   names(dat), value = TRUE),
  "Phasing"      = grep("^Phasing",            names(dat), value = TRUE),
  "FullLength"   = grep("^Full_Length",         names(dat), value = TRUE),
  "Methylation"  = grep("^Methylation",        names(dat), value = TRUE),
  "Mosaicism"    = grep("^Mosaicism",          names(dat), value = TRUE)
)

# helper: count non-"No" non-NA entries in a mechanism column for a subset
count_mech <- function(df, mech_col) {
  vals <- df[[mech_col]]
  vals <- vals[!is.na(vals) & vals != "No"]
  if (length(vals) == 0) return(NULL)
  tbl <- sort(table(vals), decreasing = TRUE)
  data.frame(Value = names(tbl), N = as.integer(tbl), stringsAsFactors = FALSE)
}

# helper: rescued subset
is_rescued <- function(outcome) {
  outcome %in% c("Definitive_Rescue", "Possible_Rescue")
}

###############################################################################
sink("step8_excluded_studies_output.txt")
cat(strrep("=", 80), "\n")
cat("STEP 8 — EXCLUDED STUDIES: DESCRIPTIVE ANALYSES\n")
cat("Systematic Review & IPD Meta-Analysis of LRS in Rare Disease\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(strrep("=", 80), "\n\n")

cat("Total records in Dataset_S1.xlsx (Tier 1 - Group B+C):", nrow(dat), "\n")
cat("Studies present:", paste(sort(unique(dat$Study_ID)), collapse = ", "), "\n\n")

###############################################################################
# 8a. PARALLEL / FIRST-LINE STUDIES
###############################################################################
cat(strrep("=", 80), "\n")
cat("8a. PARALLEL / FIRST-LINE STUDIES\n")
cat("    (Ek_2025, Smits_2025, Kamolvisit_2025)\n")
cat(strrep("=", 80), "\n\n")

parallel_ids <- c("Ek_2025", "Smits_2025", "Kamolvisit_2025")

for (sid in parallel_ids) {
  sub <- dat %>% filter(Study_ID == sid)
  n_total <- nrow(sub)
  
  cat(strrep("-", 60), "\n")
  cat(sid, " (N =", n_total, ")\n")
  cat(strrep("-", 60), "\n\n")
  
  # SR Result distribution
  cat("  Previous SR Result distribution:\n")
  sr_tbl <- table(sub$SR_Result, useNA = "ifany")
  for (nm in names(sr_tbl)) {
    cat(sprintf("    %-25s %d (%s)\n", nm, sr_tbl[nm],
                sprintf("%.1f%%", sr_tbl[nm]/n_total*100)))
  }
  
  # LRS Outcome distribution
  cat("\n  LRS Outcome distribution:\n")
  lrs_tbl <- table(sub$LRS_Outcome, useNA = "ifany")
  for (nm in names(lrs_tbl)) {
    cat(sprintf("    %-25s %d (%s)\n", nm, lrs_tbl[nm],
                sprintf("%.1f%%", lrs_tbl[nm]/n_total*100)))
  }
  
  # Concordance (SR-Positive that are Concordant by LRS)
  sr_pos <- sub %>% filter(SR_Result == "Positive")
  n_pos  <- nrow(sr_pos)
  n_conc <- sum(sr_pos$LRS_Outcome == "Concordant", na.rm = TRUE)
  if (n_pos > 0) {
    cat(sprintf("\n  Concordance (SR+ confirmed by LRS): %s\n",
                wilson_ci(n_conc, n_pos)))
  } else {
    cat("\n  No SR-Positive cases (first-line LRS study)\n")
  }
  
  # LRS-added-value: Definitive_Rescue + Possible_Rescue + Refined
  lrs_added <- sub %>% filter(LRS_Outcome %in% c("Definitive_Rescue",
                                                   "Possible_Rescue",
                                                   "Refined"))
  n_added <- nrow(lrs_added)
  cat(sprintf("\n  LRS added-value cases (Rescue + Refined): %s\n",
              wilson_ci(n_added, n_total)))
  
  # LRS-unique findings (Rescued only)
  rescued <- sub %>% filter(is_rescued(LRS_Outcome))
  n_resc  <- nrow(rescued)
  cat(sprintf("  LRS-unique rescue: %s\n", wilson_ci(n_resc, n_total)))
  
  # Refined cases
  refined <- sub %>% filter(LRS_Outcome == "Refined")
  if (nrow(refined) > 0) {
    cat(sprintf("\n  Refined cases (N = %d):\n", nrow(refined)))
  }
  
  # Missed cases
  missed <- sub %>% filter(grepl("Missed", LRS_Outcome))
  if (nrow(missed) > 0) {
    cat(sprintf("\n  Missed by LRS (N = %d):\n", nrow(missed)))
  }
  
  # SR-Negative rescue rate
  sr_neg <- sub %>% filter(SR_Result == "Negative")
  if (nrow(sr_neg) > 0) {
    n_neg_resc <- sum(is_rescued(sr_neg$LRS_Outcome))
    cat(sprintf("\n  SR-Negative rescue rate: %s\n",
                wilson_ci(n_neg_resc, nrow(sr_neg))))
  }
  
  cat("\n")
}

###############################################################################
# 8b. SECOND-HIT / MISSING ALLELE GROUP
###############################################################################
cat(strrep("=", 80), "\n")
cat("8b. SECOND-HIT / MISSING ALLELE GROUP\n")
cat("    (All studies with SR_Result == 'Second_hit_missing')\n")
cat(strrep("=", 80), "\n\n")

sh <- dat %>% filter(SR_Result == "Second_hit_missing")
n_sh <- nrow(sh)

cat("Total N:", n_sh, "\n")
cat("Studies contributing:\n")
sh_study_tbl <- table(sh$Study_ID)
for (nm in names(sh_study_tbl)) {
  cat(sprintf("  %-20s %d\n", nm, sh_study_tbl[nm]))
}

cat(sprintf("\nLRS Outcome distribution:\n"))
sh_lrs <- table(sh$LRS_Outcome, useNA = "ifany")
for (nm in names(sh_lrs)) {
  cat(sprintf("  %-25s %d (%s)\n", nm, sh_lrs[nm],
              sprintf("%.1f%%", sh_lrs[nm]/n_sh*100)))
}

rescued_sh <- sh %>% filter(is_rescued(LRS_Outcome))
n_resc_sh  <- nrow(rescued_sh)
cat(sprintf("\nRescue rate: %s\n", wilson_ci(n_resc_sh, n_sh)))

# Definitive vs Possible
n_def_sh  <- sum(rescued_sh$LRS_Outcome == "Definitive_Rescue")
n_poss_sh <- sum(rescued_sh$LRS_Outcome == "Possible_Rescue")
cat(sprintf("  Definitive: %d | Possible: %d\n", n_def_sh, n_poss_sh))

cat("\nGenes in rescued cases:\n")
gene_tbl <- sort(table(rescued_sh$Gene), decreasing = TRUE)
for (nm in names(gene_tbl)) {
  cat(sprintf("  %-15s %d\n", nm, gene_tbl[nm]))
}

cat("\nVariant types in rescued cases:\n")
vt_tbl <- sort(table(rescued_sh$Rescued_Variant_Type), decreasing = TRUE)
for (nm in names(vt_tbl)) {
  cat(sprintf("  %-15s %d\n", nm, vt_tbl[nm]))
}

cat("\nMechanisms in rescued cases:\n")
for (mname in names(mech_names)) {
  mcol <- mech_names[[mname]]
  res <- count_mech(rescued_sh, mcol)
  if (!is.null(res)) {
    for (j in seq_len(nrow(res))) {
      cat(sprintf("  [%s] %s: %d\n", mname, res$Value[j], res$N[j]))
    }
  }
}

cat("\nInheritance mode:\n")
inh_tbl <- table(sh$Inheritance_Mode, useNA = "ifany")
for (nm in names(inh_tbl)) cat(sprintf("  %-30s %d\n", nm, inh_tbl[nm]))

cat("\nZygosity:\n")
zyg_tbl <- table(sh$Zygosity, useNA = "ifany")
for (nm in names(zyg_tbl)) cat(sprintf("  %-50s %d\n", nm, zyg_tbl[nm]))

cat("\n")

###############################################################################
# 8c. DMD / MUSCULAR DYSTROPHY (Bruels_2022 + Chu_2025)
###############################################################################
cat(strrep("=", 80), "\n")
cat("8c. DMD / MUSCULAR DYSTROPHY (Bruels_2022 + Chu_2025)\n")
cat(strrep("=", 80), "\n\n")

dmd <- dat %>% filter(Study_ID %in% c("Bruels_2022", "Chu_2025"))
n_dmd <- nrow(dmd)

cat("Total N:", n_dmd, "\n")
cat("By study:\n")
dmd_study <- table(dmd$Study_ID)
for (nm in names(dmd_study)) cat(sprintf("  %-15s %d\n", nm, dmd_study[nm]))

cat("\nSR Result distribution:\n")
dmd_sr <- table(dmd$SR_Result, useNA = "ifany")
for (nm in names(dmd_sr)) {
  cat(sprintf("  %-25s %d (%s)\n", nm, dmd_sr[nm],
              sprintf("%.1f%%", dmd_sr[nm]/n_dmd*100)))
}

cat("\nLRS Outcome distribution:\n")
dmd_lrs <- table(dmd$LRS_Outcome, useNA = "ifany")
for (nm in names(dmd_lrs)) {
  cat(sprintf("  %-25s %d (%s)\n", nm, dmd_lrs[nm],
              sprintf("%.1f%%", dmd_lrs[nm]/n_dmd*100)))
}

# SRS-negative rescue
dmd_neg <- dmd %>% filter(SR_Result == "Negative")
n_neg   <- nrow(dmd_neg)
n_neg_r <- sum(is_rescued(dmd_neg$LRS_Outcome))
cat(sprintf("\nSRS-Negative subset: N = %d\n", n_neg))
cat(sprintf("SRS-Negative rescue rate: %s\n", wilson_ci(n_neg_r, n_neg)))

# Combined rescue (all SR statuses)
rescued_dmd <- dmd %>% filter(is_rescued(LRS_Outcome))
n_resc_dmd  <- nrow(rescued_dmd)
cat(sprintf("\nOverall rescue rate (all cases): %s\n", wilson_ci(n_resc_dmd, n_dmd)))
cat(sprintf("  Definitive: %d | Possible: %d\n",
            sum(rescued_dmd$LRS_Outcome == "Definitive_Rescue"),
            sum(rescued_dmd$LRS_Outcome == "Possible_Rescue")))

cat("\nGenes:\n")
gene_dmd <- sort(table(rescued_dmd$Gene), decreasing = TRUE)
for (nm in names(gene_dmd)) cat(sprintf("  %-15s %d\n", nm, gene_dmd[nm]))

cat("\nVariant types:\n")
vt_dmd <- sort(table(rescued_dmd$Rescued_Variant_Type), decreasing = TRUE)
for (nm in names(vt_dmd)) cat(sprintf("  %-15s %d\n", nm, vt_dmd[nm]))

cat("\nGenomic regions:\n")
gr_dmd <- sort(table(rescued_dmd$Genomic_Region), decreasing = TRUE)
for (nm in names(gr_dmd)) cat(sprintf("  %-50s %d\n", nm, gr_dmd[nm]))

cat("\nMechanisms in rescued cases:\n")
for (mname in names(mech_names)) {
  mcol <- mech_names[[mname]]
  res <- count_mech(rescued_dmd, mcol)
  if (!is.null(res)) {
    for (j in seq_len(nrow(res))) {
      cat(sprintf("  [%s] %s: %d\n", mname, res$Value[j], res$N[j]))
    }
  }
}

# Missed by LRS
dmd_missed <- dmd %>% filter(grepl("Missed", LRS_Outcome))
if (nrow(dmd_missed) > 0) {
  cat(sprintf("\nMissed by LRS (N = %d):\n", nrow(dmd_missed)))
}

cat("\n")

###############################################################################
# 8d. EXPERT-SELECTED (Sorrentino_2025)
###############################################################################
cat(strrep("=", 80), "\n")
cat("8d. EXPERT-SELECTED COHORT (Sorrentino_2025)\n")
cat(strrep("=", 80), "\n\n")

sorr <- dat %>% filter(Study_ID == "Sorrentino_2025")
n_sorr <- nrow(sorr)

cat("Total N:", n_sorr, "\n")
cat("Phenotype:", unique(sorr$Phenotype_Cohort_Area), "\n")

cat("\nSR Result distribution:\n")
sorr_sr <- table(sorr$SR_Result, useNA = "ifany")
for (nm in names(sorr_sr)) {
  cat(sprintf("  %-25s %d (%s)\n", nm, sorr_sr[nm],
              sprintf("%.1f%%", sorr_sr[nm]/n_sorr*100)))
}

cat("\nLRS Outcome distribution:\n")
sorr_lrs <- table(sorr$LRS_Outcome, useNA = "ifany")
for (nm in names(sorr_lrs)) {
  cat(sprintf("  %-25s %d (%s)\n", nm, sorr_lrs[nm],
              sprintf("%.1f%%", sorr_lrs[nm]/n_sorr*100)))
}

rescued_sorr <- sorr %>% filter(is_rescued(LRS_Outcome))
n_resc_sorr  <- nrow(rescued_sorr)
cat(sprintf("\nRescue rate: %s\n", wilson_ci(n_resc_sorr, n_sorr)))

cat("Mechanisms summary (rescued cases):\n")
for (mname in names(mech_names)) {
  mcol <- mech_names[[mname]]
  res <- count_mech(rescued_sorr, mcol)
  if (!is.null(res)) {
    for (j in seq_len(nrow(res))) {
      cat(sprintf("  [%s] %s: %d\n", mname, res$Value[j], res$N[j]))
    }
  }
}

cat("\nVariant types (rescued):\n")
vt_sorr <- sort(table(rescued_sorr$Rescued_Variant_Type), decreasing = TRUE)
for (nm in names(vt_sorr)) cat(sprintf("  %-15s %d\n", nm, vt_sorr[nm]))

cat("\n")

###############################################################################
# 8e. LIMITED SCOPE (Shah_2025)
###############################################################################
cat(strrep("=", 80), "\n")
cat("8e. LIMITED SCOPE — Shah_2025\n")
cat(strrep("=", 80), "\n\n")

shah <- dat %>% filter(Study_ID == "Shah_2025")
n_shah <- nrow(shah)

cat("Total N:", n_shah, "\n")
cat("Phenotype:", unique(shah$Phenotype_Cohort_Area), "\n")
cat("Platform:", unique(shah$Platform), "\n")

# Coverage stats
cov_vals <- as.numeric(gsub("'", "", shah$Patient_Coverage))
cov_vals <- cov_vals[!is.na(cov_vals)]
if (length(cov_vals) > 0) {
  cat(sprintf("Coverage: median %.1f× (range %.1f–%.1f×)\n",
              median(cov_vals), min(cov_vals), max(cov_vals)))
}

cat("\nSR Result:", paste(names(table(shah$SR_Result)), table(shah$SR_Result),
                          sep = ": ", collapse = "; "), "\n")

cat("\nLRS Outcome:\n")
shah_lrs <- table(shah$LRS_Outcome, useNA = "ifany")
for (nm in names(shah_lrs)) {
  cat(sprintf("  %-25s %d (%s)\n", nm, shah_lrs[nm],
              sprintf("%.1f%%", shah_lrs[nm]/n_shah*100)))
}

cat(sprintf("\nRescue rate: %s\n", wilson_ci(sum(is_rescued(shah$LRS_Outcome)), n_shah)))

cat("\nEXCLUSION RATIONALE:\n")
cat("  • SV-only analysis at low coverage (~7×) — no SNV, repeat expansion,\n")
cat("    or methylation calling performed\n")
cat("  • Non-syndromic ASD cohort (narrow phenotype)\n")
cat("  • Only 1 VUS identified (SNAP25-AS1 inversion, inherited from\n")
cat("    unaffected mother) — does not meet definitive diagnostic threshold\n")
cat("  • Study design limits comparability with comprehensive LRS approaches\n\n")

###############################################################################
# 8f. AlAbdi_2023 — HARDCODED DESCRIPTIVE ANALYSIS
###############################################################################
cat(strrep("=", 80), "\n")
cat("8f. AlAbdi_2023 — DESCRIPTIVE ANALYSIS (HARDCODED VALUES)\n")
cat(strrep("=", 80), "\n\n")

cat("NOTE: AlAbdi_2023 data is not in the Dataset_S1.xlsx Group B+C sheet.\n")
cat("Values below are hardcoded per the study protocol.\n\n")

# Hardcoded values from protocol
N_al      <- 34
def_al    <- 9
poss_al   <- 4
rescued_al <- def_al + poss_al  # 13
missed_al  <- 2
not_al     <- 19

cat("SRS-negative cases analyzed: N =", N_al, "\n")
cat(sprintf("  Definitive Rescue: %d\n", def_al))
cat(sprintf("  Possible Rescue:   %d\n", poss_al))
cat(sprintf("  Missed by LRS:     %d\n", missed_al))
cat(sprintf("  Not Rescued:       %d\n", not_al))

# Verify sum
stopifnot(def_al + poss_al + missed_al + not_al == N_al)
cat(sprintf("  Sum check: %d + %d + %d + %d = %d ✓\n\n",
            def_al, poss_al, missed_al, not_al, N_al))

# Overall rescue rate
cat(sprintf("Overall rescue rate (Def + Poss): %s\n", wilson_ci(rescued_al, N_al)))

# Technical rescue = 6/34
tech_resc <- 6
cat(sprintf("LRS-technical-only rescue rate:   %s\n\n", wilson_ci(tech_resc, N_al)))

cat("Mechanisms (from protocol):\n")
cat("  Reanalysis:           7 cases\n")
cat("  SV detection:         4 cases\n")
cat("  Full-length assembly: 2 cases\n")
cat(sprintf("  Sum of mechanisms:    %d (some cases may have >1)\n", 7 + 4 + 2))

cat("\nLRS-technical rescue breakdown:\n")
cat("  Of the 13 rescued, 6 required LRS-specific capabilities\n")
cat("    (SV detection = 4, full-length assembly = 2)\n")
cat("  Remaining 7 were reanalysis findings (could potentially be\n")
cat("    found on SRS re-examination)\n\n")

cat("Key contextual notes:\n")
cat("  • Consanguineous cohort (Saudi Arabia)\n")
cat("  • AR enrichment expected given consanguinity\n")
cat("  • Study excluded from Group A due to potential confounding\n")
cat("    from consanguinity effects on variant spectrum\n\n")

# Wilson CIs for key proportions
cat("Wilson CIs for key proportions:\n")
cat(sprintf("  Overall rescue (13/34):   %s\n", wilson_ci(13, 34)))
cat(sprintf("  Technical rescue (6/34):  %s\n", wilson_ci(6, 34)))
cat(sprintf("  Definitive only (9/34):   %s\n", wilson_ci(9, 34)))
cat(sprintf("  Possible only (4/34):     %s\n", wilson_ci(4, 34)))
cat(sprintf("  Missed rate (2/34):       %s\n", wilson_ci(2, 34)))

cat("\n")

###############################################################################
# CROSS-STUDY SUMMARY TABLE
###############################################################################
cat(strrep("=", 80), "\n")
cat("CROSS-STUDY SUMMARY — ALL EXCLUDED STUDIES\n")
cat(strrep("=", 80), "\n\n")

# Build summary for each study in the file
all_studies <- sort(unique(dat$Study_ID))
cat(sprintf("%-20s %5s %5s %5s %8s   %s\n",
            "Study", "N", "Def", "Poss", "Rate", "Exclusion Reason"))
cat(strrep("-", 90), "\n")

exclusion_reasons <- c(
  Bruels_2022     = "Disease-specific (DMD)",
  Chu_2025        = "Disease-specific (DMD)",
  Daida_2025a     = "Second-hit only design",
  Daida_2025b     = "Second-hit only design",
  Ek_2025         = "Parallel/first-line LRS",
  Kamolvisit_2025 = "First-line LRS (no prior SRS)",
  Sano_2022       = "Second-hit only design",
  Shah_2025       = "SV-only, low cov, narrow pheno",
  Smits_2025      = "Parallel/first-line LRS",
  Sorrentino_2025 = "Expert-selected cases"
)

for (sid in all_studies) {
  sub  <- dat %>% filter(Study_ID == sid)
  n    <- nrow(sub)
  ndef <- sum(sub$LRS_Outcome == "Definitive_Rescue", na.rm = TRUE)
  npos <- sum(sub$LRS_Outcome == "Possible_Rescue", na.rm = TRUE)
  rate <- sprintf("%.1f%%", (ndef + npos)/n * 100)
  reason <- exclusion_reasons[sid]
  cat(sprintf("%-20s %5d %5d %5d %8s   %s\n", sid, n, ndef, npos, rate, reason))
}

# Add AlAbdi row
cat(sprintf("%-20s %5d %5d %5d %8s   %s\n",
            "AlAbdi_2023*", 34, 9, 4, "38.2%", "Consanguineous cohort"))
cat("\n* AlAbdi_2023 values hardcoded (not in dataset)\n")

###############################################################################
# PLATFORM SUMMARY
###############################################################################
cat(strrep("\n", 1))
cat(strrep("-", 60), "\n")
cat("Platform distribution across excluded studies:\n")
plat_study <- dat %>% 
  group_by(Study_ID, Platform) %>% 
  summarise(n = n(), .groups = "drop")
for (i in seq_len(nrow(plat_study))) {
  cat(sprintf("  %-20s %-10s %d\n",
              plat_study$Study_ID[i], plat_study$Platform[i], plat_study$n[i]))
}

###############################################################################
# DATA QUALITY FLAGS
###############################################################################
cat("\n")
cat(strrep("-", 60), "\n")
cat("Data quality checks:\n")
cat(sprintf("  Total rows: %d\n", nrow(dat)))
cat(sprintf("  Missing Study_ID: %d\n", sum(is.na(dat$Study_ID))))
cat(sprintf("  Missing LRS_Outcome: %d\n", sum(is.na(dat$LRS_Outcome))))
cat(sprintf("  Missing SR_Result: %d\n", sum(is.na(dat$SR_Result))))
cat(sprintf("  Unique patients (by Patient_ID_Systematic_review): %d\n",
            n_distinct(dat$Patient_ID_Systematic_review)))

# Check for LRS_Outcome variants with asterisks
asterisk_outcomes <- dat %>% filter(grepl("\\*", LRS_Outcome))
if (nrow(asterisk_outcomes) > 0) {
  cat(sprintf("\n  Asterisk-flagged outcomes (N = %d):\n", nrow(asterisk_outcomes)))
  ast_tbl <- table(asterisk_outcomes$LRS_Outcome)
  for (nm in names(ast_tbl)) {
    cat(sprintf("    %-25s %d\n", nm, ast_tbl[nm]))
  }
  cat("  NOTE: These outcomes have caveats noted by study authors.\n")
  cat("    Missed** (Smits): SR-positive cases where LRS failed to detect\n")
  cat("    Possible_Rescue** (Smits): LRS finding of uncertain significance\n")
  cat("    Missed*** (Chu): LRS initially failed, later resolved\n")
}

cat("\n")
cat(strrep("=", 80), "\n")
cat("END OF STEP 8 OUTPUT\n")
cat(strrep("=", 80), "\n")

sink()

cat("✓ Output written to: step8_excluded_studies_output.txt\n")
