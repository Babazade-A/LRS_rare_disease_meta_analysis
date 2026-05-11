###############################################################################
# STEP 4: RESCUE MECHANISMS — R VERIFICATION SCRIPT
# Systematic Review: LRS Rescue in SRS-Negative Patients
# Input: Dataset_S1.xlsx, sheet "Tier 1 - Group A" (11 studies, Group A SRS-negative patients)
#
# Analysis restricted to 55 rescued cases (Definitive + Possible Rescue).
# Descriptive analysis — no formal hypothesis testing on cross-tabulations
# (cells too sparse for 9 mechanisms × multiple phenotypes/prior testing).
#
# Sections:
#   1. Mechanism frequencies with Wilson CIs & subcategory breakdowns
#   2. Mechanism burden per case (0/1/2/3+, mean, median)
#   3. Cross-tabulation: Mechanism × Prior Testing (with stratified Wilson CIs)
#   4. Cross-tabulation: Mechanism × Phenotype (with stratified Wilson CIs)
#   5. Leave-one-out sensitivity on mechanism proportions
#   6. AlAbdi_2023 descriptive (technical vs reanalysis yield)
#   7. Internal consistency diagnostics
###############################################################################

# ── 0. LIBRARIES ─────────────────────────────────────────────────────────────
for (pkg in c("readxl", "dplyr", "tidyr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(readxl)
library(dplyr)
library(tidyr)

# ── 1a. SAVE OUTPUT TO TXT FILE ──────────────────────────────────────────────
output_file <- "Step4_Rescue_Mechanisms_OUTPUT.txt"
sink(output_file, split = TRUE)  # split=TRUE prints to console AND file

# ── 1. DATA LOADING ──────────────────────────────────────────────────────────
dat <- read_excel(
  "Dataset_S1.xlsx",
  sheet     = "Tier 1 - Group A",
  skip      = 1,
  col_names = TRUE
)

cat("═══════════════════════════════════════════════════════════════════════\n")
cat("  STEP 4 — RESCUE MECHANISMS: R VERIFICATION\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

N_total <- nrow(dat)

# ── 2. RESCUE CLASSIFICATION ─────────────────────────────────────────────────
dat <- dat %>%
  mutate(
    Rescued = as.integer(
      !is.na(LRS_Outcome) &
      (LRS_Outcome == "Definitive_Rescue" | grepl("Possible_Rescue", LRS_Outcome))
    ),
    Def_Rescued  = as.integer(!is.na(LRS_Outcome) & LRS_Outcome == "Definitive_Rescue"),
    Poss_Rescued = as.integer(!is.na(LRS_Outcome) & grepl("Possible_Rescue", LRS_Outcome))
  )

# Filter to rescued cases only
resc <- dat %>% filter(Rescued == 1)
N_resc <- nrow(resc)

cat(sprintf("Total Group A: N = %d\n", N_total))
cat(sprintf("Rescued cases: N = %d (Definitive: %d, Possible: %d)\n\n",
            N_resc, sum(resc$Def_Rescued), sum(resc$Poss_Rescued)))

# ── 3. PRIOR TESTING CLASSIFICATION (for rescued cases) ─────────────────────
dat <- dat %>%
  mutate(
    wes_yes   = tolower(trimws(as.character(Prior_WES))) == "yes",
    wgs_yes   = tolower(trimws(as.character(Prior_WGS))) == "yes",
    panel_yes = tolower(trimws(as.character(Prior_GenePanel))) == "yes",
    is_ES_or_GS = (!is.na(`Previous_Tests(or simultaneously)`) &
                    `Previous_Tests(or simultaneously)` == "ES_or_GS"),
    Prior_Category = case_when(
      is_ES_or_GS              ~ "ES_or_GS",
      wes_yes & wgs_yes        ~ "WES+WGS",
      wes_yes & !wgs_yes       ~ "WES-only",
      wgs_yes & !wes_yes       ~ "WGS-only",
      panel_yes & !wes_yes & !wgs_yes ~ "Panel-only",
      TRUE                     ~ "Other/Unknown"
    )
  )

# Re-filter rescued with new columns
resc <- dat %>% filter(Rescued == 1)

# ── 4. WILSON CI FUNCTION ────────────────────────────────────────────────────
wilson_ci <- function(x, n, conf.level = 0.95) {
  if (n == 0) return(c(lower = NA_real_, upper = NA_real_))
  z      <- qnorm(1 - (1 - conf.level) / 2)
  p      <- x / n
  denom  <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  margin <- (z / denom) * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  lower  <- max(0, center - margin)
  upper  <- min(1, center + margin)
  return(c(lower = lower, upper = upper))
}

# ── 5. MECHANISM COLUMN DEFINITIONS ──────────────────────────────────────────
# 8 mechanism columns with their full header names
mech_cols <- c(
  "Reanalysis_finding",
  "SV_Detection (SV_Deletion, SV_Duplication, Inversion, InDel_Del, Indel_Ins, Complex_SV, ME_Insertion, ME_Mediated_Del, Breakpoint_Resolution_only)",
  "Dark_Region_Resolution (Pseudogene_distinguishing, Segmental_Duplication, Low_Mappability, GC_Rich, Homopolymeric, Homologous_Region)",
  "Repeat_Expansion (Repeat_Expansion_Detection, Repeat_Sizing)",
  "Phasing(cis; trans)",
  "Full_Length_Sequencing (Non_Coding",
  "Methylation_Detection (Methylation_Profile, Regional_Methylation, Imprinting)",
  "Mosaicism"
)

# Short display names for readability
mech_short <- c(
  "Reanalysis",
  "SV_Detection",
  "Dark_Region",
  "Repeat_Expansion",
  "Phasing",
  "Full_Length_Seq",
  "Methylation",
  "Mosaicism"
)

# Verify columns exist
missing_mech <- setdiff(mech_cols, colnames(dat))
if (length(missing_mech) > 0) {
  cat("WARNING: Missing mechanism columns:\n")
  for (mc in missing_mech) cat(sprintf("  %s\n", mc))
} else {
  cat(sprintf("All %d mechanism columns found.\n\n", length(mech_cols)))
}

# ── 6. HELPER: Is mechanism present? ─────────────────────────────────────────
# Mechanism is PRESENT if value is non-empty, non-NA, not "None", not "No"
is_mech_present <- function(x) {
  !is.na(x) & trimws(as.character(x)) != "" &
  trimws(as.character(x)) != "None" &
  trimws(as.character(x)) != "No"
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: MECHANISM FREQUENCIES WITH WILSON CIs & SUBCATEGORIES
# ══════════════════════════════════════════════════════════════════════════════
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  SECTION 1: MECHANISM FREQUENCIES & SUBCATEGORIES\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

mech_summary <- data.frame(
  Mechanism = character(),
  N_pos     = integer(),
  Rate_pct  = numeric(),
  Wilson_lo = numeric(),
  Wilson_hi = numeric(),
  stringsAsFactors = FALSE
)

for (idx in seq_along(mech_cols)) {
  col <- mech_cols[idx]
  short <- mech_short[idx]

  present <- is_mech_present(resc[[col]])
  n_pos   <- sum(present)
  rate    <- 100 * n_pos / N_resc
  ci      <- wilson_ci(n_pos, N_resc)

  mech_summary <- rbind(mech_summary, data.frame(
    Mechanism = short,
    N_pos     = n_pos,
    Rate_pct  = rate,
    Wilson_lo = 100 * ci["lower"],
    Wilson_hi = 100 * ci["upper"],
    stringsAsFactors = FALSE
  ))

  # Subcategory breakdown
  vals <- resc[[col]][present]
  val_table <- sort(table(vals), decreasing = TRUE)

  cat(sprintf("%-20s: %2d / %d = %5.1f%%  [%5.1f%%, %5.1f%%]\n",
              short, n_pos, N_resc, rate, 100 * ci["lower"], 100 * ci["upper"]))
  if (length(val_table) > 0) {
    for (j in seq_along(val_table)) {
      subcat <- names(val_table)[j]
      subcount <- val_table[j]
      subpct <- 100 * subcount / n_pos
      cat(sprintf("    %-35s: %2d (%5.1f%% of mechanism, %5.1f%% of rescued)\n",
                  subcat, subcount, subpct, 100 * subcount / N_resc))
    }
  }
  cat("\n")
}

# Summary table
cat("\n── Mechanism Summary Table ──\n\n")
cat(sprintf("%-20s %5s %8s  [%s]\n", "Mechanism", "N", "Rate%", "95% Wilson CI"))
cat("──────────────────────────────────────────────────────────────\n")
for (i in seq_len(nrow(mech_summary))) {
  m <- mech_summary[i, ]
  cat(sprintf("%-20s %5d %7.1f%%  [%5.1f%%, %5.1f%%]\n",
              m$Mechanism, m$N_pos, m$Rate_pct, m$Wilson_lo, m$Wilson_hi))
}

# ── Phasing-only and Full-Length-Seq-only: denominator = rescued (N=55) ───────
cat("\n── Phasing & Full-Length Sequencing: Rescued Cases (N=55) ──\n\n")

phasing_col    <- mech_cols[which(mech_short == "Phasing")]
fulllength_col <- mech_cols[which(mech_short == "Full_Length_Seq")]

# Counts among rescued cases (already in Section 1, but isolated here)
phasing_resc    <- sum(is_mech_present(resc[[phasing_col]]))
fulllength_resc <- sum(is_mech_present(resc[[fulllength_col]]))

ci_phasing_resc    <- wilson_ci(phasing_resc, N_resc)
ci_fulllength_resc <- wilson_ci(fulllength_resc, N_resc)

cat(sprintf("  Phasing (rescued):          %d / %d = %.1f%%  [%.1f%%, %.1f%%]\n",
            phasing_resc, N_resc, 100 * phasing_resc / N_resc,
            100 * ci_phasing_resc["lower"], 100 * ci_phasing_resc["upper"]))
cat(sprintf("  Full_Length_Seq (rescued):   %d / %d = %.1f%%  [%.1f%%, %.1f%%]\n",
            fulllength_resc, N_resc, 100 * fulllength_resc / N_resc,
            100 * ci_fulllength_resc["lower"], 100 * ci_fulllength_resc["upper"]))

# Subcategory breakdown among rescued
cat("\n  Phasing subcategories (rescued):\n")
phasing_vals <- resc[[phasing_col]][is_mech_present(resc[[phasing_col]])]
if (length(phasing_vals) > 0) {
  pt <- sort(table(phasing_vals), decreasing = TRUE)
  for (j in seq_along(pt)) cat(sprintf("    %s: %d\n", names(pt)[j], pt[j]))
} else {
  cat("    (none)\n")
}

cat("\n  Full_Length_Seq subcategories (rescued):\n")
fl_vals <- resc[[fulllength_col]][is_mech_present(resc[[fulllength_col]])]
if (length(fl_vals) > 0) {
  ft <- sort(table(fl_vals), decreasing = TRUE)
  for (j in seq_along(ft)) cat(sprintf("    %s: %d\n", names(ft)[j], ft[j]))
} else {
  cat("    (none)\n")
}
cat("\n")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2: MECHANISM BURDEN PER CASE
# ══════════════════════════════════════════════════════════════════════════════
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  SECTION 2: MECHANISM BURDEN PER CASE\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# Count mechanisms per rescued case
resc$mech_count <- 0L
for (col in mech_cols) {
  resc$mech_count <- resc$mech_count + as.integer(is_mech_present(resc[[col]]))
}

burden_table <- table(resc$mech_count)
cat("Mechanism burden distribution:\n")
cat("──────────────────────────────────────────────────────────────\n")
for (b in sort(as.integer(names(burden_table)))) {
  n_b <- burden_table[as.character(b)]
  ci_b <- wilson_ci(n_b, N_resc)
  cat(sprintf("  %d mechanism(s): %2d / %d = %5.1f%%  [%5.1f%%, %5.1f%%]\n",
              b, n_b, N_resc, 100 * n_b / N_resc,
              100 * ci_b["lower"], 100 * ci_b["upper"]))
}

mean_burden   <- mean(resc$mech_count)
median_burden <- median(resc$mech_count)
total_instances <- sum(resc$mech_count)

cat(sprintf("\n  Mean mechanisms per case:   %.1f\n", mean_burden))
cat(sprintf("  Median mechanisms per case: %.0f\n", median_burden))
cat(sprintf("  Total mechanism instances:  %d (across %d rescued cases)\n",
            total_instances, N_resc))

# Summation check: total instances should equal sum of mechanism N_pos values
sum_mech_pos <- sum(mech_summary$N_pos)
cat(sprintf("\n  Summation check: sum(N_pos) = %d, sum(mech_count) = %d → %s\n",
            sum_mech_pos, total_instances,
            ifelse(sum_mech_pos == total_instances, "PASS", "FAIL")))

# Cases with 0 mechanisms (rescued but no mechanism identified)
zero_mech <- resc %>% filter(mech_count == 0)
if (nrow(zero_mech) > 0) {
  cat(sprintf("\n  Cases with 0 mechanisms (%d):\n", nrow(zero_mech)))
  print(as.data.frame(zero_mech %>%
    select(Patient_ID_Systematic_review, Study_ID, LRS_Outcome, Phenotype_Cohort_Area)))
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3: MECHANISM × PRIOR TESTING CROSS-TABULATION
# ══════════════════════════════════════════════════════════════════════════════
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  SECTION 3: MECHANISM × PRIOR TESTING CROSS-TABULATION\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# Prior testing distribution among rescued
prior_resc <- table(resc$Prior_Category)
cat("Prior testing of rescued cases:\n")
print(prior_resc)
cat("\n")

# Cross-tabulation: count of each mechanism within each prior testing category
prior_cats_order <- c("WES-only", "WGS-only", "WES+WGS", "Panel-only", "Other/Unknown")
# Only include categories present in rescued
prior_cats_present <- prior_cats_order[prior_cats_order %in% names(prior_resc)]

cat("Mechanism × Prior Testing (counts):\n")
cat("──────────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-20s", "Mechanism"))
for (pc in prior_cats_present) cat(sprintf(" %10s", pc))
cat("      TOTAL\n")
cat("──────────────────────────────────────────────────────────────────────\n")

cross_prior <- matrix(0L, nrow = length(mech_short), ncol = length(prior_cats_present))
rownames(cross_prior) <- mech_short
colnames(cross_prior) <- prior_cats_present

for (idx in seq_along(mech_cols)) {
  col   <- mech_cols[idx]
  short <- mech_short[idx]

  for (pc_idx in seq_along(prior_cats_present)) {
    pc <- prior_cats_present[pc_idx]
    sub <- resc %>% filter(Prior_Category == pc)
    cross_prior[short, pc] <- sum(is_mech_present(sub[[col]]))
  }

  cat(sprintf("%-20s", short))
  for (pc in prior_cats_present) {
    cat(sprintf(" %10d", cross_prior[short, pc]))
  }
  cat(sprintf("    %5d\n", mech_summary$N_pos[idx]))
}

# Denominators row
cat(sprintf("%-20s", "N (rescued)"))
for (pc in prior_cats_present) cat(sprintf(" %10d", prior_resc[pc]))
cat(sprintf("    %5d\n", N_resc))

# Stratified rates with Wilson CIs
cat("\n── Stratified Rates (%) with Wilson CIs ──\n\n")
for (idx in seq_along(mech_cols)) {
  short <- mech_short[idx]
  cat(sprintf("  %s:\n", short))
  for (pc in prior_cats_present) {
    n_pc <- as.integer(prior_resc[pc])
    r_pc <- cross_prior[short, pc]
    if (n_pc > 0) {
      rate <- 100 * r_pc / n_pc
      ci <- wilson_ci(r_pc, n_pc)
      cat(sprintf("    %-15s: %2d/%2d = %5.1f%%  [%5.1f%%, %5.1f%%]\n",
                  pc, r_pc, n_pc, rate, 100 * ci["lower"], 100 * ci["upper"]))
    }
  }
  cat("\n")
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4: MECHANISM × PHENOTYPE CROSS-TABULATION
# ══════════════════════════════════════════════════════════════════════════════
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  SECTION 4: MECHANISM × PHENOTYPE CROSS-TABULATION\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# Phenotype distribution of rescued cases
pheno_resc <- sort(table(resc$Phenotype_Cohort_Area), decreasing = TRUE)
cat("Phenotype distribution of rescued cases:\n")
for (j in seq_along(pheno_resc)) {
  p <- names(pheno_resc)[j]
  n <- pheno_resc[j]
  cat(sprintf("  %-35s: %2d (%5.1f%%)\n", p, n, 100 * n / N_resc))
}
cat(sprintf("  %-35s: %2d\n", "TOTAL", N_resc))

# Phenotypes with N >= 3 for cross-tab (more lenient than testing threshold)
pheno_for_crosstab <- names(pheno_resc[pheno_resc >= 2])

cat(sprintf("\nPhenotypes in cross-tabulation (N >= 2): %d\n\n", length(pheno_for_crosstab)))

# Cross-tabulation counts
cat("Mechanism × Phenotype (counts):\n")
cat("──────────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-20s", "Mechanism"))
for (ph in pheno_for_crosstab) cat(sprintf(" %8s", substr(ph, 1, 8)))
cat("\n")
cat("──────────────────────────────────────────────────────────────────────\n")

cross_pheno <- matrix(0L, nrow = length(mech_short), ncol = length(pheno_for_crosstab))
rownames(cross_pheno) <- mech_short
colnames(cross_pheno) <- pheno_for_crosstab

for (idx in seq_along(mech_cols)) {
  col   <- mech_cols[idx]
  short <- mech_short[idx]

  for (ph_idx in seq_along(pheno_for_crosstab)) {
    ph <- pheno_for_crosstab[ph_idx]
    sub <- resc %>% filter(Phenotype_Cohort_Area == ph)
    cross_pheno[short, ph] <- sum(is_mech_present(sub[[col]]))
  }

  cat(sprintf("%-20s", short))
  for (ph in pheno_for_crosstab) {
    cat(sprintf(" %8d", cross_pheno[short, ph]))
  }
  cat("\n")
}

# Denominator row
cat(sprintf("%-20s", "N (rescued)"))
for (ph in pheno_for_crosstab) cat(sprintf(" %8d", pheno_resc[ph]))
cat("\n")

# Stratified rates with Wilson CIs for phenotypes with N >= 5
cat("\n── Stratified Rates by Phenotype (N >= 5 only) ──\n\n")
pheno_for_ci <- names(pheno_resc[pheno_resc >= 5])
for (idx in seq_along(mech_cols)) {
  short <- mech_short[idx]
  cat(sprintf("  %s:\n", short))
  for (ph in pheno_for_ci) {
    n_ph <- as.integer(pheno_resc[ph])
    r_ph <- cross_pheno[short, ph]
    rate <- 100 * r_ph / n_ph
    ci <- wilson_ci(r_ph, n_ph)
    cat(sprintf("    %-30s: %2d/%2d = %5.1f%%  [%5.1f%%, %5.1f%%]\n",
                ph, r_ph, n_ph, rate, 100 * ci["lower"], 100 * ci["upper"]))
  }
  cat("\n")
}

cat("NOTE: No formal hypothesis testing on mechanism × phenotype\n")
cat("cross-tabulation (cells too sparse for 9 × 14 table).\n\n")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5: LEAVE-ONE-OUT SENSITIVITY ON MECHANISM PROPORTIONS
# ══════════════════════════════════════════════════════════════════════════════
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  SECTION 5: LEAVE-ONE-OUT SENSITIVITY\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# Get studies with rescued cases
studies_with_rescue <- resc %>%
  group_by(Study_ID) %>%
  summarise(N_rescued = n(), .groups = "drop") %>%
  arrange(desc(N_rescued))

cat("Studies contributing rescued cases:\n")
for (i in seq_len(nrow(studies_with_rescue))) {
  s <- studies_with_rescue[i, ]
  cat(sprintf("  %-18s: %d rescued\n", s$Study_ID, s$N_rescued))
}
cat("\n")

# For each mechanism, do LOO by study
for (idx in seq_along(mech_cols)) {
  col   <- mech_cols[idx]
  short <- mech_short[idx]
  baseline_n   <- mech_summary$N_pos[idx]
  baseline_rate <- mech_summary$Rate_pct[idx]

  cat(sprintf("── %s (baseline: %d/%d = %.1f%%) ──\n", short, baseline_n, N_resc, baseline_rate))

  for (j in seq_len(nrow(studies_with_rescue))) {
    sid   <- studies_with_rescue$Study_ID[j]
    n_sid <- studies_with_rescue$N_rescued[j]

    # Mechanism count in this study's rescued cases
    study_resc <- resc %>% filter(Study_ID == sid)
    mech_in_study <- sum(is_mech_present(study_resc[[col]]))

    loo_resc <- N_resc - n_sid
    loo_pos  <- baseline_n - mech_in_study
    loo_rate <- 100 * loo_pos / loo_resc
    change   <- loo_rate - baseline_rate
    ci_loo   <- wilson_ci(loo_pos, loo_resc)

    cat(sprintf("  Drop %-18s: %2d/%2d = %5.1f%% (%+5.1fpp)  [%5.1f%%, %5.1f%%]\n",
                sid, loo_pos, loo_resc, loo_rate, change,
                100 * ci_loo["lower"], 100 * ci_loo["upper"]))
  }
  cat("\n")
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6: AlAbdi_2023 DESCRIPTIVE (hardcoded from protocol)
# ══════════════════════════════════════════════════════════════════════════════
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  SECTION 6: AlAbdi_2023 DESCRIPTIVE (reported values)\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# AlAbdi_2023: N=34 SRS-negative
# 9 Definitive + 4 Possible = 13 rescued, 2 Missed, 19 Not Rescued
alabdi_n      <- 34
alabdi_resc   <- 13   # 9 Def + 4 Poss
alabdi_missed <- 2
alabdi_notr   <- 19

cat(sprintf("AlAbdi_2023: N = %d SRS-negative\n", alabdi_n))
cat(sprintf("  Definitive: 9, Possible: 4, Rescued total: %d\n", alabdi_resc))
cat(sprintf("  Missed: %d, Not Rescued: %d\n", alabdi_missed, alabdi_notr))
cat(sprintf("  Summation: %d + %d + %d = %d (expected %d) → %s\n\n",
            alabdi_resc, alabdi_missed, alabdi_notr,
            alabdi_resc + alabdi_missed + alabdi_notr, alabdi_n,
            ifelse(alabdi_resc + alabdi_missed + alabdi_notr == alabdi_n, "PASS", "FAIL")))

# Mechanism breakdown
alabdi_reanalysis <- 7
alabdi_sv         <- 4
alabdi_fulllength <- 2

# Technical rescue = SV + Full_Length
alabdi_technical <- alabdi_sv + alabdi_fulllength  # = 6

ci_resc_al  <- wilson_ci(alabdi_resc, alabdi_n)
ci_tech_al  <- wilson_ci(alabdi_technical, alabdi_n)
ci_reana_al <- wilson_ci(alabdi_reanalysis, alabdi_n)

cat("Rescue yield:\n")
cat(sprintf("  Overall rescue:      %d/%d = %.1f%%  [%.1f%%, %.1f%%]\n",
            alabdi_resc, alabdi_n, 100 * alabdi_resc / alabdi_n,
            100 * ci_resc_al["lower"], 100 * ci_resc_al["upper"]))

cat(sprintf("\nMechanism breakdown:\n"))
cat(sprintf("  Reanalysis/reinterpretation: %d/%d = %.1f%%  [%.1f%%, %.1f%%]\n",
            alabdi_reanalysis, alabdi_n, 100 * alabdi_reanalysis / alabdi_n,
            100 * ci_reana_al["lower"], 100 * ci_reana_al["upper"]))
cat(sprintf("  Technical rescue (SV + Full-length): %d/%d = %.1f%%  [%.1f%%, %.1f%%]\n",
            alabdi_technical, alabdi_n, 100 * alabdi_technical / alabdi_n,
            100 * ci_tech_al["lower"], 100 * ci_tech_al["upper"]))
cat(sprintf("    - SV Detection:         %d\n", alabdi_sv))
cat(sprintf("    - Full-length assembly: %d\n", alabdi_fulllength))

cat(sprintf("\n  Technical vs Reanalysis ratio: %.1f : %.1f\n",
            alabdi_technical, alabdi_reanalysis))
cat(sprintf("  Note: Some cases may have overlapping mechanisms\n"))
cat(sprintf("  (reanalysis + technical), so counts may not sum to 13.\n\n"))


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 7: INTERNAL CONSISTENCY DIAGNOSTICS
# ══════════════════════════════════════════════════════════════════════════════
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  SECTION 7: INTERNAL CONSISTENCY DIAGNOSTICS\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# Check 1: Rescued count
check_rescued <- N_resc == 55

# Check 2: Mechanism instances sum
check_mech_sum <- sum_mech_pos == total_instances

# Check 3: All mechanism columns found
check_cols <- length(missing_mech) == 0

# Check 4: Burden distribution sums to N_resc
check_burden <- sum(burden_table) == N_resc

# Check 5: Prior testing of rescued sums to N_resc
check_prior_resc <- sum(prior_resc) == N_resc

# Check 6: Phenotype of rescued sums to N_resc
check_pheno_resc <- sum(pheno_resc) == N_resc

# Check 7: AlAbdi summation
check_alabdi <- (alabdi_resc + alabdi_missed + alabdi_notr) == alabdi_n

# Check 8: NaN in mechanism summary
check_nan <- !any(is.nan(mech_summary$Rate_pct)) &
             !any(is.nan(mech_summary$Wilson_lo)) &
             !any(is.nan(mech_summary$Wilson_hi))

checks <- c(
  rescued_N       = check_rescued,
  mech_instances  = check_mech_sum,
  cols_found      = check_cols,
  burden_sum      = check_burden,
  prior_resc_sum  = check_prior_resc,
  pheno_resc_sum  = check_pheno_resc,
  alabdi_sum      = check_alabdi,
  no_nan          = check_nan
)

cat("Diagnostic checks:\n")
for (nm in names(checks)) {
  cat(sprintf("  %-20s: %s\n", nm, ifelse(checks[nm], "PASS", "FAIL")))
}
cat(sprintf("\nOverall: %s\n",
            ifelse(all(checks), "ALL CHECKS PASSED", "ISSUES DETECTED")))


# ══════════════════════════════════════════════════════════════════════════════
# EXACT VALUES FOR MANUSCRIPT CROSS-CHECK
# ══════════════════════════════════════════════════════════════════════════════
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  EXACT VALUES (UNROUNDED) FOR MANUSCRIPT CROSS-CHECK\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

cat("── Mechanism Frequencies ──\n")
for (i in seq_len(nrow(mech_summary))) {
  m <- mech_summary[i, ]
  cat(sprintf("  %-20s: %d/%d = %.10f  CI=[%.10f, %.10f]\n",
              m$Mechanism, m$N_pos, N_resc,
              m$Rate_pct / 100, m$Wilson_lo / 100, m$Wilson_hi / 100))
}

cat("\n── Burden Distribution ──\n")
for (b in sort(as.integer(names(burden_table)))) {
  n_b <- burden_table[as.character(b)]
  cat(sprintf("  %d mechanisms: %d/%d = %.10f\n", b, n_b, N_resc, n_b / N_resc))
}
cat(sprintf("  Mean: %.10f\n", mean_burden))
cat(sprintf("  Median: %.1f\n", median_burden))

cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat("  STEP 4 VERIFICATION COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════════════\n")

# ── Close output file ────────────────────────────────────────────────────────
sink()
cat(sprintf("\nOutput saved to: %s\n", output_file))
