###############################################################################
# STEP 1: OVERALL RESCUE RATES — R VERIFICATION SCRIPT
# Systematic Review: LRS Rescue in SRS-Negative Patients
# Input: Dataset_S1.xlsx, sheet "Tier 1 - Group A"
#        (11 studies, N=538 Group A SRS-negative patients after dedup)
# Methods: Wilson Score CIs, Leave-One-Out sensitivity, descriptive counts
###############################################################################

# ── 0. LIBRARIES ─────────────────────────────────────────────────────────────
if (!requireNamespace("readxl", quietly = TRUE))   install.packages("readxl")
if (!requireNamespace("dplyr", quietly = TRUE))     install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE))     install.packages("tidyr")

library(readxl)
library(dplyr)
library(tidyr)

# ── 1a. SAVE OUTPUT TO TXT FILE ──────────────────────────────────────────────
output_file <- "Step1_Overall_Rescue_Rates_OUTPUT.txt"
sink(output_file, split = TRUE)  # split=TRUE prints to console AND file

# ── 1. DATA LOADING ──────────────────────────────────────────────────────────
# Dataset_S1.xlsx, sheet "Tier 1 - Group A":
#   Row 1 = section banner ("Study info"), Row 2 = column names, data from row 3+
# Header structure: row 1 banner + row 2 column names, so skip=1 still works.
dat <- read_excel(
  "Dataset_S1.xlsx",
  sheet     = "Tier 1 - Group A",
  skip      = 1,          # skip the group-header row (row 1)
  col_names = TRUE        # row 2 becomes header
)

cat("═══════════════════════════════════════════════════════════════════\n")
cat("  STEP 1 — OVERALL RESCUE RATES: R VERIFICATION\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

# ── 2. INTERNAL CONSISTENCY DIAGNOSTICS (Section 5) ─────────────────────────
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  SECTION 5: INTERNAL CONSISTENCY DIAGNOSTICS\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

N_total <- nrow(dat)
cat(sprintf("Total rows loaded: %d\n", N_total))

# NaN / NA verification in LRS_Outcome
na_outcome <- sum(is.na(dat$LRS_Outcome) | trimws(as.character(dat$LRS_Outcome)) == "")
cat(sprintf("Missing/empty LRS_Outcome: %d\n", na_outcome))

# Verify all Previous_SR_Result == "Negative"
sr_vals <- table(dat$`Previous_SR_Result (or simultaneous)`, useNA = "ifany")
cat("\nPrevious_SR_Result distribution:\n")
print(sr_vals)

# Outcome summation check
outcome_table <- table(dat$LRS_Outcome, useNA = "ifany")
cat("\nLRS_Outcome distribution (raw counts):\n")
print(outcome_table)
cat(sprintf("\nSum of all outcome categories: %d\n", sum(outcome_table)))
cat(sprintf("Matches total N (%d): %s\n\n", N_total, sum(outcome_table) == N_total))

# ── 3. CREATE RESCUE CLASSIFICATIONS (Section 3) ────────────────────────────
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  SECTION 3: QUANTITATIVE RESCUE CLASSIFICATIONS\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

dat <- dat %>%
  mutate(
    # Rescued: Definitive OR any Possible (including Possible_Rescue**)
    Rescued = as.integer(
      !is.na(LRS_Outcome) &
      (LRS_Outcome == "Definitive_Rescue" | grepl("Possible_Rescue", LRS_Outcome))
    ),
    # Definitive only
    Def_Rescued = as.integer(!is.na(LRS_Outcome) & LRS_Outcome == "Definitive_Rescue"),
    # Possible only
    Poss_Rescued = as.integer(!is.na(LRS_Outcome) & grepl("Possible_Rescue", LRS_Outcome)),
    # Not rescued
    Not_Rescued = as.integer(!is.na(LRS_Outcome) & LRS_Outcome == "Not_Rescued"),
    # Missed (SRS found something LRS did not)
    Missed = as.integer(!is.na(LRS_Outcome) & LRS_Outcome == "Missed"),
    # Concordant
    Concordant = as.integer(!is.na(LRS_Outcome) & LRS_Outcome == "Concordant"),
    # Refined
    Refined = as.integer(!is.na(LRS_Outcome) & LRS_Outcome == "Refined")
  )

n_def   <- sum(dat$Def_Rescued)
n_poss  <- sum(dat$Poss_Rescued)
n_comb  <- sum(dat$Rescued)
n_notr  <- sum(dat$Not_Rescued)
n_miss  <- sum(dat$Missed)
n_conc  <- sum(dat$Concordant)
n_refi  <- sum(dat$Refined)

cat("Rescue Classification Counts (Group A, N =", N_total, "):\n")
cat(sprintf("  Definitive Rescue:  %d\n", n_def))
cat(sprintf("  Possible Rescue:    %d\n", n_poss))
cat(sprintf("  Combined Rescue:    %d  (Def + Poss = %d + %d)\n", n_comb, n_def, n_poss))
cat(sprintf("  Not Rescued:        %d\n", n_notr))
cat(sprintf("  Missed:             %d\n", n_miss))
cat(sprintf("  Concordant:         %d\n", n_conc))
cat(sprintf("  Refined:            %d\n", n_refi))
cat(sprintf("\n  Summation check: %d + %d + %d + %d + %d = %d (expected %d) → %s\n\n",
            n_comb, n_notr, n_miss, n_conc, n_refi,
            n_comb + n_notr + n_miss + n_conc + n_refi,
            N_total,
            ifelse(n_comb + n_notr + n_miss + n_conc + n_refi == N_total, "PASS", "FAIL")))

# Percentages
cat("Rescue Rates (%):\n")
cat(sprintf("  Definitive:   %d / %d = %.6f%%\n", n_def, N_total, 100 * n_def / N_total))
cat(sprintf("  Possible:     %d / %d = %.6f%%\n", n_poss, N_total, 100 * n_poss / N_total))
cat(sprintf("  Combined:     %d / %d = %.6f%%\n", n_comb, N_total, 100 * n_comb / N_total))
cat(sprintf("  Not Rescued:  %d / %d = %.6f%%\n", n_notr, N_total, 100 * n_notr / N_total))

# ── 4. WILSON SCORE CONFIDENCE INTERVALS (Section 1) ────────────────────────
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  SECTION 1: CONFIDENCE INTERVAL PARAMETERS (WILSON SCORE)\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

#' Wilson Score Interval (two-sided, 95%)
#' Uses the standard Wilson formula — equivalent to prop.test(correct=FALSE)
wilson_ci <- function(x, n, conf.level = 0.95) {
  if (n == 0) return(c(lower = NA_real_, upper = NA_real_))
  z   <- qnorm(1 - (1 - conf.level) / 2)
  p   <- x / n
  denom <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  margin <- (z / denom) * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  lower  <- max(0, center - margin)
  upper  <- min(1, center + margin)
  return(c(lower = lower, upper = upper))
}

# ── 4a. Study-Level CIs ─────────────────────────────────────────────────────
study_summary <- dat %>%
  group_by(Study_ID) %>%
  summarise(
    N          = n(),
    Def        = sum(Def_Rescued),
    Poss       = sum(Poss_Rescued),
    Combined   = sum(Rescued),
    Not_Resc   = sum(Not_Rescued),
    .groups    = "drop"
  ) %>%
  arrange(desc(N))

# Compute Wilson CIs for each study
study_summary <- study_summary %>%
  rowwise() %>%
  mutate(
    Rate_pct    = 100 * Combined / N,
    Wilson_lo   = 100 * wilson_ci(Combined, N)["lower"],
    Wilson_hi   = 100 * wilson_ci(Combined, N)["upper"],
    Def_rate    = 100 * Def / N,
    Def_lo      = 100 * wilson_ci(Def, N)["lower"],
    Def_hi      = 100 * wilson_ci(Def, N)["upper"]
  ) %>%
  ungroup()

cat("Study-Level Rescue Rates with 95% Wilson CIs:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-18s %5s %5s %5s %5s %8s  [%s]\n",
            "Study", "N", "Def", "Poss", "Comb", "Rate%", "95% Wilson CI"))
cat("─────────────────────────────────────────────────────────────────\n")
for (i in seq_len(nrow(study_summary))) {
  s <- study_summary[i, ]
  cat(sprintf("%-18s %5d %5d %5d %5d %8.1f%%  [%5.1f%%, %5.1f%%]\n",
              s$Study_ID, s$N, s$Def, s$Poss, s$Combined,
              s$Rate_pct, s$Wilson_lo, s$Wilson_hi))
}

# ── 4b. Overall Group A CIs ─────────────────────────────────────────────────
cat("\n── Overall Group A CIs ──\n\n")

# Combined rescue CI
ci_comb <- wilson_ci(n_comb, N_total)
cat(sprintf("Combined Rescue: %d / %d = %.6f%%\n", n_comb, N_total, 100 * n_comb / N_total))
cat(sprintf("  95%% Wilson CI: [%.6f%%, %.6f%%]\n", 100 * ci_comb["lower"], 100 * ci_comb["upper"]))
cat(sprintf("  Rounded:       [%.1f%%, %.1f%%]\n\n", 100 * ci_comb["lower"], 100 * ci_comb["upper"]))

# Cross-check with prop.test (Wilson without continuity correction)
pt_comb <- prop.test(n_comb, N_total, correct = FALSE)
cat(sprintf("  prop.test cross-check: [%.6f%%, %.6f%%]\n\n",
            100 * pt_comb$conf.int[1], 100 * pt_comb$conf.int[2]))

# Definitive rescue CI
ci_def <- wilson_ci(n_def, N_total)
cat(sprintf("Definitive Rescue: %d / %d = %.6f%%\n", n_def, N_total, 100 * n_def / N_total))
cat(sprintf("  95%% Wilson CI: [%.6f%%, %.6f%%]\n", 100 * ci_def["lower"], 100 * ci_def["upper"]))
cat(sprintf("  Rounded:       [%.1f%%, %.1f%%]\n\n", 100 * ci_def["lower"], 100 * ci_def["upper"]))

# Possible rescue CI
ci_poss <- wilson_ci(n_poss, N_total)
cat(sprintf("Possible Rescue: %d / %d = %.6f%%\n", n_poss, N_total, 100 * n_poss / N_total))
cat(sprintf("  95%% Wilson CI: [%.6f%%, %.6f%%]\n\n", 100 * ci_poss["lower"], 100 * ci_poss["upper"]))

# ── 4c. Boundary Validation ─────────────────────────────────────────────────
cat("── Boundary Validation (all CIs must be within [0, 1]) ──\n\n")
all_lowers <- c(study_summary$Wilson_lo, 100 * ci_comb["lower"],
                100 * ci_def["lower"], 100 * ci_poss["lower"])
all_uppers <- c(study_summary$Wilson_hi, 100 * ci_comb["upper"],
                100 * ci_def["upper"], 100 * ci_poss["upper"])
boundary_ok <- all(all_lowers >= 0) & all(all_uppers <= 100)
cat(sprintf("  Min lower bound: %.6f%% (≥0: %s)\n", min(all_lowers), min(all_lowers) >= 0))
cat(sprintf("  Max upper bound: %.6f%% (≤100: %s)\n", max(all_uppers), max(all_uppers) <= 100))
cat(sprintf("  Boundary validation: %s\n", ifelse(boundary_ok, "PASS", "FAIL")))

# ── 4d. Study-level CI extremes ──────────────────────────────────────────────
min_study <- study_summary %>% filter(Wilson_lo == min(Wilson_lo))
max_study <- study_summary %>% filter(Wilson_hi == max(Wilson_hi))
cat(sprintf("\n  Narrowest lower bound: %s [%.1f%%, %.1f%%]\n",
            min_study$Study_ID, min_study$Wilson_lo, min_study$Wilson_hi))
cat(sprintf("  Widest upper bound:   %s [%.1f%%, %.1f%%]\n",
            max_study$Study_ID, max_study$Wilson_lo, max_study$Wilson_hi))

# ── 5. LEAVE-ONE-OUT ANALYSIS (Section 2) ────────────────────────────────────
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  SECTION 2: LEAVE-ONE-OUT (LOO) ANALYSIS\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

baseline_rate <- 100 * n_comb / N_total
cat(sprintf("Baseline: %d rescued / %d total = %.6f%%\n\n", n_comb, N_total, baseline_rate))

loo_results <- study_summary %>%
  rowwise() %>%
  mutate(
    LOO_N        = N_total - N,
    LOO_Rescued  = n_comb - Combined,
    LOO_Rate_pct = 100 * LOO_Rescued / LOO_N,
    Change_pp    = LOO_Rate_pct - baseline_rate,
    LOO_Wilson_lo = 100 * wilson_ci(LOO_Rescued, LOO_N)["lower"],
    LOO_Wilson_hi = 100 * wilson_ci(LOO_Rescued, LOO_N)["upper"]
  ) %>%
  ungroup() %>%
  arrange(desc(Change_pp))

cat(sprintf("%-18s %5s %5s %8s %9s  [%s]\n",
            "Dropped Study", "N_LOO", "Resc", "Rate%", "Δ(pp)", "95% Wilson CI"))
cat("────────────────────────────────────────────────────────────────────────\n")
for (i in seq_len(nrow(loo_results))) {
  r <- loo_results[i, ]
  cat(sprintf("%-18s %5d %5d %8.1f%% %+8.1fpp  [%5.1f%%, %5.1f%%]\n",
              r$Study_ID, r$LOO_N, r$LOO_Rescued,
              r$LOO_Rate_pct, r$Change_pp,
              r$LOO_Wilson_lo, r$LOO_Wilson_hi))
}

# LOO N-summation verification
cat("\n── LOO N-Summation Verification ──\n")
for (i in seq_len(nrow(loo_results))) {
  r <- loo_results[i, ]
  orig_n <- study_summary %>% filter(Study_ID == r$Study_ID) %>% pull(N)
  check  <- r$LOO_N + orig_n == N_total
  cat(sprintf("  Drop %-18s: %d + %d = %d (expected %d) → %s\n",
              r$Study_ID, r$LOO_N, orig_n, r$LOO_N + orig_n, N_total,
              ifelse(check, "PASS", "FAIL")))
}

# Robustness range
loo_min <- min(loo_results$LOO_Rate_pct)
loo_max <- max(loo_results$LOO_Rate_pct)
cat(sprintf("\nRobustness range of LOO rates: %.1f%% – %.1f%%\n", loo_min, loo_max))
cat(sprintf("Baseline rate: %.1f%%\n", round(baseline_rate, 1)))
cat(sprintf("Max deviation: %.1f pp\n", max(abs(loo_results$Change_pp))))

# Most influential study
most_up   <- loo_results %>% filter(Change_pp == max(Change_pp))
most_down <- loo_results %>% filter(Change_pp == min(Change_pp))
cat(sprintf("\nMost influential (dropping raises rate):  %s (Δ = %+.1f pp)\n",
            most_up$Study_ID, most_up$Change_pp))
cat(sprintf("Most influential (dropping lowers rate):  %s (Δ = %+.1f pp)\n",
            most_down$Study_ID, most_down$Change_pp))

# ── 6. MATHEMATICAL COMPARISON & WEIGHTS (Section 4) ─────────────────────────
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  SECTION 4: MATHEMATICAL COMPARISON & WEIGHTS\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# Study weights (share of denominator)
cat("Study Weights (share of total denominator):\n")
study_weights <- study_summary %>%
  mutate(Weight_pct = 100 * N / N_total) %>%
  arrange(desc(Weight_pct))

for (i in seq_len(nrow(study_weights))) {
  w <- study_weights[i, ]
  cat(sprintf("  %-18s: %5d / %d = %6.1f%%\n", w$Study_ID, w$N, N_total, w$Weight_pct))
}

# Heterogeneity range
rate_min <- min(study_summary$Rate_pct)
rate_max <- max(study_summary$Rate_pct)
cat(sprintf("\nHeterogeneity range of study-level rates: %.1f%% – %.1f%%\n",
            rate_min, rate_max))
cat(sprintf("Absolute spread: %.1f percentage points\n", rate_max - rate_min))

# Rate clustering for studies with N >= 20
cat("\nRate clustering (studies with N ≥ 20):\n")
large_studies <- study_summary %>% filter(N >= 20) %>% arrange(Rate_pct)
for (i in seq_len(nrow(large_studies))) {
  s <- large_studies[i, ]
  cat(sprintf("  %-18s: N=%3d, Rate=%.1f%%\n", s$Study_ID, s$N, s$Rate_pct))
}
cat(sprintf("  Range (N≥20): %.1f%% – %.1f%%\n",
            min(large_studies$Rate_pct), max(large_studies$Rate_pct)))

# ── 7. AlAbdi_2023 DESCRIPTIVE REPORT ────────────────────────────────────────
# NOTE: AlAbdi is NOT in this dataset (Group A only).
# We compute its metrics using the reported values from the protocol.
# If AlAbdi data were available, it would be loaded separately.
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  AlAbdi_2023 — DESCRIPTIVE (reported values, not in dataset)\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# AlAbdi reported: N=34, 13 rescued (headline), 6 true LRS-technical rescue
alabdi_n     <- 34
alabdi_resc  <- 13
alabdi_tech  <- 6

ci_alabdi <- wilson_ci(alabdi_resc, alabdi_n)
cat(sprintf("AlAbdi_2023: %d / %d = %.1f%%\n", alabdi_resc, alabdi_n,
            100 * alabdi_resc / alabdi_n))
cat(sprintf("  95%% Wilson CI: [%.1f%%, %.1f%%]\n",
            100 * ci_alabdi["lower"], 100 * ci_alabdi["upper"]))
cat(sprintf("  Exact values:  [%.6f%%, %.6f%%]\n",
            100 * ci_alabdi["lower"], 100 * ci_alabdi["upper"]))

ci_alabdi_tech <- wilson_ci(alabdi_tech, alabdi_n)
cat(sprintf("\nLRS-Technical Rescue: %d / %d = %.1f%%\n", alabdi_tech, alabdi_n,
            100 * alabdi_tech / alabdi_n))
cat(sprintf("  95%% Wilson CI: [%.1f%%, %.1f%%]\n",
            100 * ci_alabdi_tech["lower"], 100 * ci_alabdi_tech["upper"]))

# Fold enrichment: AlAbdi headline vs Group A
fold_enrichment <- (100 * alabdi_resc / alabdi_n) / baseline_rate
cat(sprintf("\nFold-enrichment (AlAbdi headline vs Group A): %.1f / %.1f = %.1f-fold\n",
            100 * alabdi_resc / alabdi_n, baseline_rate, fold_enrichment))

# ── 8. FULL STUDY-LEVEL SUMMARY TABLE ────────────────────────────────────────
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  COMPLETE STUDY-LEVEL SUMMARY TABLE\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

final_table <- study_summary %>%
  mutate(
    Weight_pct = round(100 * N / N_total, 1),
    Rate_fmt   = sprintf("%.1f%%", Rate_pct),
    CI_fmt     = sprintf("[%.1f%%, %.1f%%]", Wilson_lo, Wilson_hi),
    Def_fmt    = sprintf("%.1f%%", Def_rate),
    Def_CI_fmt = sprintf("[%.1f%%, %.1f%%]", Def_lo, Def_hi)
  ) %>%
  select(Study_ID, N, Def, Poss, Combined, Not_Resc, Weight_pct,
         Rate_fmt, CI_fmt, Def_fmt, Def_CI_fmt) %>%
  arrange(desc(N))

print(as.data.frame(final_table), right = FALSE)

# ── 10. EXACT (UNROUNDED) VALUES FOR MANUSCRIPT COMPARISON ───────────────────
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  EXACT VALUES (UNROUNDED) FOR MANUSCRIPT CROSS-CHECK\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

cat("── Group A Overall ──\n")
cat(sprintf("  Combined rate:    %d/%d = %.10f\n", n_comb, N_total, n_comb / N_total))
cat(sprintf("  Wilson CI lower:  %.10f\n", ci_comb["lower"]))
cat(sprintf("  Wilson CI upper:  %.10f\n", ci_comb["upper"]))
cat(sprintf("  Definitive rate:  %d/%d = %.10f\n", n_def, N_total, n_def / N_total))
cat(sprintf("  Wilson CI lower:  %.10f\n", ci_def["lower"]))
cat(sprintf("  Wilson CI upper:  %.10f\n", ci_def["upper"]))
cat(sprintf("  Possible rate:    %d/%d = %.10f\n", n_poss, N_total, n_poss / N_total))
cat(sprintf("  Wilson CI lower:  %.10f\n", ci_poss["lower"]))
cat(sprintf("  Wilson CI upper:  %.10f\n", ci_poss["upper"]))

cat("\n── AlAbdi_2023 ──\n")
cat(sprintf("  Headline rate:    %d/%d = %.10f\n", alabdi_resc, alabdi_n, alabdi_resc / alabdi_n))
cat(sprintf("  Wilson CI lower:  %.10f\n", ci_alabdi["lower"]))
cat(sprintf("  Wilson CI upper:  %.10f\n", ci_alabdi["upper"]))

cat("\n── Study-Level Exact CIs ──\n")
for (i in seq_len(nrow(study_summary))) {
  s <- study_summary[i, ]
  cat(sprintf("  %-18s: rate=%.10f  CI=[%.10f, %.10f]\n",
              s$Study_ID, s$Rate_pct / 100, s$Wilson_lo / 100, s$Wilson_hi / 100))
}

cat("\n── LOO Exact Values ──\n")
for (i in seq_len(nrow(loo_results))) {
  r <- loo_results[i, ]
  cat(sprintf("  Drop %-18s: rate=%.10f  Δ=%+.10f  CI=[%.10f, %.10f]\n",
              r$Study_ID, r$LOO_Rate_pct / 100, r$Change_pp / 100,
              r$LOO_Wilson_lo / 100, r$LOO_Wilson_hi / 100))
}

cat("\n═══════════════════════════════════════════════════════════════════\n")
cat("  STEP 1 VERIFICATION COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat(sprintf("  Internal consistency check: %s\n",
            ifelse(n_comb + n_notr + n_miss + n_conc + n_refi == N_total &
                   na_outcome == 0 & boundary_ok,
                   "ALL CHECKS PASSED", "ISSUES DETECTED")))
cat("═══════════════════════════════════════════════════════════════════\n")

# ── Close output file ────────────────────────────────────────────────────────
sink()
cat(sprintf("\nOutput saved to: %s\n", output_file))
