###############################################################################
# STEP 7: STUDY-LEVEL META-ANALYSIS — R VERIFICATION SCRIPT (v2)
# Systematic Review: LRS Rescue in SRS-Negative Patients
# Input: Dataset_S1.xlsx, sheet "Tier 1 - Group A" (11 studies, Group A SRS-negative patients)
#
# FIXES from v1:
#   - Back-transformation: metaprop with sm="PFT" stores TE.random on the
#     arcsine scale. Must use summary() to get back-transformed proportions.
#   - Steyaert_2025: excluded from prior testing subgroup (77/88 ambiguous ES_or_GS)
#   - Jang_2025: classified as Mixed in platform subgroup (excluded from test)
#   - Reports both FE and RE Q-between for subgroups
#
# Required R packages: meta, metafor
###############################################################################

# ── 0. LIBRARIES ─────────────────────────────────────────────────────────────
for (pkg in c("readxl", "dplyr", "meta", "metafor")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(readxl)
library(dplyr)
library(meta)
library(metafor)

# ── 1a. SAVE OUTPUT TO TXT FILE ─────────────────────────────────────────────
output_file <- "Step7_Meta_Analysis_OUTPUT.txt"
sink(output_file, split = TRUE)

# ── 1. DATA LOADING ──────────────────────────────────────────────────────────
dat <- read_excel(
  "Dataset_S1.xlsx",
  sheet     = "Tier 1 - Group A",
  skip      = 1,
  col_names = TRUE
)

cat("═══════════════════════════════════════════════════════════════════════\n")
cat("  STEP 7 — STUDY-LEVEL META-ANALYSIS: R VERIFICATION (v2)\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

# ── 2. RESCUE CLASSIFICATION ─────────────────────────────────────────────────
dat <- dat %>%
  mutate(
    Rescued = as.integer(
      !is.na(LRS_Outcome) &
      (LRS_Outcome == "Definitive_Rescue" | grepl("Possible_Rescue", LRS_Outcome))
    ),
    Def_Rescued = as.integer(!is.na(LRS_Outcome) & LRS_Outcome == "Definitive_Rescue")
  )

# ── 3. WILSON CI FUNCTION ────────────────────────────────────────────────────
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

# ── 4. HELPER: Extract back-transformed values from metaprop ─────────────────
get_bt <- function(mp_obj) {
  # metaprop with sm="PFT": all estimates are on Freeman-Tukey arcsine scale
  # Back-transform to proportion scale using sin²(theta)
  # Note: meta's internal back-transform uses a harmonic-mean-N correction
  # for the pooled estimate. We use sin²() which is the standard inverse.

  bt <- function(x) sin(x)^2  # Freeman-Tukey inverse: p = sin²(theta)

  # Clamp arcsine value to [0, pi/2] BEFORE sin², then clamp result to [0,1]
  # Negative arcsine values (PI extending below 0%) should map to 0%
  # Values above pi/2 (PI extending above 100%) should map to 100%
  bt_clamp <- function(x) {
    x_clamped <- max(0, min(pi/2, x))
    return(sin(x_clamped)^2)
  }

  pooled   <- bt_clamp(mp_obj$TE.random)
  ci_lower <- bt_clamp(mp_obj$lower.random)
  ci_upper <- bt_clamp(mp_obj$upper.random)

  # Prediction interval
  pi_lower <- tryCatch(bt_clamp(mp_obj$lower.predict), error = function(e) NA_real_)
  pi_upper <- tryCatch(bt_clamp(mp_obj$upper.predict), error = function(e) NA_real_)

  # I² — atomic in meta 8.x
  i2_val <- mp_obj$I2
  if (!is.na(i2_val) && i2_val > 1) i2_val <- i2_val / 100

  list(
    pooled    = pooled,
    ci_lower  = ci_lower,
    ci_upper  = ci_upper,
    pi_lower  = pi_lower,
    pi_upper  = pi_upper,
    I2        = i2_val,
    tau2      = mp_obj$tau2,    # reported on arcsine scale (standard for FT)
    tau       = mp_obj$tau,
    Q         = mp_obj$Q,
    df_Q      = mp_obj$df.Q,
    pval_Q    = mp_obj$pval.Q
  )
}

# ── 5. BUILD STUDY-LEVEL TABLE ───────────────────────────────────────────────
study_data <- dat %>%
  group_by(Study_ID) %>%
  summarise(
    N        = n(),
    Rescued  = sum(Rescued),
    Def_Only = sum(Def_Rescued),
    .groups  = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    Rate     = 100 * Rescued / N,
    W_lo     = 100 * wilson_ci(Rescued, N)["lower"],
    W_hi     = 100 * wilson_ci(Rescued, N)["upper"],
    Def_Rate = 100 * Def_Only / N
  ) %>%
  ungroup() %>%
  arrange(desc(N))

cat("Study-Level Data:\n")
cat("──────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-18s %4s %4s %4s %8s  [%s]\n",
            "Study", "N", "Resc", "Def", "Rate%", "95% Wilson CI"))
cat("──────────────────────────────────────────────────────────────────\n")
for (i in seq_len(nrow(study_data))) {
  s <- study_data[i, ]
  cat(sprintf("%-18s %4d %4d %4d %7.1f%%  [%5.1f%%, %5.1f%%]\n",
              s$Study_ID, s$N, s$Rescued, s$Def_Only,
              s$Rate, s$W_lo, s$W_hi))
}
cat(sprintf("%-18s %4d %4d %4d %7.1f%%\n", "TOTAL",
            sum(study_data$N), sum(study_data$Rescued), sum(study_data$Def_Only),
            100 * sum(study_data$Rescued) / sum(study_data$N)))

# Subgroup assignments
study_data <- study_data %>%
  mutate(
    Prior_Subgroup = case_when(
      Study_ID %in% c("Hiatt_2024", "Lesurf_2025", "Negi_2025", "Hiatt_2021") ~ "WGS-only",
      Study_ID %in% c("Cohen_2022", "Jang_2025", "Pauper_2021", "Redfield_2024") ~ "WES+WGS",
      Study_ID %in% c("Sinha_2025", "Fukuda_2023") ~ "WES-only",
      Study_ID == "Steyaert_2025" ~ NA_character_,  # excluded: 77/88 ambiguous ES_or_GS
      TRUE ~ NA_character_
    ),
    Platform_Subgroup = case_when(
      Study_ID %in% c("Fukuda_2023", "Negi_2025", "Sinha_2025") ~ "ONT",
      Study_ID == "Jang_2025" ~ "Mixed",
      TRUE ~ "PacBio"
    )
  )

cat("\n\nSubgroup assignments:\n")
for (i in seq_len(nrow(study_data))) {
  s <- study_data[i, ]
  prior_label <- ifelse(is.na(s$Prior_Subgroup),
                        "EXCLUDED (77/88 ambiguous ES_or_GS)", s$Prior_Subgroup)
  jang_note <- ifelse(s$Study_ID == "Jang_2025", " *PacBio 17 + ONT 12", "")
  cat(sprintf("  %-18s: Prior=%-16s Platform=%s%s\n",
              s$Study_ID, prior_label, s$Platform_Subgroup, jang_note))
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: PRIMARY META-ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  SECTION 1: PRIMARY META-ANALYSIS (Combined Rescue)\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

mp <- metaprop(
  event = study_data$Rescued, n = study_data$N, studlab = study_data$Study_ID,
  sm = "PFT", method.tau = "REML", method.ci = "WS", prediction = TRUE
)
bt <- get_bt(mp)

cat(sprintf("Pooled estimate (back-transformed): %.4f (%.2f%%)\n", bt$pooled, 100*bt$pooled))
cat(sprintf("95%% CI:  [%.4f, %.4f] ([%.2f%%, %.2f%%])\n",
            bt$ci_lower, bt$ci_upper, 100*bt$ci_lower, 100*bt$ci_upper))
cat(sprintf("95%% PI:  [%.4f, %.4f] ([%.2f%%, %.2f%%])\n",
            bt$pi_lower, bt$pi_upper, 100*bt$pi_lower, 100*bt$pi_upper))
cat(sprintf("Simple pooled: %d/%d = %.2f%%\n",
            sum(study_data$Rescued), sum(study_data$N),
            100*sum(study_data$Rescued)/sum(study_data$N)))
cat(sprintf("Arcsine-scale TE.random: %.6f\n", mp$TE.random))

cat(sprintf("\nHeterogeneity:\n"))
cat(sprintf("  Q=%.4f (df=%d, p=%.6f) | I²=%.1f%% | τ²=%.6f (arcsine scale) | τ=%.6f\n",
            bt$Q, bt$df_Q, bt$pval_Q, bt$I2*100, bt$tau2, bt$tau))
cat("  NOTE: τ² is on the Freeman-Tukey arcsine scale (standard for sm='PFT').\n")
cat("  It cannot be directly interpreted as a variance in proportions.\n")

i2_label <- ifelse(bt$I2*100<=25,"Low",ifelse(bt$I2*100<=50,"Moderate",
            ifelse(bt$I2*100<=75,"Substantial","Considerable")))
cat(sprintf("  Interpretation: %s\n", i2_label))

cat("\n── Study Weights ──\n")
w <- mp$w.random; w_pct <- 100*w/sum(w)
for (i in seq_len(nrow(study_data)))
  cat(sprintf("  %-18s: %.1f%%\n", study_data$Study_ID[i], w_pct[i]))


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2: LEAVE-ONE-OUT
# ══════════════════════════════════════════════════════════════════════════════
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  SECTION 2: LEAVE-ONE-OUT SENSITIVITY\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

cat(sprintf("Baseline: %.2f%%\n\n", 100*bt$pooled))
cat(sprintf("%-18s %8s %10s %10s %8s %8s %10s\n",
            "Dropped","Pooled%","CI_lo","CI_hi","Δ(pp)","I²%","τ²"))
cat("──────────────────────────────────────────────────────────────────────────────\n")

for (i in seq_len(nrow(study_data))) {
  idx <- seq_len(nrow(study_data))[-i]
  loo <- metaprop(event=study_data$Rescued[idx], n=study_data$N[idx],
                  studlab=study_data$Study_ID[idx], sm="PFT",
                  method.tau="REML", method.ci="WS", prediction=TRUE)
  lb <- get_bt(loo)
  cat(sprintf("%-18s %7.2f%% %9.2f%% %9.2f%% %+7.2fpp %7.1f%% %10.6f\n",
              study_data$Study_ID[i], 100*lb$pooled, 100*lb$ci_lower, 100*lb$ci_upper,
              100*(lb$pooled - bt$pooled), lb$I2*100, lb$tau2))
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3: SENSITIVITY — N≥20
# ══════════════════════════════════════════════════════════════════════════════
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  SECTION 3: SENSITIVITY — N≥20 RESTRICTION\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

lg <- study_data %>% filter(N >= 20)
cat(sprintf("Included: %d | Excluded: %s\n\n",
            nrow(lg), paste(study_data$Study_ID[study_data$N<20], collapse=", ")))

mp_lg <- metaprop(event=lg$Rescued, n=lg$N, studlab=lg$Study_ID,
                  sm="PFT", method.tau="REML", method.ci="WS", prediction=TRUE)
bt_lg <- get_bt(mp_lg)
cat(sprintf("Pooled: %.2f%% [%.2f%%, %.2f%%]\n", 100*bt_lg$pooled, 100*bt_lg$ci_lower, 100*bt_lg$ci_upper))
cat(sprintf("PI: [%.2f%%, %.2f%%] | I²=%.1f%% | τ²=%.6f\n",
            100*bt_lg$pi_lower, 100*bt_lg$pi_upper, bt_lg$I2*100, bt_lg$tau2))

cat("\n── Protocol-specified exclusion (Hiatt_2021, Fukuda_2023, Pauper_2021) ──\n\n")
nf <- study_data %>% filter(!Study_ID %in% c("Hiatt_2021","Fukuda_2023","Pauper_2021"))
mp_nf <- metaprop(event=nf$Rescued, n=nf$N, studlab=nf$Study_ID,
                  sm="PFT", method.tau="REML", method.ci="WS", prediction=TRUE)
bt_nf <- get_bt(mp_nf)
cat(sprintf("Pooled: %.2f%% [%.2f%%, %.2f%%]\n", 100*bt_nf$pooled, 100*bt_nf$ci_lower, 100*bt_nf$ci_upper))
cat(sprintf("PI: [%.2f%%, %.2f%%] | I²=%.1f%% | τ²=%.6f\n",
            100*bt_nf$pi_lower, 100*bt_nf$pi_upper, bt_nf$I2*100, bt_nf$tau2))


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4: SENSITIVITY — DEFINITIVE ONLY
# ══════════════════════════════════════════════════════════════════════════════
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  SECTION 4: SENSITIVITY — DEFINITIVE RESCUE ONLY\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

mp_def <- metaprop(event=study_data$Def_Only, n=study_data$N, studlab=study_data$Study_ID,
                   sm="PFT", method.tau="REML", method.ci="WS", prediction=TRUE)
bt_def <- get_bt(mp_def)
cat(sprintf("Pooled: %.2f%% [%.2f%%, %.2f%%]\n", 100*bt_def$pooled, 100*bt_def$ci_lower, 100*bt_def$ci_upper))
cat(sprintf("PI: [%.2f%%, %.2f%%] | I²=%.1f%% | τ²=%.6f\n",
            100*bt_def$pi_lower, 100*bt_def$pi_upper, bt_def$I2*100, bt_def$tau2))
cat(sprintf("Combined %.2f%% vs Def-only %.2f%% (Δ=%.2fpp)\n",
            100*bt$pooled, 100*bt_def$pooled, 100*(bt$pooled - bt_def$pooled)))


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5: PUBLICATION BIAS
# ══════════════════════════════════════════════════════════════════════════════
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  SECTION 5: PUBLICATION BIAS — EGGER'S TEST\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

egger <- metabias(mp, method.bias = "linreg")
cat(sprintf("Intercept: %.4f | t=%.4f | p=%.6f\n", egger$estimate[1], egger$statistic, egger$p.value))
cat(sprintf("Interpretation: %s\n",
            ifelse(egger$p.value < 0.05, "Significant asymmetry", "No significant asymmetry")))


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6: SUBGROUP — PRIOR TESTING
# ══════════════════════════════════════════════════════════════════════════════
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  SECTION 6: SUBGROUP — PRIOR TESTING TYPE\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

prior_sub <- study_data %>% filter(!is.na(Prior_Subgroup))
cat(sprintf("Studies in subgroup analysis: %d (Steyaert_2025 excluded: 77/88 ambiguous ES_or_GS)\n\n", nrow(prior_sub)))
cat("Composition:\n")
for (sg in unique(prior_sub$Prior_Subgroup)) {
  sids <- prior_sub$Study_ID[prior_sub$Prior_Subgroup == sg]
  cat(sprintf("  %s (%d): %s\n", sg, length(sids), paste(sids, collapse=", ")))
}

invisible(mp_prior <- metaprop(event=prior_sub$Rescued, n=prior_sub$N, studlab=prior_sub$Study_ID,
                     subgroup=prior_sub$Prior_Subgroup, sm="PFT",
                     method.tau="REML", method.ci="WS", prediction=TRUE))

cat("\n── Subgroup Estimates ──\n\n")
for (sg in unique(prior_sub$Prior_Subgroup)) {
  sg_idx <- which(prior_sub$Prior_Subgroup == sg)
  if (length(sg_idx) >= 3) {
    sg_mp <- metaprop(event=prior_sub$Rescued[sg_idx], n=prior_sub$N[sg_idx],
                      studlab=prior_sub$Study_ID[sg_idx], sm="PFT",
                      method.tau="REML", method.ci="WS", prediction=TRUE)
    sg_bt <- get_bt(sg_mp)
    cat(sprintf("  %s (%d): %.2f%% [%.2f%%, %.2f%%] | I²=%.1f%%\n",
                sg, length(sg_idx), 100*sg_bt$pooled, 100*sg_bt$ci_lower,
                100*sg_bt$ci_upper, sg_bt$I2*100))
  } else {
    cat(sprintf("  %s (%d): DESCRIPTIVE (<3 studies)\n", sg, length(sg_idx)))
    for (j in sg_idx)
      cat(sprintf("    %s: %d/%d = %.1f%%\n",
                  prior_sub$Study_ID[j], prior_sub$Rescued[j], prior_sub$N[j], prior_sub$Rate[j]))
  }
}

cat(sprintf("\nQ-between (FE):  %.4f (df=%d, p=%.6f)\n",
            mp_prior$Q.b.fixed, mp_prior$df.Q.b, mp_prior$pval.Q.b.fixed))
cat(sprintf("Q-between (RE):  %.4f (df=%d, p=%.6f)\n",
            mp_prior$Q.b.random, mp_prior$df.Q.b, mp_prior$pval.Q.b.random))


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 7: SUBGROUP — PLATFORM
# ══════════════════════════════════════════════════════════════════════════════
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  SECTION 7: SUBGROUP — PLATFORM\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

plat_sub <- study_data %>% filter(Platform_Subgroup != "Mixed")
cat("Composition (Jang excluded as Mixed):\n")
for (sg in c("ONT","PacBio")) {
  sids <- plat_sub$Study_ID[plat_sub$Platform_Subgroup == sg]
  cat(sprintf("  %s (%d): %s\n", sg, length(sids), paste(sids, collapse=", ")))
}

invisible(mp_plat <- metaprop(event=plat_sub$Rescued, n=plat_sub$N, studlab=plat_sub$Study_ID,
                    subgroup=plat_sub$Platform_Subgroup, sm="PFT",
                    method.tau="REML", method.ci="WS", prediction=TRUE))

cat("\n── Platform Estimates ──\n\n")
for (sg in c("ONT","PacBio")) {
  sg_idx <- which(plat_sub$Platform_Subgroup == sg)
  if (length(sg_idx) >= 3) {
    sg_mp <- metaprop(event=plat_sub$Rescued[sg_idx], n=plat_sub$N[sg_idx],
                      studlab=plat_sub$Study_ID[sg_idx], sm="PFT",
                      method.tau="REML", method.ci="WS", prediction=TRUE)
    sg_bt <- get_bt(sg_mp)
    cat(sprintf("  %s (%d): %.2f%% [%.2f%%, %.2f%%] | I²=%.1f%%\n",
                sg, length(sg_idx), 100*sg_bt$pooled, 100*sg_bt$ci_lower,
                100*sg_bt$ci_upper, sg_bt$I2*100))
  }
}

cat(sprintf("\nQ-between (FE):  %.4f (df=%d, p=%.6f)\n",
            mp_plat$Q.b.fixed, mp_plat$df.Q.b, mp_plat$pval.Q.b.fixed))
cat(sprintf("Q-between (RE):  %.4f (df=%d, p=%.6f)\n",
            mp_plat$Q.b.random, mp_plat$df.Q.b, mp_plat$pval.Q.b.random))


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 8: DIAGNOSTICS & EXACT VALUES
# ══════════════════════════════════════════════════════════════════════════════
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  SECTION 8: DIAGNOSTICS\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

checks <- c(
  total_N    = sum(study_data$N)==nrow(dat), total_resc = sum(study_data$Rescued)==55,
  total_def  = sum(study_data$Def_Only)==36, n_studies = nrow(study_data)==11,
  pooled_ok  = bt$pooled>=0 & bt$pooled<=1,
  ci_ok      = bt$ci_lower>=0 & bt$ci_upper<=1,
  i2_ok      = bt$I2>=0 & bt$I2<=1
)
for (nm in names(checks))
  cat(sprintf("  %-15s: %s\n", nm, ifelse(checks[nm], "PASS", "FAIL")))
cat(sprintf("\nOverall: %s\n", ifelse(all(checks), "ALL PASSED", "ISSUES")))

cat("\n── Exact Values ──\n")
cat(sprintf("  Pooled (BT):  %.15f\n", bt$pooled))
cat(sprintf("  CI:           [%.15f, %.15f]\n", bt$ci_lower, bt$ci_upper))
cat(sprintf("  PI:           [%.15f, %.15f]\n", bt$pi_lower, bt$pi_upper))
cat(sprintf("  I²:           %.15f\n", bt$I2))
cat(sprintf("  τ²:           %.15f\n", bt$tau2))
cat(sprintf("  Q:            %.15f (p=%.15f)\n", bt$Q, bt$pval_Q))
cat(sprintf("  Egger p:      %.15f\n", egger$p.value))

cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat("  STEP 7 VERIFICATION COMPLETE (v2)\n")
cat("═══════════════════════════════════════════════════════════════════════\n")

sink()
cat(sprintf("\nOutput saved to: %s\n", output_file))
