###############################################################################
# STEP 2 — RESCUE BY PRIOR TESTING TYPE: ANALYSIS FUNCTIONS  (v3)
#
# Refactored so Figure 4 (and any other downstream script) can:
#   source("Step2_Prior_Testing_v2.R")
#   res <- run_prior_testing_analysis("Dataset_S1.xlsx")
#
# Standalone use prints a full text report:
#   Rscript Step2_Prior_Testing_v2.R
#
# Returns a named list with:
#   $N_total, $N_class, $n_ambiguous, $rescued_total, $rescued_class
#   $cat_results          one-vs-rest Fisher results (3 categories)
#   $panel_only           descriptive Panel-only row
#   $pairwise             3 pre-specified pairwise contrasts
#   $depth_summary        depth-stratum rates with Wilson CIs
#   $ca_chisq, $ca_chisq_p   prop.trend.test
#   $ca_Z, $ca_Z_p           DescTools CochranArmitage (NA if not installed)
#   $platform             ONT vs PacBio
#   $ambiguous_breakdown  per-study ambiguity counts
###############################################################################

# ── 0. LIBRARIES ─────────────────────────────────────────────────────────────
.required <- c("readxl", "dplyr", "tidyr")
for (pkg in .required) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
if (!requireNamespace("DescTools", quietly = TRUE)) {
  message("DescTools not installed — Cochran-Armitage Z cross-check will be skipped.")
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
})
HAVE_DESCTOOLS <- requireNamespace("DescTools", quietly = TRUE)
if (HAVE_DESCTOOLS) suppressPackageStartupMessages(library(DescTools))

# ── 1. WILSON CI HELPER ─────────────────────────────────────────────────────
wilson_ci <- function(x, n, conf.level = 0.95) {
  if (n == 0) return(c(lower = NA_real_, upper = NA_real_))
  z      <- qnorm(1 - (1 - conf.level) / 2)
  p      <- x / n
  denom  <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  margin <- (z / denom) * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  c(lower = max(0, center - margin), upper = min(1, center + margin))
}

# ── 2. MAIN ANALYSIS FUNCTION ────────────────────────────────────────────────
run_prior_testing_analysis <- function(input_file  = "Dataset_S1.xlsx",
                                       input_sheet = "Tier 1 - Group A") {

  # ── Load with two-row header handling ──────────────────────────────────────
  raw <- read_excel(input_file, sheet = input_sheet,
                    col_names = FALSE, .name_repair = "minimal")
  header_row <- as.character(unlist(raw[2, ]))
  dat <- raw[-(1:2), ]
  colnames(dat) <- header_row

  panel_col <- if ("Prior_GenePanel" %in% colnames(dat)) {
    "Prior_GenePanel"
  } else if ("Prior_NGS_GenePanel" %in% colnames(dat)) {
    "Prior_NGS_GenePanel"
  } else {
    stop("Cannot find a Prior_GenePanel / Prior_NGS_GenePanel column.")
  }

  required_cols <- c("Study_ID", "LRS_Outcome", "Prior_WES", "Prior_WGS",
                     "Previous_Tests(or simultaneously)", "Platform")
  missing <- setdiff(required_cols, colnames(dat))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }

  N_total <- nrow(dat)

  # ── Rescue classification (asterisk-tolerant) ──────────────────────────────
  norm_outcome <- function(x) gsub("\\*+$", "", tolower(trimws(as.character(x))))
  dat <- dat %>%
    mutate(
      .outcome_norm = norm_outcome(LRS_Outcome),
      Def_Rescued   = as.integer(.outcome_norm == "definitive_rescue"),
      Poss_Rescued  = as.integer(.outcome_norm == "possible_rescue"),
      Rescued       = pmax(Def_Rescued, Poss_Rescued),
      Not_Rescued   = as.integer(.outcome_norm == "not_rescued")
    )

  # ── Prior-testing flags ────────────────────────────────────────────────────
  dat <- dat %>%
    mutate(
      wes_yes   = tolower(trimws(as.character(Prior_WES))) == "yes",
      wgs_yes   = tolower(trimws(as.character(Prior_WGS))) == "yes",
      panel_yes = tolower(trimws(as.character(.data[[panel_col]]))) == "yes",
      wes_unrep = tolower(trimws(as.character(Prior_WES))) == "not_reported",
      wgs_unrep = tolower(trimws(as.character(Prior_WGS))) == "not_reported"
    )

  # ── Ambiguous-prior-testing definition (CORRECTED) ─────────────────────────
  # Steyaert 77: Previous_Tests = "ES_or_GS"
  # Negi 8:      both Prior_WES and Prior_WGS = Not_reported
  dat <- dat %>%
    mutate(
      is_ES_or_GS = !is.na(`Previous_Tests(or simultaneously)`) &
                    `Previous_Tests(or simultaneously)` == "ES_or_GS",
      is_negi_ambiguous = (Study_ID == "Negi_2025" & wes_unrep & wgs_unrep),
      is_ambiguous = is_ES_or_GS | is_negi_ambiguous
    )

  ambiguous_breakdown <- dat %>%
    filter(is_ambiguous) %>%
    mutate(amb_type = ifelse(is_ES_or_GS, "Steyaert ES_or_GS", "Negi unreported")) %>%
    count(Study_ID, amb_type, name = "N") %>%
    arrange(desc(N))  # larger group first in figure caption

  # ── Categorize ─────────────────────────────────────────────────────────────
  dat <- dat %>%
    mutate(
      Prior_Category = case_when(
        is_ambiguous                     ~ "Ambiguous",
        wes_yes &  wgs_yes               ~ "WES+WGS",
        wes_yes & !wgs_yes               ~ "WES-only",
       !wes_yes &  wgs_yes               ~ "WGS-only",
        panel_yes & !wes_yes & !wgs_yes  ~ "Panel-only",
        TRUE                             ~ "Other/Unknown"
      ),
      Testing_Depth = ifelse(
        is_ambiguous, NA_integer_,
        as.integer(wes_yes) + as.integer(wgs_yes) + as.integer(panel_yes)
      )
    )

  dat_class       <- dat %>% filter(!is_ambiguous)
  N_class         <- nrow(dat_class)
  n_class_rescued <- sum(dat_class$Rescued)

  # ── One-vs-rest Fisher's exact (Panel-only excluded as N=1) ────────────────
  class_cats_test <- dat_class %>%
    filter(Prior_Category != "Panel-only") %>%
    group_by(Prior_Category) %>%
    summarise(N = n(), Rescued = sum(Rescued), .groups = "drop")

  N_test         <- sum(class_cats_test$N)
  n_rescued_test <- sum(class_cats_test$Rescued)

  ovr <- data.frame(
    Category = character(), N = integer(), Rescued = integer(),
    Rate_pct = numeric(), CI_lo = numeric(), CI_hi = numeric(),
    OR = numeric(), OR_lo = numeric(), OR_hi = numeric(),
    p_raw = numeric(), stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(class_cats_test))) {
    cn <- class_cats_test$Prior_Category[i]
    nc <- class_cats_test$N[i]
    rc <- class_cats_test$Rescued[i]
    nr <- N_test - nc
    rr <- n_rescued_test - rc
    mat <- matrix(c(rc, nc - rc, rr, nr - rr), nrow = 2, byrow = TRUE)
    ft  <- fisher.test(mat)
    wci <- wilson_ci(rc, nc)

    ovr <- rbind(ovr, data.frame(
      Category = cn, N = nc, Rescued = rc,
      Rate_pct = 100 * rc / nc,
      CI_lo = 100 * wci["lower"], CI_hi = 100 * wci["upper"],
      OR    = unname(ft$estimate),
      OR_lo = ft$conf.int[1], OR_hi = ft$conf.int[2],
      p_raw = ft$p.value,
      stringsAsFactors = FALSE
    ))
  }
  ovr$p_BH <- p.adjust(ovr$p_raw, method = "BH")
  ovr <- ovr %>% arrange(p_raw)
  rownames(ovr) <- NULL  # strip 'lower'/'lower1'/'lower2' inherited from wilson_ci()

  # ── Panel-only descriptive ─────────────────────────────────────────────────
  pn_dat <- dat_class %>% filter(Prior_Category == "Panel-only")
  panel_only <- if (nrow(pn_dat) > 0) {
    pn <- nrow(pn_dat); pr <- sum(pn_dat$Rescued)
    pci <- wilson_ci(pr, pn)
    df_po <- data.frame(
      Category = "Panel-only", N = pn, Rescued = pr,
      Rate_pct = 100 * pr / pn,
      CI_lo = 100 * pci["lower"], CI_hi = 100 * pci["upper"],
      stringsAsFactors = FALSE
    )
    rownames(df_po) <- NULL
    df_po
  } else NULL

  # ── Pre-specified pairwise contrasts ───────────────────────────────────────
  run_pairwise <- function(cat1, cat2) {
    sub <- dat_class %>% filter(Prior_Category %in% c(cat1, cat2))
    n1 <- sum(sub$Prior_Category == cat1)
    r1 <- sum(sub$Rescued[sub$Prior_Category == cat1])
    n2 <- sum(sub$Prior_Category == cat2)
    r2 <- sum(sub$Rescued[sub$Prior_Category == cat2])
    mat <- matrix(c(r1, n1 - r1, r2, n2 - r2), nrow = 2, byrow = TRUE)
    ft  <- fisher.test(mat)
    data.frame(
      Comparison = paste(cat1, "vs", cat2),
      cat1 = cat1, cat2 = cat2,
      n1 = n1, r1 = r1, n2 = n2, r2 = r2,
      OR    = unname(ft$estimate),
      OR_lo = ft$conf.int[1], OR_hi = ft$conf.int[2],
      p     = ft$p.value, stringsAsFactors = FALSE
    )
  }
  pairwise <- rbind(
    run_pairwise("WES-only", "WGS-only"),
    run_pairwise("WES-only", "WES+WGS"),
    run_pairwise("WGS-only", "WES+WGS")
  )

  # ── Cochran-Armitage trend test ────────────────────────────────────────────
  depth_tab <- dat_class %>%
    filter(!is.na(Testing_Depth)) %>%
    group_by(Testing_Depth) %>%
    summarise(rescued = sum(Rescued),
              not_rescued = sum(1 - Rescued), .groups = "drop") %>%
    arrange(Testing_Depth)

  ca_test    <- prop.trend.test(depth_tab$rescued,
                                depth_tab$rescued + depth_tab$not_rescued,
                                depth_tab$Testing_Depth)
  ca_chisq   <- as.numeric(ca_test$statistic)
  ca_chisq_p <- ca_test$p.value

  if (HAVE_DESCTOOLS) {
    ca_mat <- as.matrix(depth_tab[, c("rescued", "not_rescued")])
    rownames(ca_mat) <- depth_tab$Testing_Depth
    ca_desc <- DescTools::CochranArmitageTest(ca_mat)
    ca_Z   <- as.numeric(ca_desc$statistic)
    ca_Z_p <- ca_desc$p.value
  } else {
    ca_Z   <- NA_real_
    ca_Z_p <- NA_real_
  }

  depth_summary <- dat_class %>%
    group_by(Testing_Depth) %>%
    summarise(N = n(), Rescued = sum(Rescued), .groups = "drop") %>%
    rowwise() %>%
    mutate(
      Rate_pct  = 100 * Rescued / N,
      Wilson_lo = 100 * wilson_ci(Rescued, N)["lower"],
      Wilson_hi = 100 * wilson_ci(Rescued, N)["upper"]
    ) %>%
    ungroup() %>%
    arrange(Testing_Depth)

  # ── Platform analysis ──────────────────────────────────────────────────────
  plat <- dat %>%
    group_by(Platform) %>%
    summarise(N = n(), Def = sum(Def_Rescued),
              Poss = sum(Poss_Rescued), Rescued = sum(Rescued),
              .groups = "drop") %>%
    rowwise() %>%
    mutate(
      Rate_pct  = 100 * Rescued / N,
      Wilson_lo = 100 * wilson_ci(Rescued, N)["lower"],
      Wilson_hi = 100 * wilson_ci(Rescued, N)["upper"]
    ) %>%
    ungroup()

  # ── Return ─────────────────────────────────────────────────────────────────
  list(
    N_total           = N_total,
    N_class           = N_class,
    n_ambiguous       = sum(dat$is_ambiguous),
    rescued_total     = sum(dat$Rescued),
    rescued_class     = n_class_rescued,
    cat_results       = ovr,
    panel_only        = panel_only,
    pairwise          = pairwise,
    depth_summary     = depth_summary,
    ca_chisq          = ca_chisq,
    ca_chisq_p        = ca_chisq_p,
    ca_Z              = ca_Z,
    ca_Z_p            = ca_Z_p,
    platform          = plat,
    ambiguous_breakdown = ambiguous_breakdown,
    dat               = dat,
    dat_class         = dat_class
  )
}

# ── 3. STANDALONE REPORT (only when run via Rscript, not when sourced) ──────
if (sys.nframe() == 0L) {

  # ── SAVE OUTPUT TO TXT FILE (standalone runs only) ────────────────────────
  output_file <- "Step2_Prior_Testing_OUTPUT.txt"
  sink(output_file, split = TRUE)  # split=TRUE prints to console AND file

  res <- run_prior_testing_analysis()

  cat("═══════════════════════════════════════════════════════════════════════\n")
  cat("  STEP 2 — RESCUE BY PRIOR TESTING TYPE: R VERIFICATION (v3)\n")
  cat("═══════════════════════════════════════════════════════════════════════\n\n")

  cat(sprintf("Group A total:   N = %d\n", res$N_total))
  cat(sprintf("Ambiguous total: N = %d (excluded from prior-testing tests)\n",
              res$n_ambiguous))
  cat(sprintf("Classifiable:    N = %d\n",  res$N_class))
  cat(sprintf("Rescued (full):  %d\n",      res$rescued_total))
  cat(sprintf("Rescued (class): %d\n\n",    res$rescued_class))

  cat("── Ambiguous breakdown ──\n")
  for (i in seq_len(nrow(res$ambiguous_breakdown))) {
    a <- res$ambiguous_breakdown[i, ]
    cat(sprintf("  %s (%s): N = %d\n", a$Study_ID, a$amb_type, a$N))
  }
  cat("\n")

  cat("── One-vs-rest Fisher's exact (Panel A) ──\n")
  cat(sprintf("%-12s %5s %5s %8s %14s %10s %10s %10s\n",
              "Category", "N", "Resc", "Rate%", "Wilson CI",
              "OR", "p_raw", "p_BH"))
  for (i in seq_len(nrow(res$cat_results))) {
    r <- res$cat_results[i, ]
    cat(sprintf("%-12s %5d %5d %7.1f%% [%4.1f–%4.1f] %10.3f %10.6f %10.6f\n",
                r$Category, r$N, r$Rescued, r$Rate_pct,
                r$CI_lo, r$CI_hi, r$OR, r$p_raw, r$p_BH))
  }
  cat("\n")

  if (!is.null(res$panel_only)) {
    p <- res$panel_only
    cat(sprintf("Panel-only descriptive: %d/%d = %.1f%% [%.1f–%.1f] (excluded from formal tests)\n\n",
                p$Rescued, p$N, p$Rate_pct, p$CI_lo, p$CI_hi))
  }

  cat("── Pairwise contrasts (Panel B) ──\n")
  for (i in seq_len(nrow(res$pairwise))) {
    pw <- res$pairwise[i, ]
    cat(sprintf("  %-22s OR=%.3f [%.3f, %.3f]  p=%.4f\n",
                pw$Comparison, pw$OR, pw$OR_lo, pw$OR_hi, pw$p))
  }
  cat("\n")

  cat("── Testing-depth strata (Panel C) ──\n")
  for (i in seq_len(nrow(res$depth_summary))) {
    d <- res$depth_summary[i, ]
    cat(sprintf("  Depth %d: %d/%d = %.1f%% [%.1f–%.1f]\n",
                d$Testing_Depth, d$Rescued, d$N, d$Rate_pct,
                d$Wilson_lo, d$Wilson_hi))
  }
  cat(sprintf("\nCochran-Armitage:\n"))
  cat(sprintf("  prop.trend.test  : chi-sq = %.4f, p = %.4f\n",
              res$ca_chisq, res$ca_chisq_p))
  if (!is.na(res$ca_Z)) {
    cat(sprintf("  DescTools (Z)    : Z = %.4f, p = %.4f\n",
                res$ca_Z, res$ca_Z_p))
  }
  cat("\n")

  cat("── Platform comparison ──\n")
  for (i in seq_len(nrow(res$platform))) {
    p <- res$platform[i, ]
    cat(sprintf("  %-7s: %d/%d = %.1f%% [%.1f–%.1f]\n",
                p$Platform, p$Rescued, p$N, p$Rate_pct,
                p$Wilson_lo, p$Wilson_hi))
  }
  cat("\n═══════════════════════════════════════════════════════════════════════\n")

  # ── Close output file ──────────────────────────────────────────────────────
  sink()
  cat(sprintf("\nOutput saved to: %s\n", output_file))
}
