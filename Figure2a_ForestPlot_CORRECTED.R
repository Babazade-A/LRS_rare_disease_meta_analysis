###############################################################################
# FIGURE 2: Forest Plot — LRS Rescue Yield Meta-Analysis
#
# VISUAL:  forestploter — Study | WIDE CI graphic | Rate [CI] | Weight
# VALUES:  Hardcoded from verified Step 7 OUTPUT v2 (meta::metaprop + sin2 back-tx)
#
# All hardcoded values verified against Step7_Meta_Analysis_OUTPUT.txt (v2):
#   N=538, 55 events, k=11
#   Pooled 16.41% [6.87%, 29.03%], PI [0.00%, 66.78%]
#   I²=85.2%, τ²=0.052598, Q=67.48 (df=10), Egger p=0.0426
#   Sensitivity (protocol-specified exclusion of Hiatt_2021, Fukuda_2023,
#     Pauper_2021): 10.87% [5.62%, 17.59%], k=8, N=522
#
# RUN: Open in RStudio -> Source (Ctrl+Shift+S)
###############################################################################

output_dir <- tryCatch(
  dirname(rstudioapi::getSourceEditorContext()$path),
  error = function(e) getwd()
)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ── Packages ─────────────────────────────────────────────────────────────────
for (pkg in c("forestploter", "grid")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(forestploter)
library(grid)

# ── Study data (CORRECTED to match Step 7 v2 output) ────────────────────────
study_data <- data.frame(
  Study  = c("Hiatt 2021",  "Pauper 2021",  "Cohen 2022",
             "Fukuda 2023", "Hiatt 2024",   "Redfield 2024",
             "Jang 2025",   "Lesurf 2025",  "Negi 2025",
             "Sinha 2025",  "Steyaert 2025"),
  N      = c(6, 5, 176, 5, 96, 10, 29, 42, 31, 50, 88),
  Events = c(2, 0,   2, 5, 16,  1,  7,  2,  5,  4, 11),
  stringsAsFactors = FALSE
)
stopifnot(sum(study_data$N) == 538, sum(study_data$Events) == 55)

###############################################################################
# VERIFIED VALUES FROM STEP 7 OUTPUT v2 (meta::metaprop + sin2 back-transform)
###############################################################################
POOLED_EST <- 0.1641; POOLED_LO <- 0.0687; POOLED_HI <- 0.2903
PI_LO      <- 0.0000; PI_HI     <- 0.6678
I2_PCT     <- 85.2;   TAU2      <- 0.052598
Q_VAL      <- 67.48;  Q_DF      <- 10;     Q_P <- "< 0.001"
EGGER_P    <- 0.0426
# Sensitivity: protocol-specified exclusion (drop Hiatt_2021, Fukuda_2023, Pauper_2021)
SENS_EST   <- 0.1087; SENS_LO   <- 0.0562; SENS_HI <- 0.1759
SENS_K     <- 8;      SENS_N    <- 522
RE_WEIGHTS <- c(6.6, 6.1, 11.1, 6.1, 10.9, 7.9, 9.8, 10.3, 9.9, 10.4, 10.8)

# ── Wilson CI ────────────────────────────────────────────────────────────────
wilson_ci <- function(x, n) {
  if (n == 0) return(c(0, 0))
  z <- qnorm(0.975); p <- x/n; d <- 1 + z^2/n
  ctr <- (p + z^2/(2*n))/d
  mg <- (z/d)*sqrt(p*(1-p)/n + z^2/(4*n^2))
  c(max(0, ctr-mg), min(1, ctr+mg))
}

wci <- t(sapply(1:nrow(study_data), function(i)
  wilson_ci(study_data$Events[i], study_data$N[i])))
study_data$rate <- study_data$Events / study_data$N
study_data$lo   <- wci[,1]
study_data$hi   <- wci[,2]
study_data$wt   <- RE_WEIGHTS

# ── Verification ─────────────────────────────────────────────────────────────
cat("=== Step 7 v2 verified values ===\n")
cat(sprintf("  Pooled: %.2f%% [%.2f, %.2f%%]  PI: [%.2f, %.2f%%]\n",
            100*POOLED_EST, 100*POOLED_LO, 100*POOLED_HI, 100*PI_LO, 100*PI_HI))
cat(sprintf("  I2=%.1f%%, tau2=%.6f, Q=%.2f, Egger p=%.4f\n",
            I2_PCT, TAU2, Q_VAL, EGGER_P))
cat(sprintf("  Sensitivity (protocol): %.2f%% [%.2f, %.2f%%], k=%d, N=%d\n",
            100*SENS_EST, 100*SENS_LO, 100*SENS_HI, SENS_K, SENS_N))
cat(sprintf("  Weight range: %.1f%%–%.1f%%\n", min(RE_WEIGHTS), max(RE_WEIGHTS)))

# Cross-check per-study against Step 7 v2 output
cat("\n── Per-study verification ──\n")
step7_N <- c(6, 5, 176, 5, 96, 10, 29, 42, 31, 50, 88)
step7_R <- c(2, 0,   2, 5, 16,  1,  7,  2,  5,  4, 11)
step7_W <- c(6.6, 6.1, 11.1, 6.1, 10.9, 7.9, 9.8, 10.3, 9.9, 10.4, 10.8)
all_ok <- TRUE
for (i in 1:nrow(study_data)) {
  n_ok <- study_data$N[i] == step7_N[i]
  r_ok <- study_data$Events[i] == step7_R[i]
  w_ok <- study_data$wt[i] == step7_W[i]
  status <- ifelse(n_ok & r_ok & w_ok, "OK", "MISMATCH")
  if (status != "OK") all_ok <- FALSE
  cat(sprintf("  %-15s N=%3d R=%2d W=%.1f%% : %s\n",
              study_data$Study[i], study_data$N[i], study_data$Events[i],
              study_data$wt[i], status))
}
cat(sprintf("  Overall: %s\n", ifelse(all_ok, "ALL MATCHED", "ERRORS FOUND")))
cat("==============================\n\n")

###############################################################################
# BUILD TABLE — short Study labels only (no footnotes in table!)
###############################################################################

dt <- data.frame(
  Study      = study_data$Study,
  est        = study_data$rate,
  lo         = study_data$lo,
  hi         = study_data$hi,
  `Rate % [95% CI]` = sprintf("     %.1f [%.1f, %.1f]",
                               100*study_data$rate,
                               100*study_data$lo,
                               100*study_data$hi),
  `Weight %` = sprintf("%.1f", study_data$wt),
  is_summary = FALSE,
  check.names = FALSE, stringsAsFactors = FALSE
)

blank <- data.frame(
  Study="", est=NA, lo=NA, hi=NA,
  `Rate % [95% CI]`="", `Weight %`="",
  is_summary=FALSE, check.names=FALSE, stringsAsFactors=FALSE
)

pooled_r <- data.frame(
  Study = sprintf("RE Model (k=11, N=%d)", sum(study_data$N)),
  est=POOLED_EST, lo=POOLED_LO, hi=POOLED_HI,
  `Rate % [95% CI]` = sprintf("     %.1f [%.1f, %.1f]", 100*POOLED_EST, 100*POOLED_LO, 100*POOLED_HI),
  `Weight %`="",
  is_summary=TRUE, check.names=FALSE, stringsAsFactors=FALSE
)

sens_r <- data.frame(
  Study = sprintf("Sensitivity (k=%d, N=%d)", SENS_K, SENS_N),
  est=SENS_EST, lo=SENS_LO, hi=SENS_HI,
  `Rate % [95% CI]` = sprintf("     %.1f [%.1f, %.1f]", 100*SENS_EST, 100*SENS_LO, 100*SENS_HI),
  `Weight %`="",
  is_summary=TRUE, check.names=FALSE, stringsAsFactors=FALSE
)

pi_r <- data.frame(
  Study = sprintf("95%% PI: [%.1f%%, %.1f%%]", 100*PI_LO, 100*PI_HI),
  est=NA, lo=NA, hi=NA,
  `Rate % [95% CI]`="", `Weight %`="",
  is_summary=FALSE, check.names=FALSE, stringsAsFactors=FALSE
)

# NO footnote rows in table — they go in footnote= parameter
dt_all <- rbind(dt, blank, pooled_r, sens_r, pi_r)

# WIDE CI column — 40 spaces
dt_all$`  ` <- paste(rep(" ", 40), collapse = " ")

# Sizes
sizes_vec <- c(
  0.2 + 0.4 * sqrt(study_data$wt / max(study_data$wt)),  # 11 studies
  NA, 0.6, 0.5, NA  # blank, pooled, sens, PI
)

###############################################################################
# THEME
###############################################################################
tm <- forest_theme(
  base_size    = 11,
  ci_pch       = 15,
  ci_col       = "#2C5F8A",
  ci_lty       = 1,
  ci_lwd       = 2,
  ci_Theight   = 0.15,
  summary_pch  = 18,
  summary_col  = c("#B22222", "#4A7C59"),
  summary_lwd  = 2.5,
  refline_lty  = 2,
  refline_lwd  = 0.8,
  refline_col  = "#B22222",
  vertline_lty = 3,
  vertline_lwd = 0.7,
  vertline_col = "#D4A574",
  footnote_gp  = gpar(fontsize = 8, fontface = "italic", col = "gray40"),
  core = list(
    fg_params = list(hjust = 0, x = 0.04),
    bg_params = list(fill = c("#FFFFFF", "#F5F5F5")),
    padding   = unit(c(5, 3), "mm")
  ),
  colhead = list(
    fg_params = list(hjust = 0, x = 0.04, fontface = "bold")
  )
)

###############################################################################
# PLOT — no footnote here, added manually via grid.text at save time
###############################################################################
display_cols <- c("Study", "  ", "Rate % [95% CI]", "Weight %")

# Note: PI_LO=0 means the left PI vertical line sits at x=0 (the y-axis).
# We only show the right PI bound as a vertical line to avoid visual clutter.
p <- forest(
  data       = dt_all[, display_cols],
  est        = dt_all$est,
  lower      = dt_all$lo,
  upper      = dt_all$hi,
  sizes      = sizes_vec,
  is_summary = dt_all$is_summary,
  ci_column  = 2,
  ref_line   = POOLED_EST,
  vert_line  = c(PI_HI),
  xlim       = c(0, 1),
  ticks_at   = c(0, 0.25, 0.50, 0.75, 1.00),
  ticks_digits = 2,
  xlab       = "LRS Rescue Rate (proportion)",
  theme      = tm
)

# Footnote text (will be drawn via grid.text below the plot)
fn_line1 <- sprintf("Heterogeneity: I2 = %.1f%%, Q = %.2f (df = %d, p %s), tau2 = %.4f (arcsine scale)",
                    I2_PCT, Q_VAL, Q_DF, Q_P, TAU2)
fn_line2 <- sprintf("Egger's regression: p = %.3f  |  Small studies (N <= 6) in gray  |  Per-study N in Table 2",
                    EGGER_P)

# ── Post-editing ─────────────────────────────────────────────────────────────
n_s <- nrow(study_data)  # 11

p <- edit_plot(p, row = n_s + 2, gp = gpar(fontface="bold", col="#B22222"))  # pooled
p <- edit_plot(p, row = n_s + 3, gp = gpar(fontface="bold", col="#4A7C59"))  # sens
p <- edit_plot(p, row = n_s + 4, gp = gpar(fontface="italic", col="#D4A574", fontsize=10))  # PI

# Gray for small studies
for (r in which(study_data$N <= 6)) {
  p <- edit_plot(p, row = r, col = 1, gp = gpar(col = "#999999"))
}

# Borders
p <- add_border(p, part = "header", where = "bottom", gp = gpar(lwd = 1.5))
p <- add_border(p, row = n_s + 2, where = "top", gp = gpar(lwd = 1))

# Title
p <- insert_text(p,
  text = "Figure 2. Random-Effects Meta-Analysis of LRS Rescue Yield in SRS-Negative Patients",
  part = "header",
  gp = gpar(fontsize = 13, fontface = "bold")
)

###############################################################################
# SAVE
###############################################################################
f_pdf <- file.path(output_dir, "Figure2_ForestPlot.pdf")
f_png <- file.path(output_dir, "Figure2_ForestPlot.png")

cat("Saving PDF...\n")
pdf(f_pdf, width = 14, height = 10)
plot(p)
# Footnotes — centered under CI column, tight below x-axis label
grid.text(fn_line1, x = 0.22, y = 0.27,
          just = "left",
          gp = gpar(fontsize = 8, fontface = "italic", col = "gray40"))
grid.text(fn_line2, x = 0.22, y = 0.25,
          just = "left",
          gp = gpar(fontsize = 8, fontface = "italic", col = "gray40"))
# Legend for vertical lines
grid.segments(x0 = 0.22, x1 = 0.245, y0 = 0.23, y1 = 0.23,
              gp = gpar(col = "#B22222", lty = 2, lwd = 1.5))
grid.text("Pooled estimate (16.4%)", x = 0.25, y = 0.23,
          just = "left",
          gp = gpar(fontsize = 7.5, fontface = "italic", col = "gray40"))
grid.segments(x0 = 0.45, x1 = 0.475, y0 = 0.23, y1 = 0.23,
              gp = gpar(col = "#D4A574", lty = 3, lwd = 1.5))
grid.text("95% PI upper bound (66.8%); lower bound at 0.0%", x = 0.48, y = 0.23,
          just = "left",
          gp = gpar(fontsize = 7.5, fontface = "italic", col = "gray40"))
dev.off()
if (file.exists(f_pdf)) cat(sprintf("  OK: %s (%.0f KB)\n", f_pdf, file.size(f_pdf)/1024))

cat("Saving PNG...\n")
png(f_png, width = 14, height = 10, units = "in", res = 600)
plot(p)
grid.text(fn_line1, x = 0.22, y = 0.27,
          just = "left",
          gp = gpar(fontsize = 8, fontface = "italic", col = "gray40"))
grid.text(fn_line2, x = 0.22, y = 0.25,
          just = "left",
          gp = gpar(fontsize = 8, fontface = "italic", col = "gray40"))
grid.segments(x0 = 0.22, x1 = 0.245, y0 = 0.23, y1 = 0.23,
              gp = gpar(col = "#B22222", lty = 2, lwd = 1.5))
grid.text("Pooled estimate (16.4%)", x = 0.25, y = 0.23,
          just = "left",
          gp = gpar(fontsize = 7.5, fontface = "italic", col = "gray40"))
grid.segments(x0 = 0.45, x1 = 0.475, y0 = 0.23, y1 = 0.23,
              gp = gpar(col = "#D4A574", lty = 3, lwd = 1.5))
grid.text("95% PI upper bound (66.8%); lower bound at 0.0%", x = 0.48, y = 0.23,
          just = "left",
          gp = gpar(fontsize = 7.5, fontface = "italic", col = "gray40"))
dev.off()
if (file.exists(f_png)) cat(sprintf("  OK: %s (%.0f KB)\n", f_png, file.size(f_png)/1024))

cat("\n==========================================================\n")
cat("  DONE. All values from verified Step 7 v2 OUTPUT.\n")
cat(sprintf("  Pooled: %.2f%% [%.2f, %.2f%%]\n", 100*POOLED_EST, 100*POOLED_LO, 100*POOLED_HI))
cat(sprintf("  PI: [%.1f, %.2f%%]\n", 100*PI_LO, 100*PI_HI))
cat(sprintf("  Sensitivity: %.2f%% [%.2f, %.2f%%] (k=%d, N=%d)\n",
            100*SENS_EST, 100*SENS_LO, 100*SENS_HI, SENS_K, SENS_N))
cat(sprintf("  Output: %s\n", output_dir))
cat("==========================================================\n")
