#!/usr/bin/env Rscript
# =============================================================================
# Figure 2b — HPO Phenotype Associations with LRS Rescue (UPDATED)
# Top 5 Neurological + Top 5 Non-Neurological categories (deduplicated)
# Matches style of forest_plot_primary.R
# Data source: patient_hpo_categories.xlsx (corrected pipeline)
#
# CORRECTED-PIPELINE NUMBERS (all verified against Step 3 OUTPUT):
#   156 HPO categories tested (not 81)
#   65 FDR-significant (not 48)
#   n = 538 patients, 55 rescued
#   51 within-neurology categories (272 neurological patients)
#   Seizure q = 0.052 → no longer FDR-significant
#
# Column names match the corrected Excel (underscores, not spaces):
#   HPO_ID, Category_Name, Depth, N_In, Rescued_In, Odds_Ratio,
#   CI_Lower, CI_Upper, Q_value
# =============================================================================

script_dir <- tryCatch(
  dirname(rstudioapi::getSourceEditorContext()$path),
  error = function(e) getwd()
)

library(grid)
library(forestploter)
library(readxl)
library(dplyr)

# =============================================================================
# 1. READ DATA FROM BOTH SHEETS
# =============================================================================
xlsx_path <- file.path(script_dir, "patient_hpo_categories.xlsx")

# Primary analysis — 156 categories (corrected from the buggy 81)
raw <- read_excel(xlsx_path, sheet = "Statistical Analysis")

dat <- data.frame(
  hpo_id   = raw$HPO_ID,
  category = raw$Category_Name,
  depth    = raw$Depth,
  n_in     = raw$N_In,
  rescued  = raw$Rescued_In,
  or       = raw$Odds_Ratio,
  ci_lo    = raw$CI_Lower,
  ci_hi    = raw$CI_Upper,
  qval     = raw$Q_value,
  stringsAsFactors = FALSE
)

dat <- dat[!is.na(dat$or), ]
cat("Read", nrow(dat), "categories from Statistical Analysis sheet\n")

# Sanity check against corrected pipeline (Step 3 OUTPUT)
stopifnot(nrow(dat) == 156)
stopifnot(sum(dat$qval < 0.05, na.rm = TRUE) == 65)

# Within-Neurology sheet — to identify which HPO IDs are neuro
neuro_raw <- read_excel(xlsx_path, sheet = "Within-Neurology Analysis", skip = 3)
neuro_ids <- neuro_raw$HPO_ID
neuro_ids <- neuro_ids[!is.na(neuro_ids)]
stopifnot(length(neuro_ids) == 51)  # 51 neuro sub-categories in corrected pipeline

# Tag each category as neuro or non-neuro
dat$is_neuro <- dat$hpo_id %in% neuro_ids

cat("Neuro categories:", sum(dat$is_neuro), "\n")
cat("Non-neuro categories:", sum(!dat$is_neuro), "\n")

# =============================================================================
# 2. DEDUPLICATE PARENT/CHILD TERMS WITH IDENTICAL STATS
#    HPO hierarchy causes e.g. "Spasticity" & "Upper motor neuron dysfunction"
#    to have identical (n_in, rescued, OR). Keep deepest (most specific) term.
# =============================================================================
# First dedup within same neuro/non-neuro group
dat <- dat %>%
  group_by(n_in, rescued, or, is_neuro) %>%
  arrange(desc(depth)) %>%
  slice(1) %>%
  ungroup()

# Also dedup across groups (e.g. Microcephaly=neuro vs Decreased head circ=non-neuro)
# Keep the deepest term overall when stats are identical
dat <- dat %>%
  group_by(n_in, rescued, or) %>%
  arrange(desc(depth)) %>%
  slice(1) %>%
  ungroup()

cat("After deduplication:", nrow(dat), "categories\n")

# =============================================================================
# 3. SELECT TOP 5 FROM EACH GROUP (BY OR, DESCENDING)
# =============================================================================
top_neuro <- dat %>%
  filter(is_neuro) %>%
  arrange(desc(or)) %>%
  head(5)

top_nonneuro <- dat %>%
  filter(!is_neuro) %>%
  arrange(desc(or)) %>%
  head(5)

cat("\nTop 5 Neuro picks:\n")
for (i in 1:nrow(top_neuro)) {
  cat(sprintf("  %s: OR=%.2f, n=%d, rescued=%d, q=%.4f\n",
              top_neuro$category[i], top_neuro$or[i],
              top_neuro$n_in[i], top_neuro$rescued[i], top_neuro$qval[i]))
}

cat("\nTop 5 Non-Neuro picks:\n")
for (i in 1:nrow(top_nonneuro)) {
  cat(sprintf("  %s: OR=%.2f, n=%d, rescued=%d, q=%.4f\n",
              top_nonneuro$category[i], top_nonneuro$or[i],
              top_nonneuro$n_in[i], top_nonneuro$rescued[i], top_nonneuro$qval[i]))
}

# Combine: neuro on top, then non-neuro, each sorted by OR desc
sel <- bind_rows(
  top_neuro %>% mutate(group = "Neurological"),
  top_nonneuro %>% mutate(group = "Non-Neurological")
)

cat("\nTotal categories for plot:", nrow(sel), "\n")

# =============================================================================
# 4. BUILD DISPLAY COLUMNS
# =============================================================================
sel$sig <- sel$qval < 0.05

# Format category name with HPO code
sel$Category <- sapply(1:nrow(sel), function(i) {
  nm <- sel$category[i]
  hpo <- sel$hpo_id[i]
  # Only abbreviate CNS to keep readable
  nm <- gsub("central nervous system", "CNS", nm, ignore.case = TRUE)
  paste0(nm, " (", hpo, ")")
})
# Pad all Category strings to fixed width so forestploter allocates enough space
max_len <- max(nchar(sel$Category)) + 8
sel$Category <- formatC(sel$Category, width = -max_len, flag = "-")

sel$`Rescue\n(n/N)` <- paste0(sel$rescued, "/", sel$n_in)

sel$` ` <- paste(rep(" ", 30), collapse = " ")

sel$`OR (95% CI)` <- sprintf("%.2f (%.2f\u2013%.2f)", sel$or, sel$ci_lo, sel$ci_hi)

sel$`q-value` <- ifelse(sel$qval < 0.001,
                        formatC(sel$qval, format = "e", digits = 1),
                        sprintf("%.3f", sel$qval))

# =============================================================================
# 5. FORESTPLOTER
# =============================================================================
disp_cols <- c("Category", "Rescue\n(n/N)", " ", "OR (95% CI)", "q-value")
plot_data <- sel[, disp_cols]

tm <- forest_theme(
  base_size     = 11,
  ci_pch        = 15,
  ci_col        = "black",
  ci_fill       = "black",
  ci_lwd        = 1.5,
  ci_Theight    = unit(0.18, "inches"),
  refline_gp    = gpar(lwd = 1, lty = "dashed", col = "grey40"),
  xaxis_gp      = gpar(fontsize = 10, lwd = 0.5),
  footnote_gp   = gpar(fontsize = 8, fontface = "italic", col = "grey30"),
  title_gp      = gpar(fontsize = 14, fontface = "bold"),
  core          = list(
    bg_params = list(fill = c("white", "#f5f5f5"), col = c("white", "#f5f5f5")),
    padding   = unit(c(4, 3.5), "mm")
  ),
  colhead       = list(
    bg_params = list(fill = "white", col = "white"),
    fg_params = list(fontface = 2L, fontsize = 11, hjust = 0, x = 0.05)
  )
)

p <- forest(
  data      = plot_data,
  est       = sel$or,
  lower     = sel$ci_lo,
  upper     = sel$ci_hi,
  sizes     = ifelse(sel$n_in > 100, 0.7, ifelse(sel$n_in > 30, 0.5, 0.35)),
  ci_column = 3,
  ref_line  = 1,
  x_trans   = "log",
  xlim      = c(0.5, 100),
  ticks_at  = c(0.5, 1, 2, 5, 10, 20, 50, 100),
  xlab      = "Odds Ratio (log scale)",
  arrow_lab = NULL,
  title     = "Figure 2b: HPO Phenotype Associations with Diagnostic Rescue",
  footnote  = paste(
    "",
    "Red = FDR-significant (q < 0.05).",
    "n = 538 patients, 55 rescued. 156 HPO categories tested; 65 FDR-significant.",
    sep = "\n"
  ),
  theme     = tm
)

# =============================================================================
# 6. COLOR BY SIGNIFICANCE
# =============================================================================
sig_rows  <- which(sel$sig)
nsig_rows <- which(!sel$sig)

if (length(sig_rows) > 0) {
  p <- edit_plot(p, row = sig_rows, col = 3, which = "ci",
                 gp = gpar(col = "#C0392B", fill = "#C0392B"))
  p <- edit_plot(p, row = sig_rows, col = 1,
                 gp = gpar(col = "#C0392B", fontface = "bold"))
}

if (length(nsig_rows) > 0) {
  p <- edit_plot(p, row = nsig_rows, col = 3, which = "ci",
                 gp = gpar(col = "#95A5A6", fill = "#95A5A6"))
  p <- edit_plot(p, row = nsig_rows, col = 1,
                 gp = gpar(col = "#7F8C8D"))
}

# Add group separator: border between neuro and non-neuro
n_neuro_rows <- nrow(top_neuro)
p <- add_border(p, part = "body", row = n_neuro_rows, where = "bottom",
                gp = gpar(lwd = 0.8, lty = "dotted", col = "grey50"))

# Add group labels via insert_text
p <- insert_text(p,
                 text = "Abnormality of the nervous system",
                 row = 1,
                 just = "left",
                 gp = gpar(fontsize = 10, fontface = "bold.italic",
                           col = "#2C3E50"))

p <- insert_text(p,
                 text = "Other Phenotypes",
                 row = n_neuro_rows + 2,
                 just = "left",
                 gp = gpar(fontsize = 10, fontface = "bold.italic",
                           col = "#2C3E50"))

p <- add_border(p, part = "header", row = 1, where = "bottom",
                gp = gpar(lwd = 1))

# =============================================================================
# 7. EXPORT
# =============================================================================
out_path <- file.path(script_dir, "figure_2b_NEWONE.pdf")

pdf(out_path, width = 18, height = 8)
plot(p)
dev.off()

cat("\n\u2713 Saved:", out_path, "\n")
cat("  Dimensions: 18 x 8 in (vector PDF)\n")
cat("  Categories plotted:", nrow(sel), "\n")
