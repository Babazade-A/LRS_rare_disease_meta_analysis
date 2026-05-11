###############################################################################
# STEP 6: VARIANT CHARACTERISTICS — R VERIFICATION SCRIPT
# Systematic Review: LRS Rescue in SRS-Negative Patients
# Input: Dataset_S1.xlsx, sheet "Tier 1 - Group A" (11 studies, Group A SRS-negative patients)
#
# Analysis restricted to 55 rescued cases (Definitive + Possible Rescue).
# All descriptive — no formal hypothesis testing.
#
# Key rules:
#   - Multi-value entries (e.g., "1; 33984") = separate variants → count each
#   - Variant sizes: each value counted individually (not max per patient)
#   - Genomic regions: each annotation counted individually
#   - ">" values: strip ">", use the number
#   - P520_Hiatt2021: excluded from size calculations (complex chromothripsis)
#   - "Expansion", "NA", text-only: excluded from size calculations
#   - ACMG: prioritize highest tier per patient (P > LP > VUS)
###############################################################################

# ── 0. LIBRARIES ─────────────────────────────────────────────────────────────
for (pkg in c("readxl", "dplyr", "tidyr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(readxl)
library(dplyr)
library(tidyr)

# ── 1a. SAVE OUTPUT TO TXT FILE ─────────────────────────────────────────────
output_file <- "Step6_Variant_Characteristics_OUTPUT.txt"
sink(output_file, split = TRUE)

# ── 1. DATA LOADING ──────────────────────────────────────────────────────────
dat <- read_excel(
  "Dataset_S1.xlsx",
  sheet     = "Tier 1 - Group A",
  skip      = 1,
  col_names = TRUE
)

cat("═══════════════════════════════════════════════════════════════════════\n")
cat("  STEP 6 — VARIANT CHARACTERISTICS: R VERIFICATION\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

N_total <- nrow(dat)

# ── 2. RESCUE CLASSIFICATION & FILTER ────────────────────────────────────────
dat <- dat %>%
  mutate(
    Rescued = as.integer(
      !is.na(LRS_Outcome) &
      (LRS_Outcome == "Definitive_Rescue" | grepl("Possible_Rescue", LRS_Outcome))
    )
  )

resc <- dat %>% filter(Rescued == 1)
N_resc <- nrow(resc)
cat(sprintf("Rescued cases: N = %d\n\n", N_resc))

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


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: VARIANT TYPE DISTRIBUTION
# ══════════════════════════════════════════════════════════════════════════════
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  SECTION 1: VARIANT TYPE DISTRIBUTION\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

vtype_raw <- as.character(resc$Rescued_Variant_Type)

# Detect involvement of each broad class per patient
has_sv    <- grepl("SV|MEI", vtype_raw, ignore.case = TRUE)
has_indel <- grepl("INDEL", vtype_raw, ignore.case = TRUE)
has_str   <- grepl("STR|Repeat|Expansion", vtype_raw, ignore.case = TRUE)
has_snv   <- grepl("\\bSNV\\b", vtype_raw, ignore.case = TRUE)
has_ch3   <- grepl("CH3|Methyl", vtype_raw, ignore.case = TRUE)

# Total involvement (any mention)
total_sv    <- sum(has_sv)
total_indel <- sum(has_indel)
total_str   <- sum(has_str)
total_snv   <- sum(has_snv)
total_ch3   <- sum(has_ch3)

# Pure (single class only)
n_classes <- as.integer(has_sv) + as.integer(has_indel) + as.integer(has_str) +
             as.integer(has_snv) + as.integer(has_ch3)
pure_sv    <- sum(has_sv & n_classes == 1)
pure_indel <- sum(has_indel & n_classes == 1)
pure_str   <- sum(has_str & n_classes == 1)
pure_snv   <- sum(has_snv & n_classes == 1)
pure_ch3   <- sum(has_ch3 & n_classes == 1)

is_compound <- n_classes > 1
n_compound  <- sum(is_compound)

cat("── Broad Variant Categories ──\n\n")
cat(sprintf("%-20s %5s %8s %5s %8s  [%s]\n",
            "Category", "Total", "Total%", "Pure", "Pure%", "Total 95% CI"))
cat("──────────────────────────────────────────────────────────────────────\n")

cat_names <- c("SV", "INDEL", "STR/Repeat", "SNV", "CH3/Methylation")
t_vals <- c(total_sv, total_indel, total_str, total_snv, total_ch3)
p_vals <- c(pure_sv, pure_indel, pure_str, pure_snv, pure_ch3)

for (i in seq_along(cat_names)) {
  ci <- wilson_ci(t_vals[i], N_resc)
  cat(sprintf("%-20s %5d %7.1f%% %5d %7.1f%%  [%5.1f%%, %5.1f%%]\n",
              cat_names[i], t_vals[i], 100 * t_vals[i] / N_resc,
              p_vals[i], 100 * p_vals[i] / N_resc,
              100 * ci["lower"], 100 * ci["upper"]))
}

cat(sprintf("\nPure (single class): %d | Compound (multi-class): %d | Total: %d\n",
            sum(p_vals), n_compound, N_resc))

sum_check_vtype <- sum(p_vals) + n_compound
cat(sprintf("Summation check: %d pure + %d compound = %d (expected %d) → %s\n\n",
            sum(p_vals), n_compound, sum_check_vtype, N_resc,
            ifelse(sum_check_vtype == N_resc, "PASS", "FAIL")))

# SV subtypes (patient-level)
cat("── SV Subtypes (among patients with SV involvement) ──\n\n")
sv_types_raw <- vtype_raw[has_sv]
sv_del     <- sum(grepl("DEL|Del", sv_types_raw) & !grepl("INDEL", sv_types_raw, ignore.case = TRUE))
sv_dup     <- sum(grepl("DUP|Dup", sv_types_raw))
sv_inv     <- sum(grepl("INV|Inv", sv_types_raw))
sv_ins     <- sum(grepl("SV_INS", sv_types_raw, ignore.case = TRUE))
sv_complex <- sum(grepl("Complex", sv_types_raw, ignore.case = TRUE))
sv_mei     <- sum(grepl("MEI", sv_types_raw, ignore.case = TRUE))

cat(sprintf("  SV Deletions:    %d\n", sv_del))
cat(sprintf("  SV Duplications: %d\n", sv_dup))
cat(sprintf("  SV Inversions:   %d\n", sv_inv))
cat(sprintf("  SV Insertions:   %d\n", sv_ins))
cat(sprintf("  SV Complex:      %d\n", sv_complex))
cat(sprintf("  SV MEI:          %d\n", sv_mei))
cat(sprintf("  Total SV patients: %d\n\n", total_sv))

# Compound variant combinations
cat("── Compound Variant Combinations ──\n\n")
compound_patients <- resc[is_compound, ]
cat(sprintf("Compound cases: %d\n", n_compound))
if (n_compound > 0) {
  comp_types <- sort(table(compound_patients$Rescued_Variant_Type), decreasing = TRUE)
  for (j in seq_along(comp_types)) {
    cat(sprintf("  '%s': %d\n", names(comp_types)[j], comp_types[j]))
  }
}
cat("\n")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2: VARIANT SIZE STATISTICS
# ══════════════════════════════════════════════════════════════════════════════
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  SECTION 2: VARIANT SIZE STATISTICS\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# Parse variant sizes — each semicolon-separated value is a SEPARATE variant
# UPDATED RULE (decisions.txt, Apr 2026):
#   - P520_Hiatt2021 chromothripsis: INCLUDED, manually parsed as two SV
#     components [126,000,000 bp pericentric inv; 9,300,000 bp CGR]
#   - "NA" / "Expansion" / non-numeric tokens: still excluded as text-only

parse_variant_sizes <- function(patient_ids, size_strings, study_ids) {
  # The chromothripsis special case requires BOTH Patient_ID_Original ==
  # "Proband_4" AND Study_ID == "Hiatt_2021" to fire, so a future "Proband_4"
  # from any other study cannot accidentally trigger the hand-parse.
  all_sizes <- numeric(0)
  all_pids  <- character(0)
  n_excluded_text  <- 0
  n_excluded_na    <- 0
  n_p520_manual    <- 0

  for (i in seq_along(size_strings)) {
    pid <- patient_ids[i]
    s   <- as.character(size_strings[i])
    study <- study_ids[i]

    # Hiatt_2021 / Proband_4 chromothripsis: hand-parse as [126e6, 9.3e6]
    # (Defensive: requires both pid AND study match.)
    if (!is.na(pid) && !is.na(study) &&
        pid == "Proband_4" && study == "Hiatt_2021") {
      all_sizes <- c(all_sizes, 126000000, 9300000)
      all_pids  <- c(all_pids, pid, pid)
      n_p520_manual <- n_p520_manual + 1
      next
    }

    if (is.na(s) || trimws(s) %in% c("NA", "N/A", "", "None")) {
      n_excluded_na <- n_excluded_na + 1
      next
    }

    # Split on semicolon and comma
    parts <- unlist(strsplit(s, "[;,]"))

    for (part in parts) {
      part <- trimws(part)
      # Skip text-only entries
      if (part %in% c("Expansion", "NA", "N/A", "")) {
        n_excluded_text <- n_excluded_text + 1
        next
      }
      # Strip ~, >, quotes
      part <- gsub("[\u0027\u2018\u2019\u0022\u201C\u201D\u0060~>]", "", part)
      part <- trimws(part)
      # Skip if still non-numeric
      val <- suppressWarnings(as.numeric(part))
      if (is.na(val)) {
        n_excluded_text <- n_excluded_text + 1
      } else {
        all_sizes <- c(all_sizes, val)
        all_pids  <- c(all_pids, pid)
      }
    }
  }

  list(sizes = all_sizes, pids = all_pids,
       n_excluded_na   = n_excluded_na,
       n_excluded_text = n_excluded_text,
       n_p520_manual   = n_p520_manual)
}

size_result <- parse_variant_sizes(resc$Patient_ID_Original, resc$Rescued_Variant_Size, resc$Study_ID)
all_sizes <- size_result$sizes
all_pids  <- size_result$pids

cat(sprintf("Variant size parsing:\n"))
cat(sprintf("  Patients: %d\n", N_resc))
cat(sprintf("  Individual variant sizes parsed: %d\n", length(all_sizes)))
cat(sprintf("  Excluded — NA/missing: %d patients\n", size_result$n_excluded_na))
cat(sprintf("  Excluded — text-only: %d entries\n", size_result$n_excluded_text))
cat(sprintf("  P520 chromothripsis hand-parsed: %d patient (2 SV components)\n", size_result$n_p520_manual))

# Overall descriptive stats
cat(sprintf("\n── Overall Variant Sizes (n=%d individual variants) ──\n\n", length(all_sizes)))

if (length(all_sizes) > 0) {
  cat(sprintf("  n:      %d\n", length(all_sizes)))
  cat(sprintf("  Median: %.1f bp\n", median(all_sizes)))
  cat(sprintf("  Mean:   %.1f bp\n", mean(all_sizes)))
  cat(sprintf("  IQR:    [%.1f, %.1f] bp\n", quantile(all_sizes, 0.25), quantile(all_sizes, 0.75)))
  cat(sprintf("  Range:  [%.1f, %.1f] bp\n", min(all_sizes), max(all_sizes)))
}

# Size by variant type — match each size back to its patient's variant type
cat("\n── Variant Sizes by Type ──\n\n")

# Create a mapping from patient ID to variant type flags
resc_lookup <- data.frame(
  pid     = resc$Patient_ID_Original,
  has_sv  = has_sv,
  has_indel = has_indel,
  has_str = has_str,
  has_snv = has_snv,
  has_ch3 = has_ch3,
  vtype   = vtype_raw,
  stringsAsFactors = FALSE
)

size_df <- data.frame(pid = all_pids, size_bp = all_sizes, stringsAsFactors = FALSE)
size_df <- merge(size_df, resc_lookup, by = "pid", all.x = TRUE)

# ──────────────────────────────────────────────────────────────────
# SECTION 2A: STRUCTURAL VARIANT SIZE DISTRIBUTION (HEADLINE STAT)
# Rule (decisions.txt, Apr 2026):
#   - Variant-level (each SV counted individually, not max-per-patient)
#   - Includes SV components from BOTH pure-SV and compound rescued cases
#   - P520_Hiatt2021 included as [126e6, 9.3e6]
#   - Filter ≥ 50 bp (canonical SV cutoff; drops SNV "1" tokens
#     that appear in compound-case size strings)
# ──────────────────────────────────────────────────────────────────
cat("\n── SV Size Distribution (per-variant, ≥50 bp, includes P520) ──\n\n")

sv_sizes <- size_df$size_bp[size_df$has_sv & size_df$size_bp >= 50]
if (length(sv_sizes) > 0) {
  cat(sprintf("  n SV variants:   %d\n", length(sv_sizes)))
  cat(sprintf("  n SV patients:   %d\n", length(unique(size_df$pid[size_df$has_sv & size_df$size_bp >= 50]))))
  cat(sprintf("  Median:          %.0f bp (%.1f kb)\n", median(sv_sizes), median(sv_sizes)/1000))
  cat(sprintf("  Mean:            %.0f bp\n", mean(sv_sizes)))
  cat(sprintf("  IQR:             [%.0f, %.0f] bp = [%.1f kb, %.1f kb]\n",
              quantile(sv_sizes, 0.25), quantile(sv_sizes, 0.75),
              quantile(sv_sizes, 0.25)/1000, quantile(sv_sizes, 0.75)/1000))
  cat(sprintf("  Range:           [%.0f, %.0f] bp = [%.0f bp, %.1f Mb]\n",
              min(sv_sizes), max(sv_sizes), min(sv_sizes), max(sv_sizes)/1e6))
  ord_min <- floor(log10(min(sv_sizes)))
  ord_max <- ceiling(log10(max(sv_sizes)))
  cat(sprintf("  Orders of mag:   %d (10^%d to 10^%d)\n", ord_max - ord_min, ord_min, ord_max))
}

# For pure SV patients: their sizes are SV sizes
# For compound patients: we can't definitively assign which size goes to which class
# Report overall and by patient type where unambiguous

report_size_stats <- function(sizes, label) {
  n <- length(sizes)
  if (n == 0) {
    cat(sprintf("  %-25s: n=0\n", label))
    return()
  }
  cat(sprintf("  %-25s: n=%d | Median=%.1f | IQR=[%.1f, %.1f] | Range=[%.1f, %.1f]\n",
              label, n, median(sizes),
              quantile(sizes, 0.25), quantile(sizes, 0.75),
              min(sizes), max(sizes)))
}

# Pure SV sizes
pure_sv_sizes <- size_df$size_bp[size_df$has_sv & !size_df$has_indel &
                                  !size_df$has_str & !size_df$has_snv & !size_df$has_ch3]
report_size_stats(pure_sv_sizes, "Pure SV")

# Pure INDEL sizes
pure_indel_sizes <- size_df$size_bp[size_df$has_indel & !size_df$has_sv &
                                     !size_df$has_str & !size_df$has_snv & !size_df$has_ch3]
report_size_stats(pure_indel_sizes, "Pure INDEL")

# Pure STR/Repeat sizes
pure_str_sizes <- size_df$size_bp[size_df$has_str & !size_df$has_sv &
                                   !size_df$has_indel & !size_df$has_snv & !size_df$has_ch3]
report_size_stats(pure_str_sizes, "Pure STR/Repeat")

# Pure SNV sizes
pure_snv_sizes <- size_df$size_bp[size_df$has_snv & !size_df$has_sv &
                                   !size_df$has_indel & !size_df$has_str & !size_df$has_ch3]
report_size_stats(pure_snv_sizes, "Pure SNV")

# Compound sizes (all compound patients' variants together)
compound_sizes <- size_df$size_bp[size_df$has_sv + size_df$has_indel + size_df$has_str +
                                   size_df$has_snv + size_df$has_ch3 > 1]
report_size_stats(compound_sizes, "Compound (all)")
cat("\n")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3: GENOMIC REGION
# ══════════════════════════════════════════════════════════════════════════════
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  SECTION 3: GENOMIC REGION (variant-level annotations)\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# Split on ; and , to get individual annotations
region_raw <- as.character(resc$Genomic_Region)

all_regions <- character(0)
for (i in seq_along(region_raw)) {
  s <- region_raw[i]
  if (is.na(s) || trimws(s) %in% c("NA", "N/A", "", "None")) next
  parts <- unlist(strsplit(s, "[;,]"))
  parts <- trimws(parts)
  parts <- parts[parts != "" & parts != "NA"]
  all_regions <- c(all_regions, parts)
}

n_annotations <- length(all_regions)
n_patients_with_region <- sum(!is.na(region_raw) & !region_raw %in% c("NA", "N/A", "", "None"))
n_extra <- n_annotations - n_patients_with_region

cat(sprintf("Patients with region data: %d / %d\n", n_patients_with_region, N_resc))
cat(sprintf("Total variant-level annotations: %d\n", n_annotations))
cat(sprintf("Extra annotations (multi-region patients): %d\n\n", n_extra))

# Normalize region names
normalize_region <- function(r) {
  r <- trimws(r)
  case_when(
    grepl("^Exonic$", r, ignore.case = TRUE)         ~ "Exonic",
    grepl("^3.?UTR$", r, ignore.case = TRUE)         ~ "UTR",
    grepl("^5.?UTR$", r, ignore.case = TRUE)         ~ "UTR",
    grepl("UTR.*Splice|Splice.*UTR", r, ignore.case = TRUE) ~ "UTR",
    grepl("^Deep.?Intronic$", r, ignore.case = TRUE) ~ "Deep_Intronic",
    grepl("^Intronic$", r, ignore.case = TRUE)       ~ "Intronic",
    grepl("^Non.?coding$", r, ignore.case = TRUE)    ~ "Non-coding",
    grepl("Multigene|Large.?Scale", r, ignore.case = TRUE) ~ "Multigene/Large-scale",
    grepl("^Genic$", r, ignore.case = TRUE)          ~ "Genic",
    grepl("Regulat|Promoter", r, ignore.case = TRUE) ~ "Regulatory/Promoter",
    grepl("^Complex$", r, ignore.case = TRUE)        ~ "Complex",
    TRUE ~ r
  )
}

all_regions_norm <- normalize_region(all_regions)
region_table <- sort(table(all_regions_norm), decreasing = TRUE)

cat(sprintf("%-25s %5s %8s  [%s]\n", "Region", "N", "%", "95% Wilson CI"))
cat("──────────────────────────────────────────────────────────────────\n")
for (j in seq_along(region_table)) {
  rname <- names(region_table)[j]
  rn    <- region_table[j]
  rpct  <- 100 * rn / n_annotations
  ci    <- wilson_ci(rn, n_annotations)
  cat(sprintf("%-25s %5d %7.1f%%  [%5.1f%%, %5.1f%%]\n",
              rname, rn, rpct, 100 * ci["lower"], 100 * ci["upper"]))
}
cat(sprintf("%-25s %5d\n", "TOTAL annotations", n_annotations))

# Non-exonic aggregation
non_exonic_regions <- c("UTR", "Deep_Intronic", "Intronic", "Non-coding",
                        "Multigene/Large-scale", "Genic", "Regulatory/Promoter", "Complex")
n_nonexonic <- sum(all_regions_norm %in% non_exonic_regions)
ci_nonex <- wilson_ci(n_nonexonic, n_annotations)
cat(sprintf("\nNon-exonic combined: %d / %d = %.1f%%  [%.1f%%, %.1f%%]\n\n",
            n_nonexonic, n_annotations, 100 * n_nonexonic / n_annotations,
            100 * ci_nonex["lower"], 100 * ci_nonex["upper"]))

# Summation verification
cat(sprintf("Region summation: %d annotations = %d base patients + %d extra → %s\n\n",
            n_annotations, n_patients_with_region, n_extra,
            ifelse(n_annotations == n_patients_with_region + n_extra, "PASS", "FAIL")))


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4: MOBILE ELEMENT & GENE INVOLVEMENT
# ══════════════════════════════════════════════════════════════════════════════
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  SECTION 4: MOBILE ELEMENT & GENE INVOLVEMENT\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# Mobile element involvement
me_raw <- as.character(resc$Mobile_element_involvement)
me_present <- !is.na(me_raw) & !me_raw %in% c("No", "None", "NA", "N/A", "")
n_me <- sum(me_present)
ci_me <- wilson_ci(n_me, N_resc)

cat(sprintf("Mobile element involvement: %d / %d = %.1f%%  [%.1f%%, %.1f%%]\n",
            n_me, N_resc, 100 * n_me / N_resc,
            100 * ci_me["lower"], 100 * ci_me["upper"]))
cat(sprintf("No involvement: %d\n\n", N_resc - n_me))

if (n_me > 0) {
  me_types <- table(me_raw[me_present])
  cat("Mobile element types:\n")
  for (j in seq_along(me_types)) {
    cat(sprintf("  %s: %d\n", names(me_types)[j], me_types[j]))
  }
  cat("\n")
}

# Gene involvement
cat("── Gene Involvement ──\n\n")
gene_raw <- as.character(resc$Gene)

# Split multi-gene entries on ; and ,
all_genes <- character(0)
for (g in gene_raw) {
  if (is.na(g) || trimws(g) %in% c("NA", "N/A", "", "None")) next
  parts <- unlist(strsplit(g, "[;,]"))
  parts <- trimws(parts)
  parts <- parts[parts != "" & parts != "NA"]
  all_genes <- c(all_genes, parts)
}

unique_genes <- unique(all_genes)
gene_counts  <- sort(table(all_genes), decreasing = TRUE)

cat(sprintf("Unique genes (before splitting multi-gene): %d\n",
            length(unique(gene_raw[!is.na(gene_raw) & gene_raw != "NA"]))))
cat(sprintf("Unique individual genes (after splitting): %d\n", length(unique_genes)))
cat(sprintf("Total gene entries: %d\n\n", length(all_genes)))

# Recurrent genes (n >= 2)
recurrent <- gene_counts[gene_counts >= 2]
singleton_genes <- gene_counts[gene_counts == 1]

cat(sprintf("Recurrent genes (n >= 2): %d\n", length(recurrent)))
for (j in seq_along(recurrent)) {
  cat(sprintf("  %s: %d patients\n", names(recurrent)[j], recurrent[j]))
}

cat(sprintf("\nSingleton genes (n = 1): %d\n\n", length(singleton_genes)))


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5: ZYGOSITY & INHERITANCE MODES
# ══════════════════════════════════════════════════════════════════════════════
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  SECTION 5: ZYGOSITY & INHERITANCE MODES\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# Zygosity
cat("── Zygosity ──\n\n")
zyg_raw <- as.character(resc$Zygosity)

# Normalize: "Heterozygous; Compound_Heterozygous" → assign to Compound_Heterozygous
# (prioritize compound over simple heterozygous)
zyg_norm <- case_when(
  grepl("Compound", zyg_raw, ignore.case = TRUE) ~ "Compound_Heterozygous",
  grepl("Heterozygous", zyg_raw, ignore.case = TRUE) ~ "Heterozygous",
  grepl("Hemizygous", zyg_raw, ignore.case = TRUE) ~ "Hemizygous",
  grepl("Homozygous", zyg_raw, ignore.case = TRUE) ~ "Homozygous",
  TRUE ~ zyg_raw
)

zyg_table <- sort(table(zyg_norm), decreasing = TRUE)

cat(sprintf("%-30s %5s %8s  [%s]\n", "Zygosity", "N", "%", "95% Wilson CI"))
cat("──────────────────────────────────────────────────────────────────\n")
for (j in seq_along(zyg_table)) {
  zname <- names(zyg_table)[j]
  zn    <- zyg_table[j]
  zpct  <- 100 * zn / N_resc
  ci    <- wilson_ci(zn, N_resc)
  cat(sprintf("%-30s %5d %7.1f%%  [%5.1f%%, %5.1f%%]\n",
              zname, zn, zpct, 100 * ci["lower"], 100 * ci["upper"]))
}
cat(sprintf("%-30s %5d\n", "TOTAL", sum(zyg_table)))

zyg_sum_ok <- sum(zyg_table) == N_resc
cat(sprintf("Summation check: %d (expected %d) → %s\n\n",
            sum(zyg_table), N_resc, ifelse(zyg_sum_ok, "PASS", "FAIL")))

# Compound heterozygous sub-categories
cat("── Compound Heterozygous Sub-categories ──\n\n")
comp_het_raw <- as.character(resc$Compound_Het_Structure)
comp_het_cases <- comp_het_raw[grepl("Compound", zyg_norm, ignore.case = TRUE)]
comp_het_table <- sort(table(comp_het_cases), decreasing = TRUE)
for (j in seq_along(comp_het_table)) {
  cat(sprintf("  %s: %d\n", names(comp_het_table)[j], comp_het_table[j]))
}
cat("\n")

# Inheritance mode
cat("── Inheritance Mode ──\n\n")
inh_raw <- as.character(resc$Inheritance_Mode)

# Normalize complex entries
inh_norm <- case_when(
  grepl("^AD$", inh_raw) ~ "Autosomal_Dominant",
  grepl("^AR$", inh_raw) ~ "Autosomal_Recessive",
  grepl("^XL[RD]?$", inh_raw) ~ "X-linked",
  grepl("XL", inh_raw, ignore.case = TRUE) ~ "X-linked",
  grepl("AD.*AR|AR.*AD", inh_raw) ~ "Mixed",
  grepl("Unknown|NA", inh_raw) ~ "Unknown",
  is.na(inh_raw) ~ "Unknown",
  TRUE ~ inh_raw
)

inh_table <- sort(table(inh_norm), decreasing = TRUE)

cat(sprintf("%-25s %5s %8s  [%s]\n", "Inheritance", "N", "%", "95% Wilson CI"))
cat("──────────────────────────────────────────────────────────────────\n")
for (j in seq_along(inh_table)) {
  iname <- names(inh_table)[j]
  in_   <- inh_table[j]
  ipct  <- 100 * in_ / N_resc
  ci    <- wilson_ci(in_, N_resc)
  cat(sprintf("%-25s %5d %7.1f%%  [%5.1f%%, %5.1f%%]\n",
              iname, in_, ipct, 100 * ci["lower"], 100 * ci["upper"]))
}
cat(sprintf("%-25s %5d\n", "TOTAL", sum(inh_table)))

inh_sum_ok <- sum(inh_table) == N_resc
cat(sprintf("Summation check: %d (expected %d) → %s\n\n",
            sum(inh_table), N_resc, ifelse(inh_sum_ok, "PASS", "FAIL")))

# X-linked subtypes
cat("── X-linked Subtypes ──\n\n")
xl_patients <- resc[grepl("XL", inh_raw, ignore.case = TRUE), ]
if (nrow(xl_patients) > 0) {
  xl_detail <- table(xl_patients$Inheritance_Mode)
  for (j in seq_along(xl_detail)) {
    cat(sprintf("  %s: %d\n", names(xl_detail)[j], xl_detail[j]))
  }
} else {
  cat("  No X-linked cases\n")
}
cat("\n")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6: ACMG CLASSIFICATION
# ══════════════════════════════════════════════════════════════════════════════
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  SECTION 6: ACMG CLASSIFICATION\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

acmg_raw <- as.character(resc$ACMG_Classification)

# Prioritize highest tier per patient: Pathogenic > Likely_Pathogenic > VUS
acmg_priority <- case_when(
  grepl("Pathogenic", acmg_raw) & !grepl("Likely", acmg_raw) ~ "Pathogenic",
  grepl("Likely_Pathogenic", acmg_raw) ~ ifelse(
    grepl("\\bPathogenic\\b", gsub("Likely_Pathogenic", "", acmg_raw)), "Pathogenic", "Likely_Pathogenic"
  ),
  grepl("VUS", acmg_raw) ~ "VUS",
  is.na(acmg_raw) ~ "Not_reported",
  TRUE ~ acmg_raw
)

# Simpler approach: check for "Pathogenic" after removing "Likely_Pathogenic"
acmg_clean <- sapply(acmg_raw, function(x) {
  if (is.na(x)) return("Not_reported")
  # Split on ; and ,
  parts <- unlist(strsplit(x, "[;,]"))
  parts <- trimws(parts)
  if (any(parts == "Pathogenic")) return("Pathogenic")
  if (any(grepl("Likely_Pathogenic", parts))) return("Likely_Pathogenic")
  if (any(grepl("VUS", parts))) return("VUS")
  return(x)
})

acmg_table <- sort(table(acmg_clean), decreasing = TRUE)

cat(sprintf("%-25s %5s %8s  [%s]\n", "ACMG Class", "N", "%", "95% Wilson CI"))
cat("──────────────────────────────────────────────────────────────────\n")
for (j in seq_along(acmg_table)) {
  aname <- names(acmg_table)[j]
  an    <- acmg_table[j]
  apct  <- 100 * an / N_resc
  ci    <- wilson_ci(an, N_resc)
  cat(sprintf("%-25s %5d %7.1f%%  [%5.1f%%, %5.1f%%]\n",
              aname, an, apct, 100 * ci["lower"], 100 * ci["upper"]))
}
cat(sprintf("%-25s %5d\n", "TOTAL", sum(acmg_table)))

acmg_sum_ok <- sum(acmg_table) == N_resc
cat(sprintf("Summation check: %d (expected %d) → %s\n\n",
            sum(acmg_table), N_resc, ifelse(acmg_sum_ok, "PASS", "FAIL")))

# Pathogenic + Likely_Pathogenic combined
n_p  <- ifelse("Pathogenic" %in% names(acmg_table), acmg_table["Pathogenic"], 0)
n_lp <- ifelse("Likely_Pathogenic" %in% names(acmg_table), acmg_table["Likely_Pathogenic"], 0)
n_plp <- n_p + n_lp
ci_plp <- wilson_ci(n_plp, N_resc)

cat(sprintf("Pathogenic + Likely_Pathogenic: %d / %d = %.1f%%  [%.1f%%, %.1f%%]\n\n",
            n_plp, N_resc, 100 * n_plp / N_resc,
            100 * ci_plp["lower"], 100 * ci_plp["upper"]))


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 7: INTERNAL CONSISTENCY DIAGNOSTICS
# ══════════════════════════════════════════════════════════════════════════════
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  SECTION 7: INTERNAL CONSISTENCY DIAGNOSTICS\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# Boundary verification
all_ci_vals <- numeric(0)
for (i in seq_along(t_vals)) {
  ci <- wilson_ci(t_vals[i], N_resc)
  all_ci_vals <- c(all_ci_vals, ci)
}
# Add region, zygosity, inheritance, ACMG CIs
for (tab in list(zyg_table, inh_table, acmg_table)) {
  for (j in seq_along(tab)) {
    ci <- wilson_ci(tab[j], N_resc)
    all_ci_vals <- c(all_ci_vals, ci)
  }
}
boundary_ok <- all(all_ci_vals >= 0, na.rm = TRUE) & all(all_ci_vals <= 1, na.rm = TRUE)

checks <- c(
  rescued_N        = N_resc == 55,
  vtype_sum        = sum_check_vtype == N_resc,
  zyg_sum          = zyg_sum_ok,
  inh_sum          = inh_sum_ok,
  acmg_sum         = acmg_sum_ok,
  region_sum       = (n_annotations == n_patients_with_region + n_extra),
  boundary_ok      = boundary_ok
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

cat("── Variant Type Proportions ──\n")
for (i in seq_along(cat_names)) {
  ci <- wilson_ci(t_vals[i], N_resc)
  cat(sprintf("  %-20s: %d/%d = %.10f  CI=[%.10f, %.10f]\n",
              cat_names[i], t_vals[i], N_resc,
              t_vals[i] / N_resc, ci["lower"], ci["upper"]))
}

cat("\n── Variant Size Summary (all parsed sizes) ──\n")
if (length(all_sizes) > 0) {
  cat(sprintf("  n=%d | Median=%.10f | Mean=%.10f | IQR=[%.10f, %.10f]\n",
              length(all_sizes), median(all_sizes), mean(all_sizes),
              quantile(all_sizes, 0.25), quantile(all_sizes, 0.75)))
}

cat("\n── SV Size Summary (per-variant, ≥50 bp, includes P520) ──\n")
sv_sizes_exact <- size_df$size_bp[size_df$has_sv & size_df$size_bp >= 50]
if (length(sv_sizes_exact) > 0) {
  cat(sprintf("  n=%d | Median=%.0f bp | Mean=%.1f bp | IQR=[%.0f, %.0f] bp | Range=[%.0f, %.0f] bp\n",
              length(sv_sizes_exact), median(sv_sizes_exact), mean(sv_sizes_exact),
              quantile(sv_sizes_exact, 0.25), quantile(sv_sizes_exact, 0.75),
              min(sv_sizes_exact), max(sv_sizes_exact)))
}

cat("\n── ACMG ──\n")
cat(sprintf("  P+LP: %d/%d = %.10f  CI=[%.10f, %.10f]\n",
            n_plp, N_resc, n_plp / N_resc, ci_plp["lower"], ci_plp["upper"]))

cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat("  STEP 6 VERIFICATION COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════════════\n")

# ── Close output file ────────────────────────────────────────────────────────
sink()
cat(sprintf("\nOutput saved to: %s\n", output_file))
