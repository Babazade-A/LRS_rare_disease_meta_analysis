###############################################################################
# generate_corrected_excel.R
#
# Regenerates patient_hpo_categories.xlsx from source data (Dataset_S1.xlsx
# Tier 1 - Group A sheet, and hp.json) with corrected DAG traversal and
# verified statistics.
#
# Category selection: bottom-up algorithm (all HPO ancestors including self),
# minimum 5 patients per category, depth >= 2 (at or below organ-system level).
# Categories are determined algorithmically — not hardcoded.
#
# Produces 3 analysis sheets:
#   1. Statistical Analysis       (Primary, all patients)
#   2. Sensitivity Analysis       (Characterized only, excluding patients w/o HPO)
#   3. Within-Neurology Analysis  (Neuro patients only)
#
# Plus 4 reference sheets: Rescue Rates, Category Summary, Patient Categories,
# HPO Term Lookup
#
# Requirements: readxl, jsonlite, openxlsx
###############################################################################

cat("Loading libraries...\n")
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")

library(readxl)
library(jsonlite)
library(openxlsx)

# Open sink for output capture (split=TRUE prints to console AND file)
output_file <- "Step3_HPO_OUTPUT.txt"
sink(output_file, split = TRUE)

###############################################################################
# S1. LOAD HPO ONTOLOGY AND BUILD DAG
###############################################################################
cat("S1. Loading HPO ontology...\n")

hpo_raw <- fromJSON("hp.json", simplifyVector = FALSE)
graph   <- hpo_raw$graphs[[1]]

# Build node lookup (HPO ID -> label)
node_labels <- list()
for (node in graph$nodes) {
  nid <- node$id
  if (grepl("HP_", nid, fixed = TRUE)) {
    hp_id <- gsub("http://purl.obolibrary.org/obo/HP_", "HP:", nid, fixed = TRUE)
    node_labels[[hp_id]] <- ifelse(is.null(node$lbl), "", node$lbl)
  }
}
cat(sprintf("  Loaded %d HPO nodes\n", length(node_labels)))

# Build parent/child adjacency lists
parents_of <- list()   # parents_of[["HP:xxx"]] = c("HP:yyy", ...)
children_of <- list()  # children_of[["HP:xxx"]] = c("HP:yyy", ...)

for (edge in graph$edges) {
  if (!is.null(edge$pred) && edge$pred == "is_a") {
    sub_id <- gsub("http://purl.obolibrary.org/obo/HP_", "HP:", edge$sub, fixed = TRUE)
    obj_id <- gsub("http://purl.obolibrary.org/obo/HP_", "HP:", edge$obj, fixed = TRUE)
    if (startsWith(sub_id, "HP:") && startsWith(obj_id, "HP:")) {
      parents_of[[sub_id]]  <- c(parents_of[[sub_id]], obj_id)
      children_of[[obj_id]] <- c(children_of[[obj_id]], sub_id)
    }
  }
}

# Deduplicate
parents_of  <- lapply(parents_of, unique)
children_of <- lapply(children_of, unique)

edge_count <- sum(vapply(parents_of, length, integer(1)))
cat(sprintf("  Built DAG: %d is_a edges\n", edge_count))

# --- Helper: get ALL ancestors of a term INCLUDING ITSELF (BFS) ---
# This matches the original algorithm: "Return set of all ancestors including self"
get_all_ancestors <- function(hp_id) {
  visited <- character(0)
  stack   <- hp_id
  while (length(stack) > 0) {
    current <- stack[length(stack)]
    stack   <- stack[-length(stack)]
    if (current %in% visited) next
    visited <- c(visited, current)
    pars <- parents_of[[current]]
    if (!is.null(pars)) {
      stack <- c(stack, pars)
    }
  }
  unique(visited)
}

# --- Helper: get ALL descendants of a term INCLUDING ITSELF (BFS) ---
get_all_descendants <- function(hp_id) {
  visited <- character(0)
  queue   <- hp_id
  head    <- 1L
  while (head <= length(queue)) {
    current <- queue[head]
    head    <- head + 1L
    if (current %in% visited) next
    visited <- c(visited, current)
    kids <- children_of[[current]]
    if (!is.null(kids)) {
      queue <- c(queue, kids)
    }
  }
  unique(visited)
}

# Compute depths via BFS from root HP:0000001
cat("  Computing node depths...\n")
depth_of <- list()
depth_of[["HP:0000001"]] <- 0L
queue <- "HP:0000001"
head  <- 1L
while (head <= length(queue)) {
  current <- queue[head]
  head    <- head + 1L
  kids    <- children_of[[current]]
  if (!is.null(kids)) {
    for (kid in kids) {
      if (is.null(depth_of[[kid]])) {
        depth_of[[kid]] <- depth_of[[current]] + 1L
        queue <- c(queue, kid)
      }
    }
  }
}
cat(sprintf("  Depths computed for %d nodes\n", length(depth_of)))

###############################################################################
# S2. LOAD PATIENT DATA
###############################################################################
cat("S2. Loading patient data...\n")

# Read source data: row 1 of xlsx is a banner, row 2 contains column names,
# data starts row 3. skip=1 makes row 2 the header row of the resulting df.
df <- read_excel("Dataset_S1.xlsx", sheet = "Tier 1 - Group A",
                 skip = 1, col_names = FALSE)
# skip=1 skips the top header row; row 2 (sub-headers) becomes row 1 of df
# Actual data starts from row 2 of df

# Column indices (1-based) in the Tier 1 - Group A sheet:
#   Col 11 = HPO_Terms_given
#   Col 12 = HPO_Terms_inferred
#   Col 28 = LRS_Outcome
col_K           <- 11   # HPO_Terms_given
col_L           <- 12   # HPO_Terms_inferred
col_LRS_Outcome <- 28   # LRS_Outcome

# Parse HPO terms from a cell value
parse_hpo <- function(val) {
  if (is.na(val) || is.null(val)) return(character(0))
  val <- as.character(val)
  if (val %in% c("Not_reported", "NA", "None", "", "N/A")) return(character(0))
  m <- regmatches(val, gregexpr("HP:\\d{7}", val))[[1]]
  return(m)
}

n_patients <- nrow(df) - 1  # subtract sub-header row
cat(sprintf("  Raw rows: %d, data rows: %d\n", nrow(df), n_patients))

# Parse all patients (skip sub-header = row 1 of df)
patient_terms   <- vector("list", n_patients)
patient_rescued <- logical(n_patients)
patient_has_hpo <- logical(n_patients)

for (i in seq_len(n_patients)) {
  row_idx <- i + 1  # skip sub-header
  terms_k <- parse_hpo(df[[col_K]][row_idx])
  terms_l <- parse_hpo(df[[col_L]][row_idx])
  all_terms <- unique(c(terms_k, terms_l))

  # Keep only terms that exist in ontology (have a computed depth)
  valid_terms <- all_terms[vapply(all_terms, function(t) !is.null(depth_of[[t]]), logical(1))]

  patient_terms[[i]]   <- valid_terms
  patient_has_hpo[i]   <- length(valid_terms) > 0

  outcome <- as.character(df[[col_LRS_Outcome]][row_idx])
  patient_rescued[i]   <- outcome %in% c("Definitive_Rescue", "Possible_Rescue")
}

n_total           <- n_patients
n_rescued         <- sum(patient_rescued)
n_characterized   <- sum(patient_has_hpo)
n_uncharacterized <- n_total - n_characterized

cat(sprintf("  Total: %d, Rescued: %d\n", n_total, n_rescued))
cat(sprintf("  Characterized: %d, Uncharacterized: %d\n", n_characterized, n_uncharacterized))

###############################################################################
# S3. ASSIGN PATIENTS TO CATEGORIES (ALL ANCESTORS INCLUDING SELF)
###############################################################################
cat("S3. Assigning patients to HPO categories...\n")

# For each patient, compute the set of ALL categories they belong to
# (their direct terms + all ancestors of each term, including the term itself)
patient_cats <- vector("list", n_patients)
for (i in seq_len(n_patients)) {
  if (!patient_has_hpo[i]) {
    patient_cats[[i]] <- character(0)
    next
  }
  cats <- character(0)
  for (term in patient_terms[[i]]) {
    ancs <- get_all_ancestors(term)  # includes term itself
    cats <- c(cats, ancs)
  }
  patient_cats[[i]] <- unique(cats)
}
cat("  Done.\n")

###############################################################################
# S4. ALGORITHMIC CATEGORY SELECTION
###############################################################################
cat("S4. Selecting categories (>= 5 patients, depth >= 2)...\n")

MIN_CATEGORY_SIZE <- 5L
DEPTH_CAP         <- 2L    # shallowest permitted depth (organ-system level)

# Count patients per category
cat_patient_counts <- list()
for (i in seq_len(n_patients)) {
  for (hpo_cat in patient_cats[[i]]) {
    d <- depth_of[[hpo_cat]]
    if (!is.null(d) && d >= DEPTH_CAP) {
      if (is.null(cat_patient_counts[[hpo_cat]])) {
        cat_patient_counts[[hpo_cat]] <- 0L
      }
      cat_patient_counts[[hpo_cat]] <- cat_patient_counts[[hpo_cat]] + 1L
    }
  }
}

# Select qualifying categories
all_cat_ids <- names(cat_patient_counts)
qualifying_mask <- vapply(all_cat_ids, function(hpo_cat) {
  cat_patient_counts[[hpo_cat]] >= MIN_CATEGORY_SIZE
}, logical(1))

primary_categories <- all_cat_ids[qualifying_mask]
cat(sprintf("  Qualifying categories: %d\n", length(primary_categories)))

###############################################################################
# S5. IDENTIFY NEUROLOGICAL PATIENTS AND SELECT NEURO CATEGORIES
###############################################################################
cat("S5. Identifying neurological patients...\n")

neuro_root <- "HP:0000707"
# Get all descendants INCLUDING the root itself
neuro_desc_set <- get_all_descendants(neuro_root)

# A patient is neurological if ANY of their terms is the neuro root
# or a descendant of the neuro root.
# We check against neuro_desc_set which includes the root.
# Note: get_all_ancestors() would miss patients whose ONLY term is
# HP:0000707 itself.
is_neuro <- logical(n_patients)
for (i in seq_len(n_patients)) {
  if (!patient_has_hpo[i]) next
  for (term in patient_terms[[i]]) {
    if (term %in% neuro_desc_set) {
      is_neuro[i] <- TRUE
      break
    }
  }
}

neuro_indices  <- which(is_neuro)
n_neuro        <- length(neuro_indices)
n_neuro_rescue <- sum(patient_rescued[neuro_indices])
cat(sprintf("  Neurological patients: %d, Rescued: %d\n", n_neuro, n_neuro_rescue))
cat(sprintf("  Baseline rescue rate: %d/%d = %.1f%%\n",
            n_neuro_rescue, n_neuro, n_neuro_rescue / n_neuro * 100))

# Select neuro sub-categories: must be under nervous system AND >= 5 neuro patients
cat("  Selecting within-neurology categories...\n")
neuro_categories <- character(0)
for (hpo_cat in primary_categories) {
  # Category must be a descendant of (or equal to) neuro root
  if (!(hpo_cat %in% neuro_desc_set)) next
  # Count neuro patients in this category
  n_in_neuro <- sum(vapply(neuro_indices, function(i) hpo_cat %in% patient_cats[[i]], logical(1)))
  if (n_in_neuro >= MIN_CATEGORY_SIZE) {
    neuro_categories <- c(neuro_categories, hpo_cat)
  }
}
cat(sprintf("  Within-neurology categories: %d\n", length(neuro_categories)))

###############################################################################
# S6. STATISTICAL TESTING FUNCTION
###############################################################################

# Haldane-Anscombe correction for zero cells
haldane_or <- function(a, b, c, d) {
  if (a * d == 0 || b * c == 0) {
    a <- a + 0.5; b <- b + 0.5; c <- c + 0.5; d <- d + 0.5
  }
  (a * d) / (b * c)
}

# 95% CI using log-OR with Haldane correction
log_or_ci <- function(a, b, c, d) {
  if (a * d == 0 || b * c == 0) {
    a <- a + 0.5; b <- b + 0.5; c <- c + 0.5; d <- d + 0.5
  }
  log_or <- log((a * d) / (b * c))
  se     <- sqrt(1/a + 1/b + 1/c + 1/d)
  c(exp(log_or - 1.96 * se), exp(log_or + 1.96 * se))
}

run_fisher_analysis <- function(patient_idx, categories, label) {
  # patient_idx: integer vector of patient indices to include
  # categories: character vector of HPO IDs to test
  # Returns a data frame with results

  n_total_sub <- length(patient_idx)
  results <- data.frame(
    HPO_ID         = character(0),
    Category_Name  = character(0),
    Depth          = integer(0),
    N_In           = integer(0),
    Rescued_In     = integer(0),
    Not_Rescued_In = integer(0),
    Rescue_Rate_In = numeric(0),
    N_Out          = integer(0),
    Rescued_Out    = integer(0),
    Not_Rescued_Out= integer(0),
    Rescue_Rate_Out= numeric(0),
    Odds_Ratio     = numeric(0),
    CI_Lower       = numeric(0),
    CI_Upper       = numeric(0),
    P_value        = numeric(0),
    stringsAsFactors = FALSE
  )

  cat(sprintf("  Testing %d categories in %s (N=%d)...\n", length(categories), label, n_total_sub))

  for (hpo_cat in categories) {
    in_idx  <- patient_idx[vapply(patient_idx, function(i) hpo_cat %in% patient_cats[[i]], logical(1))]
    out_idx <- setdiff(patient_idx, in_idx)

    n_in  <- length(in_idx)
    n_out <- length(out_idx)
    r_in  <- sum(patient_rescued[in_idx])
    r_out <- sum(patient_rescued[out_idx])
    nr_in <- n_in - r_in
    nr_out<- n_out - r_out

    # Fisher's exact test (two-sided)
    ft    <- fisher.test(matrix(c(r_in, nr_in, r_out, nr_out), nrow = 2), alternative = "two.sided")
    p_val <- ft$p.value

    # Odds ratio with Haldane correction
    or_val <- haldane_or(r_in, nr_in, r_out, nr_out)
    ci     <- log_or_ci(r_in, nr_in, r_out, nr_out)

    rr_in  <- ifelse(n_in > 0, r_in / n_in, 0)
    rr_out <- ifelse(n_out > 0, r_out / n_out, 0)

    dep <- depth_of[[hpo_cat]]
    if (is.null(dep)) dep <- NA_integer_

    results <- rbind(results, data.frame(
      HPO_ID         = hpo_cat,
      Category_Name  = ifelse(is.null(node_labels[[hpo_cat]]), hpo_cat, node_labels[[hpo_cat]]),
      Depth          = as.integer(dep),
      N_In           = n_in,
      Rescued_In     = r_in,
      Not_Rescued_In = nr_in,
      Rescue_Rate_In = rr_in,
      N_Out          = n_out,
      Rescued_Out    = r_out,
      Not_Rescued_Out= nr_out,
      Rescue_Rate_Out= rr_out,
      Odds_Ratio     = or_val,
      CI_Lower       = ci[1],
      CI_Upper       = ci[2],
      P_value        = p_val,
      stringsAsFactors = FALSE
    ))
  }

  # Benjamini-Hochberg FDR correction
  m <- nrow(results)
  ord <- order(results$P_value)
  results$Q_value <- NA_real_
  results$Q_value[ord] <- pmin(results$P_value[ord] * m / seq_len(m), 1.0)
  # Enforce monotonicity (descending through sorted order)
  for (j in (m - 1):1) {
    idx_j   <- ord[j]
    idx_j1  <- ord[j + 1]
    results$Q_value[idx_j] <- min(results$Q_value[idx_j], results$Q_value[idx_j1])
  }

  results$Significance <- ifelse(results$Q_value < 0.05,
                                 "Significant (FDR < 0.05)",
                                 "Not significant")

  # Sort by OR descending
  results <- results[order(-results$Odds_Ratio), ]
  rownames(results) <- NULL
  return(results)
}

###############################################################################
# S7. RUN ALL THREE ANALYSES
###############################################################################

# --- Primary Analysis (all patients) ---
cat(sprintf("S7a. Running Primary Analysis (N=%d)...\n", n_total))
primary_results <- run_fisher_analysis(seq_len(n_total), primary_categories, "Primary")
n_sig_primary   <- sum(primary_results$Q_value < 0.05)
cat(sprintf("  %d categories tested, %d FDR-significant\n", nrow(primary_results), n_sig_primary))

# --- Sensitivity Analysis (characterized only) ---
cat("S7b. Running Sensitivity Analysis (characterized only)...\n")
char_indices    <- which(patient_has_hpo)
char_results    <- run_fisher_analysis(char_indices, primary_categories, "Sensitivity")
n_sig_char      <- sum(char_results$Q_value < 0.05)
n_char_rescued  <- sum(patient_rescued[char_indices])
cat(sprintf("  %d characterized patients, %d FDR-significant\n", length(char_indices), n_sig_char))

# Build sensitivity comparison table
sens_df <- data.frame(
  HPO_ID        = primary_results$HPO_ID,
  Category_Name = primary_results$Category_Name,
  Depth         = primary_results$Depth,
  N_Patients    = primary_results$N_In,
  N_Rescued     = primary_results$Rescued_In,
  Rescue_Rate   = primary_results$Rescue_Rate_In,
  stringsAsFactors = FALSE
)

# Match char_results to primary order
char_match <- match(sens_df$HPO_ID, char_results$HPO_ID)

sens_df$OR_Full   <- primary_results$Odds_Ratio
sens_df$CI_Full   <- paste0("[", round(primary_results$CI_Lower, 2), ", ",
                             round(primary_results$CI_Upper, 2), "]")
sens_df$Q_Full    <- primary_results$Q_value
sens_df$Sig_Full  <- ifelse(primary_results$Q_value < 0.05, "FDR<0.05", "NS")

sens_df$OR_Char   <- char_results$Odds_Ratio[char_match]
sens_df$CI_Char   <- paste0("[", round(char_results$CI_Lower[char_match], 2), ", ",
                             round(char_results$CI_Upper[char_match], 2), "]")
sens_df$Q_Char    <- char_results$Q_value[char_match]
sens_df$Sig_Char  <- ifelse(char_results$Q_value[char_match] < 0.05, "FDR<0.05", "NS")

sens_df$OR_Change_Pct <- (sens_df$OR_Char - sens_df$OR_Full) / sens_df$OR_Full
sens_df$Tier_Changed  <- ifelse(sens_df$Sig_Full != sens_df$Sig_Char, "Yes", NA)

# Sort by OR_Full descending
sens_df <- sens_df[order(-sens_df$OR_Full), ]
rownames(sens_df) <- NULL

# --- Within-Neurology Analysis ---
cat(sprintf("S7c. Running Within-Neurology Analysis (N=%d)...\n", n_neuro))
neuro_results <- run_fisher_analysis(neuro_indices, neuro_categories, "Within-Neurology")
n_sig_neuro   <- sum(neuro_results$Q_value < 0.05)
cat(sprintf("  %d categories tested, %d FDR-significant\n", nrow(neuro_results), n_sig_neuro))

###############################################################################
# S8. BUILD REFERENCE SHEETS
###############################################################################
cat("S8. Building reference sheets...\n")

# --- Rescue Rates (All) ---
rescue_df <- data.frame(
  HPO_ID        = primary_results$HPO_ID,
  Category_Name = primary_results$Category_Name,
  N_Total       = primary_results$N_In,
  N_Rescued     = primary_results$Rescued_In,
  Rescue_Rate   = round(primary_results$Rescue_Rate_In * 100, 1),
  stringsAsFactors = FALSE
)
rescue_df <- rescue_df[order(-rescue_df$Rescue_Rate), ]
rownames(rescue_df) <- NULL

# --- Category Summary ---
cat_summary <- data.frame(
  HPO_ID        = primary_results$HPO_ID,
  Category_Name = primary_results$Category_Name,
  Depth         = primary_results$Depth,
  N_Patients    = primary_results$N_In,
  stringsAsFactors = FALSE
)
cat_summary <- cat_summary[order(-cat_summary$N_Patients), ]
rownames(cat_summary) <- NULL

# --- Patient Categories ---
pat_cat_rows <- list()
for (i in seq_len(n_patients)) {
  if (!patient_has_hpo[i]) next
  cats_in_primary <- intersect(patient_cats[[i]], primary_categories)
  cat_names <- vapply(cats_in_primary, function(hp) {
    lbl <- node_labels[[hp]]
    if (is.null(lbl)) hp else lbl
  }, character(1))
  pat_cat_rows[[length(pat_cat_rows) + 1]] <- data.frame(
    Patient_Index = i,
    Rescued       = patient_rescued[i],
    N_Terms       = length(patient_terms[[i]]),
    N_Categories  = length(cats_in_primary),
    Categories    = paste(cat_names, collapse = "; "),
    stringsAsFactors = FALSE
  )
}
pat_cat_df <- do.call(rbind, pat_cat_rows)

# --- HPO Term Lookup ---
all_used_terms <- unique(unlist(patient_terms))
lookup_rows <- list()
for (term in sort(all_used_terms)) {
  lbl <- node_labels[[term]]
  dep <- depth_of[[term]]
  lookup_rows[[length(lookup_rows) + 1]] <- data.frame(
    HPO_ID = term,
    Label  = ifelse(is.null(lbl), "", lbl),
    Depth  = ifelse(is.null(dep), NA, dep),
    N_Patients = sum(vapply(seq_len(n_patients), function(i) term %in% patient_terms[[i]], logical(1))),
    stringsAsFactors = FALSE
  )
}
lookup_df <- do.call(rbind, lookup_rows)
lookup_df <- lookup_df[order(-lookup_df$N_Patients), ]
rownames(lookup_df) <- NULL

###############################################################################
# S9. WRITE EXCEL FILE
###############################################################################
cat("S9. Writing corrected Excel file...\n")

wb <- createWorkbook()

# --- Sheet 1: Statistical Analysis ---
addWorksheet(wb, "Statistical Analysis")
writeData(wb, "Statistical Analysis", primary_results)

# --- Sheet 2: Sensitivity Analysis ---
addWorksheet(wb, "Sensitivity Analysis")

title_sens <- paste0("Sensitivity Analysis: Full Cohort (N=", n_total,
                      ") vs Characterized Only (N=", n_characterized, ")")
subtitle_sens <- paste0("Excluding ", n_uncharacterized,
                         " patients without HPO terms. Base rescue rate: ",
                         n_char_rescued, "/", n_characterized, " = ",
                         round(n_char_rescued / n_characterized * 100, 1), "%.")

writeData(wb, "Sensitivity Analysis", data.frame(x = title_sens), startRow = 1, colNames = FALSE)
writeData(wb, "Sensitivity Analysis", data.frame(x = subtitle_sens), startRow = 2, colNames = FALSE)
writeData(wb, "Sensitivity Analysis", sens_df, startRow = 4)

# --- Sheet 3: Within-Neurology Analysis ---
addWorksheet(wb, "Within-Neurology Analysis")

title_neuro <- paste0("Within-Neurology Analysis: ", n_neuro,
                       " patients with \u22651 HPO term under Abnormality of the nervous system (HP:0000707)")
subtitle_neuro <- paste0("Baseline rescue rate: ", n_neuro_rescue, "/", n_neuro,
                          " = ", round(n_neuro_rescue / n_neuro * 100, 1),
                          "%. Each category tested against other neuro patients.")

writeData(wb, "Within-Neurology Analysis", data.frame(x = title_neuro), startRow = 1, colNames = FALSE)
writeData(wb, "Within-Neurology Analysis", data.frame(x = subtitle_neuro), startRow = 2, colNames = FALSE)
writeData(wb, "Within-Neurology Analysis", neuro_results, startRow = 4)

# --- Sheet 4: Rescue Rates ---
addWorksheet(wb, "Rescue Rates (All)")
writeData(wb, "Rescue Rates (All)", rescue_df)

# --- Sheet 5: Category Summary ---
addWorksheet(wb, "Category Summary")
writeData(wb, "Category Summary", cat_summary)

# --- Sheet 6: Patient Categories ---
addWorksheet(wb, "Patient Categories")
writeData(wb, "Patient Categories", pat_cat_df)

# --- Sheet 7: HPO Term Lookup ---
addWorksheet(wb, "HPO Term Lookup")
writeData(wb, "HPO Term Lookup", lookup_df)

saveWorkbook(wb, "patient_hpo_categories.xlsx", overwrite = TRUE)
cat("  Saved: patient_hpo_categories.xlsx\n")


###############################################################################
# S10. VERIFICATION PRINTOUT
###############################################################################
cat("\n")
cat("================================================================\n")
cat("VERIFICATION SUMMARY\n")
cat("================================================================\n")
cat(sprintf("Total patients:          %d\n", n_total))
cat(sprintf("Total rescued:           %d\n", n_rescued))
cat(sprintf("Characterized:           %d\n", n_characterized))
cat(sprintf("Uncharacterized:         %d\n", n_uncharacterized))
cat(sprintf("Neurological:            %d\n", n_neuro))
cat(sprintf("Neuro rescued:           %d\n", n_neuro_rescue))
cat(sprintf("Neuro rescue rate:       %.1f%%\n", n_neuro_rescue / n_neuro * 100))
cat(sprintf("\nCategories selected:     %d (min %d patients, depth >= %d)\n",
            length(primary_categories), MIN_CATEGORY_SIZE, DEPTH_CAP))
cat(sprintf("Primary Analysis:        %d tested, %d significant\n",
            nrow(primary_results), n_sig_primary))
cat(sprintf("Sensitivity Analysis:    %d tested, %d significant\n",
            nrow(char_results), n_sig_char))
cat(sprintf("Within-Neurology:        %d tested, %d significant\n",
            nrow(neuro_results), n_sig_neuro))

# Key numbers for manuscript cross-check
cat("\n--- Key numbers for manuscript cross-check ---\n")

# Neurodevelopmental delay in within-neurology
nd_idx <- which(neuro_results$Category_Name == "Neurodevelopmental delay")
if (length(nd_idx) > 0) {
  cat(sprintf("Neurodevelopmental delay (within-neuro):\n"))
  cat(sprintf("  Rescue rate IN:  %.1f%% (%d/%d)\n",
              neuro_results$Rescue_Rate_In[nd_idx] * 100,
              neuro_results$Rescued_In[nd_idx],
              neuro_results$N_In[nd_idx]))
  cat(sprintf("  Rescue rate OUT: %.1f%% (%d/%d)\n",
              neuro_results$Rescue_Rate_Out[nd_idx] * 100,
              neuro_results$Rescued_Out[nd_idx],
              neuro_results$N_Out[nd_idx]))
  cat(sprintf("  OR = %.4f, p = %.4e, q = %.6f\n",
              neuro_results$Odds_Ratio[nd_idx],
              neuro_results$P_value[nd_idx],
              neuro_results$Q_value[nd_idx]))
}

# Neurodevelopmental abnormality
nda_idx <- which(neuro_results$Category_Name == "Neurodevelopmental abnormality")
if (length(nda_idx) > 0) {
  cat(sprintf("Neurodevelopmental abnormality (within-neuro):\n"))
  cat(sprintf("  q = %.6f [%s]\n",
              neuro_results$Q_value[nda_idx],
              neuro_results$Significance[nda_idx]))
}

# Heart morphology in primary
hm_idx <- which(primary_results$HPO_ID == "HP:0001627")
if (length(hm_idx) > 0) {
  cat(sprintf("Abnormal heart morphology (primary):\n"))
  cat(sprintf("  N_In = %d, Rescue = %.1f%%, OR = %.4f, q = %.4f [%s]\n",
              primary_results$N_In[hm_idx],
              primary_results$Rescue_Rate_In[hm_idx] * 100,
              primary_results$Odds_Ratio[hm_idx],
              primary_results$Q_value[hm_idx],
              primary_results$Significance[hm_idx]))
}

# Seizure (borderline category)
sz_idx <- which(primary_results$HPO_ID == "HP:0001250")
if (length(sz_idx) > 0) {
  cat(sprintf("Seizure (primary):\n"))
  cat(sprintf("  q = %.6f [%s]\n",
              primary_results$Q_value[sz_idx],
              primary_results$Significance[sz_idx]))
}

# Nervous system
ns_idx <- which(primary_results$HPO_ID == "HP:0000707")
if (length(ns_idx) > 0) {
  cat(sprintf("Abnormality of the nervous system (primary):\n"))
  cat(sprintf("  N_In = %d\n", primary_results$N_In[ns_idx]))
}

cat("\n================================================================\n")
cat("DONE. All numbers generated from source data.\n")
cat("Categories selected algorithmically (not hardcoded).\n")
cat("================================================================\n")

# Close sink
sink()
cat(sprintf("\nOutput saved to: %s\n", output_file))
