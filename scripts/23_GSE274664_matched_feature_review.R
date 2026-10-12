# ============================================================
# Script 23
# Matched comparison of combined and hairpin-associated counts
# ============================================================
#
# Both analyses use the original Scripts 14 and 15 settings:
# - Same 12 samples
# - filterByExpr: min.count = 3, min.total.count = 10
# - TMM normalization after filtering
# - estimateDisp: robust = TRUE
# - glmQLFit: robust = TRUE
# - Responder versus Non_responder
# - FDR < 0.05 and absolute logFC >= 1
#
# Existing analysis files are NOT modified.
# Agreement does not establish that pooling is biologically valid.
# ============================================================

main <- function() {
  
  # ----------------------------------------------------------
  # 1. Paths and packages
  # ----------------------------------------------------------
  
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  
  objects <- file.path(root, "results", "R_objects")
  
  files <- c(
    Raw = file.path(
      root, "data", "raw",
      "GSE274664_CD19_EV_miRNA_readcounts.txt.gz"
    ),
    Combined = file.path(
      objects,
      "GSE274664_mature_miRNA_count_matrix.rds"
    ),
    Metadata = file.path(
      objects,
      "GSE274664_final_miRNA_metadata.rds"
    ),
    Previous = file.path(
      objects,
      "GSE274664_miRNA_DE_results.rds"
    )
  )
  
  if (any(!file.exists(files))) {
    stop(
      "Missing files:\n",
      paste(files[!file.exists(files)], collapse = "\n")
    )
  }
  
  for (pkg in c("edgeR", "statmod")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(
        "Required package is missing: ", pkg,
        ". Install it before running this script."
      )
    }
  }
  
  output_parent <- file.path(
    root, "data", "processed", "matched_feature_review"
  )
  
  dir.create(
    output_parent,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  output <- tempfile(
    pattern = paste0(
      "Review_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"
    ),
    tmpdir = output_parent
  )
  
  dir.create(output)
  
  # ----------------------------------------------------------
  # 2. Same final samples
  # ----------------------------------------------------------
  
  meta <- as.data.frame(readRDS(files[["Metadata"]]))
  previous <- as.data.frame(readRDS(files[["Previous"]]))
  combined <- as.matrix(readRDS(files[["Combined"]]))
  
  if (!all(c("Sample_ID", "Response_group") %in%
           names(meta))) {
    stop("Metadata columns are missing.")
  }
  
  meta$Sample_ID <- trimws(as.character(meta$Sample_ID))
  
  meta$Response_group <- factor(
    as.character(meta$Response_group),
    levels = c("Non_responder", "Responder")
  )
  
  if (
    anyNA(meta$Sample_ID) ||
    anyDuplicated(meta$Sample_ID) ||
    anyNA(meta$Response_group)
  ) {
    stop("Invalid sample identifiers or response groups.")
  }
  
  if (
    nrow(meta) != 12 ||
    sum(meta$Response_group == "Responder") != 8 ||
    sum(meta$Response_group == "Non_responder") != 4 ||
    "Hg19_699_CD19" %in% meta$Sample_ID
  ) {
    stop("Expected the final 12-sample set from Script 14.")
  }
  
  if (!all(meta$Sample_ID %in% colnames(combined))) {
    stop("Final samples are missing from the combined matrix.")
  }
  
  combined <- combined[, meta$Sample_ID, drop = FALSE]
  
  rownames(meta) <- meta$Sample_ID
  
  design <- model.matrix(
    ~ Response_group,
    data = meta
  )
  
  if (qr(design)$rank != ncol(design)) {
    stop("The design matrix is not full rank.")
  }
  
  # ----------------------------------------------------------
  # 3. Read hairpin-associated rows directly from the raw file
  # ----------------------------------------------------------
  
  cat("Reading the original count file...\n")
  
  con <- gzfile(files[["Raw"]], open = "rt")
  
  raw <- tryCatch(
    read.delim(
      con,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = ""
    ),
    finally = close(con)
  )
  
  if (ncol(raw) < 2 || anyDuplicated(names(raw))) {
    stop("Unexpected raw-file columns.")
  }
  
  feature <- trimws(as.character(raw[[1]]))
  title_row <- which(tolower(feature) == "title")
  
  if (length(title_row) != 1) {
    stop("Expected one title row.")
  }
  
  if (!all(meta$Sample_ID %in% names(raw))) {
    stop("Final samples are missing from the raw file.")
  }
  
  labels <- vapply(
    meta$Sample_ID,
    function(id) trimws(as.character(raw[[id]][title_row])),
    character(1)
  )
  
  raw_groups <- ifelse(
    grepl("^NR[0-9]+$", labels),
    "Non_responder",
    ifelse(
      grepl("^R[0-9]+$", labels),
      "Responder",
      NA_character_
    )
  )
  
  if (
    anyNA(raw_groups) ||
    any(raw_groups != as.character(meta$Response_group))
  ) {
    stop("Raw labels and saved response groups disagree.")
  }
  
  selected <- grepl("^hsa-", feature) &
    grepl("_hairpin$", feature)
  
  if (!any(selected)) {
    stop("No human hairpin-associated rows were found.")
  }
  
  hairpin <- matrix(
    NA_real_,
    nrow = sum(selected),
    ncol = nrow(meta),
    dimnames = list(
      sub("_hairpin$", "", feature[selected]),
      meta$Sample_ID
    )
  )
  
  for (j in seq_len(nrow(meta))) {
    hairpin[, j] <- suppressWarnings(
      as.numeric(
        as.character(raw[[meta$Sample_ID[j]]][selected])
      )
    )
  }
  
  # No mature rows are added to the hairpin matrix.
  validate_counts <- function(x, label) {
    
    if (
      !is.numeric(x) ||
      anyNA(x) ||
      any(!is.finite(x)) ||
      any(x < 0)
    ) {
      stop("Invalid counts in: ", label)
    }
    
    if (any(abs(x - round(x)) > 1e-8)) {
      stop("Non-integer read counts in: ", label)
    }
    
    if (
      is.null(rownames(x)) ||
      anyDuplicated(rownames(x)) ||
      anyDuplicated(colnames(x))
    ) {
      stop("Invalid or duplicated identifiers in: ", label)
    }
    
    if (any(colSums(x) <= 0)) {
      stop("Zero library in: ", label)
    }
  }
  
  validate_counts(combined, "Combined matrix")
  validate_counts(hairpin, "Hairpin-associated matrix")
  
  # ----------------------------------------------------------
  # 4. One identical analysis function for both matrices
  # ----------------------------------------------------------
  
  analyze <- function(counts, label) {
    
    cat("\nAnalyzing:", label, "\n")
    
    y <- edgeR::DGEList(counts = counts)
    
    # EXACT filtering call used in Script 14.
    keep <- edgeR::filterByExpr(
      y,
      group = meta$Response_group,
      min.count = 3,
      min.total.count = 10
    )
    
    filter_table <- data.frame(
      miRNA = rownames(counts),
      Retained = keep
    )
    
    if (sum(keep) < 10) {
      stop("Fewer than 10 features passed filtering: ", label)
    }
    
    y <- y[keep, , keep.lib.sizes = FALSE]
    
    if (any(colSums(y$counts) <= 0)) {
      stop("A zero filtered library was generated: ", label)
    }
    
    if ("normLibSizes" %in% getNamespaceExports("edgeR")) {
      y <- edgeR::normLibSizes(y, method = "TMM")
    } else {
      y <- edgeR::calcNormFactors(y, method = "TMM")
    }
    
    if (
      any(!is.finite(y$samples$norm.factors)) ||
      any(y$samples$norm.factors <= 0)
    ) {
      stop("Invalid normalization factors: ", label)
    }
    
    # EXACT model settings used in Script 15.
    y <- edgeR::estimateDisp(
      y,
      design = design,
      robust = TRUE
    )
    
    fit <- edgeR::glmQLFit(
      y,
      design = design,
      robust = TRUE
    )
    
    test <- edgeR::glmQLFTest(
      fit,
      coef = "Response_groupResponder"
    )
    
    tab <- edgeR::topTags(
      test,
      n = Inf,
      sort.by = "PValue"
    )$table
    
    result <- data.frame(
      miRNA = rownames(tab),
      tab,
      row.names = NULL,
      check.names = FALSE
    )
    
    result$Primary_significant <- (
      result$FDR < 0.05 & abs(result$logFC) >= 1
    )
    
    result$Nominal_significant <- (
      result$PValue < 0.05 & abs(result$logFC) >= 1
    )
    
    cat("Features tested:", nrow(result), "\n")
    cat(
      "Primary significant:",
      sum(result$Primary_significant), "\n"
    )
    
    list(
      Results = result,
      Filtering = filter_table,
      DGEList = y,
      LogCPM = edgeR::cpm(y, log = TRUE, prior.count = 2)
    )
  }
  
  combined_review <- analyze(
    combined,
    "Previous combined feature definition"
  )
  
  hairpin_review <- analyze(
    hairpin,
    "Hairpin-associated only"
  )
  
  # ----------------------------------------------------------
  # 5. Compare the four candidates
  # ----------------------------------------------------------
  
  candidates <- c(
    "hsa-miR-509-3p",
    "hsa-miR-26b-5p",
    "hsa-miR-196b-5p",
    "hsa-miR-375-3p"
  )
  
  a <- combined_review$Results
  b <- hairpin_review$Results
  
  ia <- match(candidates, a$miRNA)
  ib <- match(candidates, b$miRNA)
  
  comparison <- data.frame(
    miRNA = candidates,
    Combined_tested = !is.na(ia),
    Combined_logFC = a$logFC[ia],
    Combined_PValue = a$PValue[ia],
    Combined_FDR = a$FDR[ia],
    Hairpin_tested = !is.na(ib),
    Hairpin_logFC = b$logFC[ib],
    Hairpin_PValue = b$PValue[ib],
    Hairpin_FDR = b$FDR[ib],
    stringsAsFactors = FALSE
  )
  
  comparison$Same_direction <- NA
  
  both <- comparison$Combined_tested &
    comparison$Hairpin_tested
  
  comparison$Same_direction[both] <- (
    sign(comparison$Combined_logFC[both]) ==
      sign(comparison$Hairpin_logFC[both])
  )
  
  comparison$Review_status <- vapply(
    seq_len(nrow(comparison)),
    function(i) {
      
      if (!comparison$Hairpin_tested[i]) {
        return("Hairpin_not_retained_for_testing")
      }
      
      if (!comparison$Combined_tested[i]) {
        return("Combined_not_retained_for_testing")
      }
      
      if (!isTRUE(comparison$Same_direction[i])) {
        return("Direction_changed")
      }
      
      if (
        comparison$Hairpin_FDR[i] < 0.05 &&
        abs(comparison$Hairpin_logFC[i]) >= 1
      ) {
        return("Same_direction_primary_support")
      }
      
      if (
        comparison$Hairpin_PValue[i] < 0.05 &&
        abs(comparison$Hairpin_logFC[i]) >= 1
      ) {
        return("Same_direction_nominal_support_only")
      }
      
      "Same_direction_without_nominal_support"
    },
    character(1)
  )
  
  # ----------------------------------------------------------
  # 6. Check reproduction of the previous combined results
  # ----------------------------------------------------------
  
  if (!all(c("miRNA", "logFC", "FDR") %in%
           names(previous))) {
    stop("Previous results have unexpected columns.")
  }
  
  if (anyDuplicated(previous$miRNA)) {
    stop("Previous results contain duplicated identifiers.")
  }
  
  reproduction <- data.frame(
    Same_tested_features = setequal(a$miRNA, previous$miRNA),
    Max_logFC_difference = NA_real_,
    Max_FDR_difference = NA_real_
  )
  
  common <- intersect(a$miRNA, previous$miRNA)
  
  if (length(common) > 0) {
    
    i_new <- match(common, a$miRNA)
    i_old <- match(common, previous$miRNA)
    
    reproduction$Max_logFC_difference <- max(
      abs(a$logFC[i_new] - previous$logFC[i_old])
    )
    
    reproduction$Max_FDR_difference <- max(
      abs(a$FDR[i_new] - previous$FDR[i_old])
    )
  }
  
  # ----------------------------------------------------------
  # 7. Save separate outputs
  # ----------------------------------------------------------
  
  write.csv(
    a,
    file.path(output, "Combined_matched_DE_results.csv"),
    row.names = FALSE
  )
  
  write.csv(
    b,
    file.path(output, "Hairpin_matched_DE_results.csv"),
    row.names = FALSE
  )
  
  write.csv(
    comparison,
    file.path(output, "Four_candidate_matched_comparison.csv"),
    row.names = FALSE
  )
  
  write.csv(
    reproduction,
    file.path(output, "Previous_result_reproduction_check.csv"),
    row.names = FALSE
  )
  
  write.csv(
    combined_review$Filtering,
    file.path(output, "Combined_filtering.csv"),
    row.names = FALSE
  )
  
  write.csv(
    hairpin_review$Filtering,
    file.path(output, "Hairpin_filtering.csv"),
    row.names = FALSE
  )
  
  saveRDS(
    list(
      Combined = combined_review,
      Hairpin = hairpin_review,
      Candidate_comparison = comparison,
      Reproduction_check = reproduction,
      Metadata = meta,
      Design = design,
      Settings = list(
        min.count = 3,
        min.total.count = 10,
        Normalization = "TMM",
        Robust = TRUE,
        FDR_threshold = 0.05,
        Absolute_logFC_threshold = 1
      )
    ),
    file.path(output, "GSE274664_matched_feature_review.rds")
  )
  
  writeLines(
    capture.output(sessionInfo()),
    file.path(output, "R_session_information.txt")
  )
  
  # ----------------------------------------------------------
  # 8. Summary
  # ----------------------------------------------------------
  
  cat(
    "\n============================================\n",
    "Matched feature review completed\n",
    "============================================\n",
    sep = ""
  )
  
  cat("\nPrevious combined-result reproduction:\n")
  print(reproduction, row.names = FALSE, digits = 6)
  
  cat("\nFour-candidate matched comparison:\n")
  print(comparison, row.names = FALSE, digits = 4)
  
  cat(
    "\nImportant:\n",
    "Both analyses used the same original settings.\n",
    "No previous files were modified.\n",
    "Agreement supports sensitivity robustness, not pooling validity.\n",
    "Hairpin labels do not automatically establish precursor abundance.\n",
    "The original counting pipeline still determines whether pooling is valid.\n",
    sep = ""
  )
  
  cat(
    "\nOutput folder:\n",
    normalizePath(output, winslash = "/", mustWork = TRUE),
    "\n",
    sep = ""
  )
  
  invisible(comparison)
}

main()