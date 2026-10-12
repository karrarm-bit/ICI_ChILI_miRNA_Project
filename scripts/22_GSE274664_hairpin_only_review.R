# ============================================================
# Script 22
# GSE274664 hairpin-associated-only sensitivity review
# ============================================================
#
# Purpose:
# - Read the ORIGINAL count file.
# - Analyze human hairpin-associated features only.
# - Use the same 12 samples retained in the previous analysis.
# - Compare the four previous candidates with this analysis.
# - Preserve all previous files.
#
# Important:
# Hairpin-associated labels are database/count-file categories.
# They are not automatically measurements of precursor abundance.
# This analysis does not prove that mature and hairpin counts
# can be pooled without double counting.
# ============================================================

main <- function() {
  
  # ----------------------------------------------------------
  # 1. Paths
  # ----------------------------------------------------------
  
  project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  
  raw_file <- file.path(
    project_dir,
    "data",
    "raw",
    "GSE274664_CD19_EV_miRNA_readcounts.txt.gz"
  )
  
  object_dir <- file.path(
    project_dir,
    "results",
    "R_objects"
  )
  
  metadata_file <- file.path(
    object_dir,
    "GSE274664_final_miRNA_metadata.rds"
  )
  
  baseline_file <- file.path(
    object_dir,
    "GSE274664_miRNA_DE_results.rds"
  )
  
  review_dir <- file.path(
    project_dir,
    "data",
    "processed",
    "hairpin_only_review"
  )
  
  if (!dir.exists(project_dir)) {
    stop("Project directory does not exist.")
  }
  
  if (!file.exists(raw_file)) {
    stop("Original count file is missing: ", raw_file)
  }
  
  if (!file.exists(metadata_file)) {
    stop("Final sample metadata is missing: ", metadata_file)
  }
  
  if (!file.exists(baseline_file)) {
    stop("Previous differential-expression results are missing.")
  }
  
  dir.create(
    review_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # Create a separate folder for every run.
  run_dir <- tempfile(
    pattern = paste0(
      "Review_",
      format(Sys.time(), "%Y%m%d_%H%M%S"),
      "_"
    ),
    tmpdir = review_dir
  )
  
  dir.create(run_dir)
  
  # ----------------------------------------------------------
  # 2. Required packages
  # ----------------------------------------------------------
  
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages(
      "BiocManager",
      repos = "https://cloud.r-project.org"
    )
  }
  
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    BiocManager::install(
      "edgeR",
      ask = FALSE,
      update = FALSE
    )
  }
  
  if (!requireNamespace("statmod", quietly = TRUE)) {
    install.packages(
      "statmod",
      repos = "https://cloud.r-project.org"
    )
  }
  
  if (
    !requireNamespace("edgeR", quietly = TRUE) ||
    !requireNamespace("statmod", quietly = TRUE)
  ) {
    stop("edgeR or statmod could not be loaded.")
  }
  
  # ----------------------------------------------------------
  # 3. Read metadata and previous results
  # ----------------------------------------------------------
  
  metadata <- as.data.frame(
    readRDS(metadata_file)
  )
  
  baseline <- as.data.frame(
    readRDS(baseline_file)
  )
  
  required_metadata <- c(
    "Sample_ID",
    "Response_group"
  )
  
  if (!all(required_metadata %in% names(metadata))) {
    stop(
      "Metadata must contain Sample_ID and Response_group."
    )
  }
  
  if (!all(c("miRNA", "logFC", "FDR") %in%
           names(baseline))) {
    stop(
      "Previous results must contain miRNA, logFC and FDR."
    )
  }
  
  metadata$Sample_ID <- trimws(
    as.character(metadata$Sample_ID)
  )
  
  metadata$Response_group <- as.character(
    metadata$Response_group
  )
  
  if (
    anyNA(metadata$Sample_ID) ||
    anyDuplicated(metadata$Sample_ID)
  ) {
    stop("Missing or duplicated metadata sample identifiers.")
  }
  
  if (
    anyNA(metadata$Response_group) ||
    !all(metadata$Response_group %in%
         c("Non_responder", "Responder"))
  ) {
    stop("Unexpected or missing response groups.")
  }
  
  if ("Hg19_699_CD19" %in% metadata$Sample_ID) {
    stop(
      "NR5 remains in metadata. ",
      "Use the final metadata produced by Script 14."
    )
  }
  
  if (
    nrow(metadata) != 12 ||
    sum(metadata$Response_group == "Responder") != 8 ||
    sum(metadata$Response_group == "Non_responder") != 4
  ) {
    stop(
      "Expected the previous final sample set: ",
      "12 samples, 8 Responders and 4 Non-responders."
    )
  }
  
  # ----------------------------------------------------------
  # 4. Read the original count file
  # ----------------------------------------------------------
  
  cat("Reading the ORIGINAL count file...\n")
  
  connection <- gzfile(
    raw_file,
    open = "rt"
  )
  
  raw <- tryCatch(
    read.delim(
      connection,
      header = TRUE,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = ""
    ),
    finally = close(connection)
  )
  
  if (ncol(raw) < 2) {
    stop("The original file has an unexpected structure.")
  }
  
  if (anyDuplicated(names(raw))) {
    stop("Duplicated column names were found in the raw file.")
  }
  
  feature_id <- trimws(
    as.character(raw[[1]])
  )
  
  title_rows <- which(
    tolower(feature_id) == "title"
  )
  
  if (length(title_rows) != 1) {
    stop("Expected exactly one title row in the original file.")
  }
  
  missing_samples <- setdiff(
    metadata$Sample_ID,
    names(raw)
  )
  
  if (length(missing_samples) > 0) {
    stop(
      "Samples absent from the raw file: ",
      paste(missing_samples, collapse = ", ")
    )
  }
  
  # Check that raw labels agree with saved response groups.
  raw_labels <- vapply(
    metadata$Sample_ID,
    function(sample_id) {
      trimws(
        as.character(raw[[sample_id]][title_rows])
      )
    },
    character(1)
  )
  
  raw_groups <- ifelse(
    grepl("^NR[0-9]+$", raw_labels),
    "Non_responder",
    ifelse(
      grepl("^R[0-9]+$", raw_labels),
      "Responder",
      NA_character_
    )
  )
  
  if (
    anyNA(raw_groups) ||
    any(raw_groups != metadata$Response_group)
  ) {
    stop(
      "Raw sample labels do not agree with final metadata."
    )
  }
  
  metadata$Sample_label <- raw_labels
  
  # ----------------------------------------------------------
  # 5. Select human hairpin-associated rows
  # ----------------------------------------------------------
  
  selected_rows <- grepl(
    "^hsa-",
    feature_id
  ) & grepl(
    "_hairpin$",
    feature_id
  )
  
  if (!any(selected_rows)) {
    stop("No human hairpin-associated rows were found.")
  }
  
  hairpin_ids <- feature_id[selected_rows]
  
  cat(
    "Human hairpin-associated rows:",
    length(hairpin_ids),
    "\n"
  )
  
  counts <- matrix(
    NA_real_,
    nrow = length(hairpin_ids),
    ncol = nrow(metadata),
    dimnames = list(
      hairpin_ids,
      metadata$Sample_ID
    )
  )
  
  for (j in seq_len(nrow(metadata))) {
    
    sample_id <- metadata$Sample_ID[j]
    
    original_values <- trimws(
      as.character(raw[[sample_id]][selected_rows])
    )
    
    numeric_values <- suppressWarnings(
      as.numeric(original_values)
    )
    
    if (
      anyNA(numeric_values) ||
      any(!is.finite(numeric_values)) ||
      any(numeric_values < 0)
    ) {
      stop(
        "Invalid hairpin count values in sample: ",
        sample_id
      )
    }
    
    if (
      any(abs(numeric_values -
              round(numeric_values)) > 1e-8)
    ) {
      stop(
        "Non-integer values were found in a read-count file."
      )
    }
    
    counts[, j] <- numeric_values
  }
  
  # Remove only the category suffix.
  # Do not remove the 3p/5p arm designation.
  hairpin_names <- sub(
    "_hairpin$",
    "",
    hairpin_ids
  )
  
  if (anyDuplicated(hairpin_names)) {
    stop(
      "Duplicated hairpin-associated names were found. ",
      "Inspect these rows before aggregation."
    )
  }
  
  rownames(counts) <- hairpin_names
  
  # No mature rows are added to this matrix.
  stopifnot(
    identical(
      colnames(counts),
      metadata$Sample_ID
    )
  )
  
  # ----------------------------------------------------------
  # 6. Sample count summary
  # ----------------------------------------------------------
  
  sample_statistics <- data.frame(
    Sample_ID = metadata$Sample_ID,
    Sample_label = metadata$Sample_label,
    Response_group = metadata$Response_group,
    Hairpin_associated_library_size = colSums(counts),
    Detected_hairpin_associated_features = colSums(
      counts > 0
    ),
    stringsAsFactors = FALSE
  )
  
  cat("\nSample statistics:\n")
  
  print(
    sample_statistics,
    row.names = FALSE
  )
  
  write.csv(
    sample_statistics,
    file.path(
      run_dir,
      "Hairpin_only_sample_statistics.csv"
    ),
    row.names = FALSE
  )
  
  if (any(colSums(counts) <= 0)) {
    stop(
      "A zero hairpin-associated library was detected. ",
      "No additional samples were automatically excluded."
    )
  }
  
  # Remove completely unexpressed features.
  counts <- counts[
    rowSums(counts) > 0,
    ,
    drop = FALSE
  ]
  
  saveRDS(
    counts,
    file.path(
      run_dir,
      "Hairpin_only_counts_before_expression_filtering.rds"
    )
  )
  
  # ----------------------------------------------------------
  # 7. Statistical design
  # ----------------------------------------------------------
  
  metadata$Response_group <- factor(
    metadata$Response_group,
    levels = c(
      "Non_responder",
      "Responder"
    )
  )
  
  rownames(metadata) <- metadata$Sample_ID
  
  design <- model.matrix(
    ~ Response_group,
    data = metadata
  )
  
  if (qr(design)$rank != ncol(design)) {
    stop("The statistical design is not full rank.")
  }
  
  coefficient <- match(
    "Response_groupResponder",
    colnames(design)
  )
  
  if (is.na(coefficient)) {
    stop("Responder comparison coefficient was not found.")
  }
  
  # ----------------------------------------------------------
  # 8. Filtering and TMM normalization
  # ----------------------------------------------------------
  
  cat("\nFiltering hairpin-associated features...\n")
  
  y <- edgeR::DGEList(
    counts = counts
  )
  
  # Explicit standard filtering settings.
  # If the previous analysis used different settings,
  # differences may reflect both feature definition
  # and filtering.
  keep <- edgeR::filterByExpr(
    y,
    design = design,
    min.count = 10,
    min.total.count = 15,
    large.n = 10,
    min.prop = 0.7
  )
  
  filtering_table <- data.frame(
    miRNA = rownames(y),
    Retained_for_testing = keep,
    stringsAsFactors = FALSE
  )
  
  write.csv(
    filtering_table,
    file.path(
      run_dir,
      "Hairpin_only_expression_filter.csv"
    ),
    row.names = FALSE
  )
  
  cat("Features before filtering:", nrow(y), "\n")
  cat("Features retained:", sum(keep), "\n")
  
  if (sum(keep) < 2) {
    stop("Too few features remain for differential analysis.")
  }
  
  y <- y[
    keep,
    ,
    keep.lib.sizes = FALSE
  ]
  
  if (any(colSums(y$counts) <= 0)) {
    stop(
      "A library became zero after filtering. ",
      "No normalization was attempted."
    )
  }
  
  y <- edgeR::calcNormFactors(
    y,
    method = "TMM"
  )
  
  # ----------------------------------------------------------
  # 9. Differential-expression analysis
  # ----------------------------------------------------------
  
  cat("\nFitting the hairpin-associated-only model...\n")
  
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
    coef = coefficient
  )
  
  results <- edgeR::topTags(
    test,
    n = Inf,
    sort.by = "PValue"
  )$table
  
  results$miRNA <- rownames(results)
  
  results$Primary_significant <- (
    results$FDR < 0.05 &
      abs(results$logFC) >= 1
  )
  
  results$Nominal_candidate <- (
    results$PValue < 0.05 &
      abs(results$logFC) >= 1
  )
  
  results$Direction <- ifelse(
    results$logFC > 0,
    "Higher_in_Responder",
    ifelse(
      results$logFC < 0,
      "Higher_in_Non_responder",
      "No_change"
    )
  )
  
  results <- results[
    ,
    c(
      "miRNA",
      "logFC",
      "logCPM",
      "F",
      "PValue",
      "FDR",
      "Primary_significant",
      "Nominal_candidate",
      "Direction"
    )
  ]
  
  rownames(results) <- NULL
  
  # ----------------------------------------------------------
  # 10. Compare the four previous candidates
  # ----------------------------------------------------------
  
  candidates <- c(
    "hsa-miR-509-3p",
    "hsa-miR-26b-5p",
    "hsa-miR-196b-5p",
    "hsa-miR-375-3p"
  )
  
  baseline$miRNA <- trimws(
    as.character(baseline$miRNA)
  )
  
  if (anyDuplicated(baseline$miRNA)) {
    stop("Previous results contain duplicated miRNA names.")
  }
  
  baseline_index <- match(
    candidates,
    baseline$miRNA
  )
  
  if (anyNA(baseline_index)) {
    stop(
      "At least one previous candidate is missing ",
      "from baseline results."
    )
  }
  
  hairpin_index <- match(
    candidates,
    results$miRNA
  )
  
  comparison <- data.frame(
    miRNA = candidates,
    Previous_logFC = baseline$logFC[baseline_index],
    Previous_FDR = baseline$FDR[baseline_index],
    Present_in_hairpin_matrix = (
      candidates %in% rownames(counts)
    ),
    Tested_in_hairpin_analysis = !is.na(hairpin_index),
    Hairpin_logFC = results$logFC[hairpin_index],
    Hairpin_PValue = results$PValue[hairpin_index],
    Hairpin_FDR = results$FDR[hairpin_index],
    stringsAsFactors = FALSE
  )
  
  comparison$Same_direction <- NA
  
  tested <- comparison$Tested_in_hairpin_analysis
  
  comparison$Same_direction[tested] <- (
    sign(comparison$Previous_logFC[tested]) ==
      sign(comparison$Hairpin_logFC[tested])
  )
  
  comparison$Review_status <- vapply(
    seq_len(nrow(comparison)),
    function(i) {
      
      if (!comparison$Present_in_hairpin_matrix[i]) {
        return("Absent_or_unexpressed_in_hairpin_matrix")
      }
      
      if (!comparison$Tested_in_hairpin_analysis[i]) {
        return("Not_retained_after_expression_filtering")
      }
      
      if (!isTRUE(comparison$Same_direction[i])) {
        return("Direction_changed")
      }
      
      if (
        comparison$Hairpin_FDR[i] < 0.05 &&
        abs(comparison$Hairpin_logFC[i]) >= 1
      ) {
        return("Same_direction_and_primary_threshold_met")
      }
      
      if (
        comparison$Hairpin_PValue[i] < 0.05 &&
        abs(comparison$Hairpin_logFC[i]) >= 1
      ) {
        return("Same_direction_with_nominal_support_only")
      }
      
      return("Same_direction_without_nominal_support")
    },
    character(1)
  )
  
  # ----------------------------------------------------------
  # 11. Save separate review outputs
  # ----------------------------------------------------------
  
  write.csv(
    results,
    file.path(
      run_dir,
      "Hairpin_only_differential_expression.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    comparison,
    file.path(
      run_dir,
      "Four_candidate_hairpin_comparison.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    list(
      Metadata = metadata,
      Design = design,
      DGEList = y,
      Results = results,
      Candidate_comparison = comparison,
      Sample_statistics = sample_statistics,
      Filtering = filtering_table,
      Analysis_settings = list(
        Feature_category = "Human_hairpin_associated_only",
        Normalization = "TMM",
        Robust_estimation = TRUE,
        min.count = 10,
        min.total.count = 15,
        large.n = 10,
        min.prop = 0.7,
        Comparison = "Responder_vs_Non_responder"
      )
    ),
    file.path(
      run_dir,
      "GSE274664_hairpin_only_review.rds"
    )
  )
  
  writeLines(
    capture.output(sessionInfo()),
    file.path(
      run_dir,
      "R_session_information.txt"
    )
  )
  
  # ----------------------------------------------------------
  # 12. Simple candidate comparison figure
  # ----------------------------------------------------------
  
  plot_file <- file.path(
    run_dir,
    "Four_candidate_logFC_comparison.pdf"
  )
  
  grDevices::pdf(
    plot_file,
    width = 8,
    height = 5
  )
  
  tryCatch({
    
    plot_values <- rbind(
      Previous_combined = comparison$Previous_logFC,
      Hairpin_only = comparison$Hairpin_logFC
    )
    
    colnames(plot_values) <- sub(
      "^hsa-",
      "",
      comparison$miRNA
    )
    
    finite_values <- plot_values[
      is.finite(plot_values)
    ]
    
    plot_range <- range(
      c(0, finite_values)
    )
    
    padding <- max(
      1,
      diff(plot_range) * 0.15
    )
    
    graphics::par(
      mar = c(7, 5, 4, 2)
    )
    
    graphics::barplot(
      plot_values,
      beside = TRUE,
      col = c("#777777", "#377EB8"),
      ylim = plot_range + c(-padding, padding),
      las = 2,
      ylab = "log2 fold change",
      main = "Candidate comparison",
      sub = "Positive values indicate higher abundance in Responders"
    )
    
    graphics::abline(
      h = 0,
      col = "#000000"
    )
    
    graphics::legend(
      "topright",
      legend = c(
        "Previous combined analysis",
        "Hairpin-associated only"
      ),
      fill = c("#777777", "#377EB8"),
      bty = "n",
      cex = 0.85
    )
    
  }, finally = grDevices::dev.off())
  
  # ----------------------------------------------------------
  # 13. Console summary
  # ----------------------------------------------------------
  
  cat(
    "\n============================================\n",
    "Hairpin-associated-only review completed\n",
    "============================================\n",
    sep = ""
  )
  
  cat("Samples analyzed:", nrow(metadata), "\n")
  
  cat(
    "Hairpin-associated features tested:",
    nrow(results),
    "\n"
  )
  
  cat(
    "Features meeting FDR and fold-change thresholds:",
    sum(results$Primary_significant),
    "\n"
  )
  
  cat(
    "Common dispersion:",
    y$common.dispersion,
    "\n"
  )
  
  cat("\nFour-candidate comparison:\n")
  
  print(
    comparison,
    row.names = FALSE,
    digits = 4
  )
  
  cat(
    "\nInterpretation:\n",
    "Same direction and primary threshold met: ",
    "candidate remains supported in this sensitivity analysis.\n",
    "Nominal support only: exploratory evidence, not FDR support.\n",
    "Not tested: absence of a test is not evidence of no effect.\n",
    "Direction changed: the earlier interpretation needs review.\n",
    sep = ""
  )
  
  cat(
    "\nImportant:\n",
    "This is not independent validation.\n",
    "This does not prove that mature/hairpin pooling is valid.\n",
    "Filtering and robust-model settings may differ from ",
    "the previous script and can affect the comparison.\n",
    "No previous results or Word documents were modified.\n",
    sep = ""
  )
  
  cat(
    "\nOutput folder:\n",
    normalizePath(
      run_dir,
      winslash = "/",
      mustWork = TRUE
    ),
    "\n",
    sep = ""
  )
  
  invisible(comparison)
}

main()