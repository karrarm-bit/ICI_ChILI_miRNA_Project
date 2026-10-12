# ==============================================================================
# Script 20: Test enrichment of miRNA targets among primary mRNA DEGs
#
# Background: unique, valid gene symbols with finite mRNA test results.
# Main tests: three exact-name miRNA target sets.
# Provisional miR-375 alias: separate exploratory test.
#
# This script performs no online queries and requires no additional packages.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------------------------

mrna_file <- "results/R_objects/GSE287540_DE_results.rds"

integration_file <- paste0(
  "results/R_objects/",
  "GSE274664_GSE287540_target_integration.rds"
)

if (!all(file.exists(c(mrna_file, integration_file)))) {
  stop("Input files were not found. Open the correct RStudio project.")
}

mrna_results <- readRDS(mrna_file)
integration_results <- readRDS(integration_file)

required_columns <- c("Gene_symbol", "logFC", "PValue", "FDR")

if (!is.data.frame(mrna_results) ||
    !all(required_columns %in% names(mrna_results))) {
  stop("The mRNA results table has an unexpected structure.")
}

if (!all(vapply(
  mrna_results[, c("logFC", "PValue", "FDR"), drop = FALSE],
  is.numeric,
  logical(1)
))) {
  stop("mRNA logFC, PValue, and FDR must be numeric.")
}

# ------------------------------------------------------------------------------
# 2. Define the tested-gene background and DEG set
# ------------------------------------------------------------------------------

mrna_results$Match_symbol <- toupper(
  trimws(as.character(mrna_results$Gene_symbol))
)

valid_test <- (
  !is.na(mrna_results$Match_symbol) &
    nzchar(mrna_results$Match_symbol) &
    is.finite(mrna_results$logFC) &
    is.finite(mrna_results$PValue) &
    is.finite(mrna_results$FDR)
)

tested_results <- mrna_results[
  valid_test,
  ,
  drop = FALSE
]

background_genes <- sort(unique(tested_results$Match_symbol))

deg_rows <- (
  tested_results$FDR < 0.05 &
    abs(tested_results$logFC) >= 1
)

deg_genes <- sort(unique(
  tested_results$Match_symbol[deg_rows]
))

N <- length(background_genes)
D <- length(deg_genes)

if (N == 0L || D == 0L || D >= N) {
  stop(
    "A valid background containing both DEG and non-DEG genes is required. ",
    "Confirm that the mRNA file contains all tested genes, not only DEGs."
  )
}

cat("Background gene symbols:", N, "\n")
cat("Primary DEG symbols:", D, "\n")

# Duplicate symbols are counted once.
duplicate_symbol_rows <- (
  nrow(tested_results) - length(background_genes)
)

cat("Additional rows sharing a gene symbol:", duplicate_symbol_rows, "\n")

# ------------------------------------------------------------------------------
# 3. Function for one-sided Fisher overlap-enrichment testing
# ------------------------------------------------------------------------------

test_target_set <- function(target_table, miRNA_name, status) {
  
  required_target_columns <- c("miRNA", "Match_symbol")
  
  if (!all(required_target_columns %in% names(target_table))) {
    stop("Required target-table columns are missing.")
  }
  
  selected <- target_table[
    !is.na(target_table$miRNA) &
      target_table$miRNA == miRNA_name,
    ,
    drop = FALSE
  ]
  
  targets <- toupper(trimws(as.character(selected$Match_symbol)))
  
  targets <- sort(unique(
    targets[!is.na(targets) & nzchar(targets)]
  ))
  
  # Only targets belonging to the tested-gene background are eligible.
  eligible_targets <- intersect(targets, background_genes)
  
  K <- length(eligible_targets)
  overlap_genes <- intersect(eligible_targets, deg_genes)
  A <- length(overlap_genes)
  
  # Contingency table:
  #                         DEG       Non-DEG
  # Target                  A         K - A
  # Not target              D - A     N - K - D + A
  
  contingency_table <- matrix(
    c(
      A,
      K - A,
      D - A,
      N - K - D + A
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      Target_status = c("Target", "Not_target"),
      DEG_status = c("DEG", "Non_DEG")
    )
  )
  
  expected_overlap <- D * K / N
  
  if (K == 0L || K == N) {
    
    p_value <- NA_real_
    odds_ratio <- NA_real_
    
  } else {
    
    fisher_result <- fisher.test(
      contingency_table,
      alternative = "greater"
    )
    
    p_value <- fisher_result$p.value
    odds_ratio <- unname(fisher_result$estimate)
  }
  
  fold_enrichment <- if (expected_overlap > 0) {
    A / expected_overlap
  } else {
    NA_real_
  }
  
  data.frame(
    miRNA = miRNA_name,
    Mapping_status = status,
    Background_genes = N,
    DEG_genes = D,
    Total_database_targets = length(targets),
    Targets_in_background = K,
    Observed_overlap = A,
    Expected_overlap = expected_overlap,
    Fold_enrichment = fold_enrichment,
    Fisher_odds_ratio = odds_ratio,
    Enrichment_PValue = p_value,
    Overlap_genes = if (A > 0L) {
      paste(sort(overlap_genes), collapse = "; ")
    } else {
      ""
    },
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# 4. Test the three main target sets
# ------------------------------------------------------------------------------

main_targets <- integration_results$Main_target_genes

main_miRNAs <- c(
  "hsa-miR-509-3p",
  "hsa-miR-26b-5p",
  "hsa-miR-196b-5p"
)

main_results <- do.call(
  rbind,
  lapply(
    main_miRNAs,
    function(miRNA_name) {
      test_target_set(
        main_targets,
        miRNA_name,
        "Exact_name_query"
      )
    }
  )
)

# Correct across the three main miRNA overlap tests.
main_results$Enrichment_FDR <- p.adjust(
  main_results$Enrichment_PValue,
  method = "BH"
)

# ------------------------------------------------------------------------------
# 5. Test the provisional miR-375 set separately
# ------------------------------------------------------------------------------

provisional_results <- test_target_set(
  integration_results$Provisional_375_targets,
  "hsa-miR-375-3p",
  "Alias_pending_verification"
)

# Not included in the main multiple-testing family.
provisional_results$Enrichment_FDR <- NA_real_

# ------------------------------------------------------------------------------
# 6. Save outputs
# ------------------------------------------------------------------------------

csv_directory <- "data/processed/target_integration"

dir.create(
  csv_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  main_results,
  file.path(
    csv_directory,
    "Main_miRNA_target_overlap_enrichment.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  provisional_results,
  file.path(
    csv_directory,
    "Provisional_miR375_target_overlap_enrichment.csv"
  ),
  row.names = FALSE,
  na = ""
)

saveRDS(
  list(
    Main_tests = main_results,
    Provisional_test = provisional_results,
    Background_genes = background_genes,
    DEG_genes = deg_genes
  ),
  "results/R_objects/miRNA_target_overlap_enrichment.rds"
)

# ------------------------------------------------------------------------------
# 7. Completion report
# ------------------------------------------------------------------------------

cat("\n============================================\n")
cat("Target-overlap enrichment testing completed\n")
cat("============================================\n")

display_columns <- c(
  "miRNA",
  "Targets_in_background",
  "Observed_overlap",
  "Expected_overlap",
  "Fold_enrichment",
  "Enrichment_PValue",
  "Enrichment_FDR"
)

cat("\nMain enrichment results:\n")
print(main_results[, display_columns], row.names = FALSE)

cat("\nProvisional miR-375 result:\n")
print(provisional_results[, display_columns], row.names = FALSE)

cat("\nInterpretation:\n")
cat("Fold enrichment > 1 indicates more overlap than expected.\n")
cat("Main Enrichment_FDR < 0.05 supports statistical enrichment.\n")
cat("A non-significant result does not invalidate individual candidates.\n")
cat("This test does not establish direct regulation or causality.\n")
cat("Database coverage and gene-symbol mapping can influence results.\n")
cat("The provisional miR-375 result must not be treated as confirmed.\n")