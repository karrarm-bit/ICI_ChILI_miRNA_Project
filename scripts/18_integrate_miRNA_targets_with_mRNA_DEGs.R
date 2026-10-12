# ==============================================================================
# Script 18: Exploratory integration of miRNA targets with mRNA DEGs
#
# mRNA comparison: ChILI_TOL versus Control_AC
# miRNA comparison: Responder versus Non_responder
#
# These are different cohorts and biological comparisons.
# Overlap does not demonstrate direct regulation or inverse correlation.
# miR-375 alias results remain provisional.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Define input files
# ------------------------------------------------------------------------------

mrna_file <- "results/R_objects/GSE287540_DE_results.rds"

target_directory <- "results/R_objects/miRNA_targets"

target_files <- c(
  file.path(
    target_directory,
    "hsa_miR_509_3p_miRTarBase_unique_genes.rds"
  ),
  file.path(
    target_directory,
    "hsa_miR_26b_5p_miRTarBase_unique_genes.rds"
  ),
  file.path(
    target_directory,
    "hsa_miR_196b_5p_miRTarBase_unique_genes.rds"
  )
)

provisional_file <- file.path(
  target_directory,
  "hsa_miR_375_3p_alias_miRTarBase_unique_genes.rds"
)

input_files <- c(mrna_file, target_files, provisional_file)

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0L) {
  stop(
    "Missing input files:\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 2. Create output directories
# ------------------------------------------------------------------------------

rds_directory <- "results/R_objects"
csv_directory <- "data/processed/target_integration"
table_directory <- "results/tables"

dir.create(rds_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(csv_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 3. Read and check mRNA results
# ------------------------------------------------------------------------------

cat("Loading mRNA differential-expression results...\n")

mrna_results <- readRDS(mrna_file)

required_mrna_columns <- c(
  "Gene_symbol",
  "Ensembl_ID",
  "Entrez_ID",
  "logFC",
  "PValue",
  "FDR"
)

if (!is.data.frame(mrna_results)) {
  stop("The mRNA results object must be a data frame.")
}

missing_columns <- setdiff(
  required_mrna_columns,
  names(mrna_results)
)

if (length(missing_columns) > 0L) {
  stop(
    "Missing mRNA columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (!is.numeric(mrna_results$logFC) ||
    !is.numeric(mrna_results$FDR)) {
  stop("mRNA logFC and FDR columns must be numeric.")
}

mrna_results$Match_symbol <- toupper(
  trimws(as.character(mrna_results$Gene_symbol))
)

# Reproduce the primary DEG thresholds.
significant_rows <- (
  is.finite(mrna_results$logFC) &
    is.finite(mrna_results$FDR) &
    mrna_results$FDR < 0.05 &
    abs(mrna_results$logFC) >= 1
)

mrna_degs <- mrna_results[
  significant_rows,
  ,
  drop = FALSE
]

cat("Primary DEG rows:", nrow(mrna_degs), "\n")

if (nrow(mrna_degs) != 203L) {
  warning(
    "The DEG count differs from the previously reported 203. ",
    "Review the input file and thresholds."
  )
}

mrna_degs$ChILI_direction <- ifelse(
  mrna_degs$logFC > 0,
  "Higher_in_ChILI",
  "Lower_in_ChILI"
)

# ------------------------------------------------------------------------------
# 4. Read miRNA target tables
# ------------------------------------------------------------------------------

read_target_table <- function(file, mapping_status) {
  
  target_table <- readRDS(file)
  
  required_columns <- c(
    "miRNA",
    "Target_symbol",
    "Database",
    "Evidence_records"
  )
  
  if (!is.data.frame(target_table) ||
      !all(required_columns %in% names(target_table))) {
    stop("Unexpected target table structure in: ", file)
  }
  
  target_table$Match_symbol <- toupper(
    trimws(as.character(target_table$Target_symbol))
  )
  
  valid_rows <- (
    !is.na(target_table$Match_symbol) &
      nzchar(target_table$Match_symbol)
  )
  
  target_table <- target_table[
    valid_rows,
    ,
    drop = FALSE
  ]
  
  target_table$Mapping_status <- mapping_status
  
  # A consistent structure allows tables to be combined safely.
  target_table[, c(
    "miRNA",
    "Target_symbol",
    "Match_symbol",
    "Database",
    "Evidence_records",
    "Mapping_status"
  )]
}

main_target_list <- lapply(
  target_files,
  read_target_table,
  mapping_status = "Exact_name_query"
)

main_targets <- do.call(rbind, main_target_list)
rownames(main_targets) <- NULL

provisional_targets <- read_target_table(
  provisional_file,
  mapping_status = "Alias_pending_verification"
)

all_targets <- rbind(
  main_targets,
  provisional_targets
)

# ------------------------------------------------------------------------------
# 5. Prepare mRNA information for matching
# ------------------------------------------------------------------------------

mrna_match_table <- data.frame(
  Match_symbol = mrna_degs$Match_symbol,
  mRNA_Gene_symbol = as.character(mrna_degs$Gene_symbol),
  mRNA_Ensembl_ID = as.character(mrna_degs$Ensembl_ID),
  mRNA_Entrez_ID = as.character(mrna_degs$Entrez_ID),
  mRNA_logFC = mrna_degs$logFC,
  mRNA_PValue = mrna_degs$PValue,
  mRNA_FDR = mrna_degs$FDR,
  ChILI_direction = mrna_degs$ChILI_direction,
  stringsAsFactors = FALSE
)

valid_mrna_symbols <- (
  !is.na(mrna_match_table$Match_symbol) &
    nzchar(mrna_match_table$Match_symbol)
)

mrna_match_table <- mrna_match_table[
  valid_mrna_symbols,
  ,
  drop = FALSE
]

# ------------------------------------------------------------------------------
# 6. Match target symbols to primary DEGs
# ------------------------------------------------------------------------------

match_targets <- function(target_table) {
  
  matched_table <- merge(
    target_table,
    mrna_match_table,
    by = "Match_symbol",
    all = FALSE,
    sort = FALSE
  )
  
  matched_table$miRNA_comparison <- rep(
    "Responder_vs_Non_responder",
    nrow(matched_table)
  )
  
  matched_table$mRNA_comparison <- rep(
    "ChILI_TOL_vs_Control_AC",
    nrow(matched_table)
  )
  
  matched_table$Interpretation <- rep(
    "Exploratory_cross_dataset_target_overlap",
    nrow(matched_table)
  )
  
  matched_table <- matched_table[
    order(
      matched_table$miRNA,
      matched_table$mRNA_FDR,
      matched_table$mRNA_Gene_symbol
    ),
    ,
    drop = FALSE
  ]
  
  rownames(matched_table) <- NULL
  
  matched_table
}

main_overlap <- match_targets(main_targets)

provisional_overlap <- match_targets(
  provisional_targets
)

main_overlap_up <- main_overlap[
  main_overlap$ChILI_direction == "Higher_in_ChILI",
  ,
  drop = FALSE
]

main_overlap_down <- main_overlap[
  main_overlap$ChILI_direction == "Lower_in_ChILI",
  ,
  drop = FALSE
]

# ------------------------------------------------------------------------------
# 7. Summarize overlap for each miRNA
# ------------------------------------------------------------------------------

all_overlap <- rbind(
  main_overlap,
  provisional_overlap
)

miRNA_names <- unique(all_targets$miRNA)

summary_list <- lapply(
  miRNA_names,
  function(current_miRNA) {
    
    current_targets <- all_targets[
      all_targets$miRNA == current_miRNA,
      ,
      drop = FALSE
    ]
    
    current_overlap <- all_overlap[
      all_overlap$miRNA == current_miRNA,
      ,
      drop = FALSE
    ]
    
    data.frame(
      miRNA = current_miRNA,
      Mapping_status = current_targets$Mapping_status[1],
      Unique_target_genes = length(
        unique(current_targets$Match_symbol)
      ),
      Overlapping_DEG_symbols = length(
        unique(current_overlap$Match_symbol)
      ),
      Higher_in_ChILI = length(
        unique(
          current_overlap$Match_symbol[
            current_overlap$ChILI_direction == "Higher_in_ChILI"
          ]
        )
      ),
      Lower_in_ChILI = length(
        unique(
          current_overlap$Match_symbol[
            current_overlap$ChILI_direction == "Lower_in_ChILI"
          ]
        )
      ),
      stringsAsFactors = FALSE
    )
  }
)

overlap_summary <- do.call(rbind, summary_list)
rownames(overlap_summary) <- NULL

# ------------------------------------------------------------------------------
# 8. Save results
# ------------------------------------------------------------------------------

output_tables <- list(
  Target_summary = overlap_summary,
  Main_target_genes = main_targets,
  Main_DEG_overlap = main_overlap,
  Main_ChILI_up = main_overlap_up,
  Main_ChILI_down = main_overlap_down,
  Provisional_375_targets = provisional_targets,
  Provisional_375_overlap = provisional_overlap
)

saveRDS(
  output_tables,
  file.path(
    rds_directory,
    "GSE274664_GSE287540_target_integration.rds"
  )
)

for (table_name in names(output_tables)) {
  
  write.csv(
    output_tables[[table_name]],
    file.path(
      csv_directory,
      paste0("Target_integration_", table_name, ".csv")
    ),
    row.names = FALSE,
    na = ""
  )
}

# Excel export is optional; RDS and CSV files are always saved.
if (requireNamespace("openxlsx", quietly = TRUE)) {
  
  openxlsx::write.xlsx(
    output_tables,
    file = file.path(
      table_directory,
      "GSE274664_GSE287540_target_integration.xlsx"
    ),
    overwrite = TRUE
  )
  
  cat("\nExcel workbook saved successfully.\n")
  
} else {
  
  cat("\nopenxlsx is unavailable; RDS and CSV outputs were saved.\n")
}

# ------------------------------------------------------------------------------
# 9. Completion report
# ------------------------------------------------------------------------------

cat("\n============================================\n")
cat("Exploratory target integration completed\n")
cat("============================================\n")

cat("Primary DEG rows:", nrow(mrna_degs), "\n")

cat(
  "Unique overlapping genes for the three exact-name miRNAs:",
  length(unique(main_overlap$Match_symbol)),
  "\n"
)

cat(
  "Provisional miR-375 overlapping genes:",
  length(unique(provisional_overlap$Match_symbol)),
  "\n"
)

cat("\nPer-miRNA overlap summary:\n")
print(overlap_summary, row.names = FALSE)

cat("\nMain overlapping genes:\n")

if (nrow(main_overlap) > 0L) {
  
  print(
    main_overlap[, c(
      "miRNA",
      "mRNA_Gene_symbol",
      "mRNA_logFC",
      "mRNA_FDR",
      "ChILI_direction"
    )],
    row.names = FALSE
  )
  
} else {
  
  cat("No exact-symbol overlap was found for the three main miRNAs.\n")
}

cat("\nImportant limitations:\n")
cat("Matching used gene symbols, not verified gene aliases.\n")
cat("No overlap-enrichment significance test was performed.\n")
cat("miR-375 alias results remain separate and provisional.\n")
cat("The two datasets represent different biological comparisons.\n")
cat("Overlap does not establish direct regulation or causality.\n")

cat("\nOutput directory:", csv_directory, "\n")