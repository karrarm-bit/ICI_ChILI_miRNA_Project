# ==============================================================================
# Script 17D: Retrieve miRTarBase targets for hsa-miR-196b-5p
# Run from the root of the ICI_ChILI_miRNA_Project RStudio project.
# ==============================================================================

options(timeout = 600)

# ------------------------------------------------------------------------------
# 1. Check the required package and project folders
# ------------------------------------------------------------------------------

if (!requireNamespace("multiMiR", quietly = TRUE)) {
  stop(
    "multiMiR is unavailable. Run BiocManager::install('multiMiR') first.",
    call. = FALSE
  )
}

if (!dir.exists("data") || !dir.exists("results")) {
  stop(
    "Open the correct RStudio project before running this script.",
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 2. Define the miRNA and output locations
# ------------------------------------------------------------------------------

selected_miRNA <- "hsa-miR-196b-5p"
file_prefix <- "hsa_miR_196b_5p"

rds_directory <- "results/R_objects/miRNA_targets"
csv_directory <- "data/processed/miRNA_targets"

dir.create(
  rds_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  csv_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

checkpoint_file <- file.path(
  rds_directory,
  paste0(file_prefix, "_miRTarBase_query.rds")
)

evidence_rds_file <- file.path(
  rds_directory,
  paste0(file_prefix, "_miRTarBase_evidence.rds")
)

targets_rds_file <- file.path(
  rds_directory,
  paste0(file_prefix, "_miRTarBase_unique_genes.rds")
)

evidence_csv_file <- file.path(
  csv_directory,
  paste0(file_prefix, "_miRTarBase_evidence.csv")
)

targets_csv_file <- file.path(
  csv_directory,
  paste0(file_prefix, "_miRTarBase_unique_genes.csv")
)

# ------------------------------------------------------------------------------
# 3. Retrieve targets or load the saved checkpoint
# ------------------------------------------------------------------------------

cat("============================================\n")
cat("Single-miRNA miRTarBase target extraction\n")
cat("============================================\n")
cat("Selected miRNA:", selected_miRNA, "\n")

if (file.exists(checkpoint_file)) {
  
  cat("\nLoading the saved query checkpoint...\n")
  
  target_query <- readRDS(checkpoint_file)
  
} else {
  
  cat("\nSearching miRTarBase for one miRNA only...\n")
  cat("Do not run the script again while the query is active.\n")
  
  target_query <- tryCatch(
    {
      multiMiR::get_multimir(
        org = "hsa",
        mirna = selected_miRNA,
        table = "mirtarbase",
        summary = FALSE,
        add.link = FALSE,
        use.tibble = FALSE,
        legacy.out = FALSE
      )
    },
    error = function(e) {
      stop(
        "The online query failed:\n",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  # Save before processing to avoid repeating a successful download.
  saveRDS(target_query, checkpoint_file)
  
  cat("\nQuery checkpoint saved successfully.\n")
}

# ------------------------------------------------------------------------------
# 4. Extract the evidence table
# ------------------------------------------------------------------------------

target_evidence <- tryCatch(
  as.data.frame(target_query@data),
  error = function(e) {
    stop(
      "The saved query could not be read:\n",
      conditionMessage(e),
      call. = FALSE
    )
  }
)

cat("Downloaded evidence records:", nrow(target_evidence), "\n")

if (nrow(target_evidence) == 0L) {
  stop(
    "No records were returned. Check identifier and database coverage.",
    call. = FALSE
  )
}

required_columns <- c(
  "database",
  "mature_mirna_id",
  "target_symbol"
)

missing_columns <- setdiff(
  required_columns,
  names(target_evidence)
)

if (length(missing_columns) > 0L) {
  stop(
    "Missing result columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 5. Keep exact-match miRTarBase evidence
# ------------------------------------------------------------------------------

database_values <- tolower(
  trimws(as.character(target_evidence$database))
)

miRNA_values <- tolower(
  trimws(as.character(target_evidence$mature_mirna_id))
)

keep_records <- (
  !is.na(database_values) &
    !is.na(miRNA_values) &
    database_values == "mirtarbase" &
    miRNA_values == tolower(selected_miRNA)
)

target_evidence <- target_evidence[
  keep_records,
  ,
  drop = FALSE
]

if (nrow(target_evidence) == 0L) {
  stop(
    "No exact-match miRTarBase records remained.",
    call. = FALSE
  )
}

rownames(target_evidence) <- NULL

# Preserve all evidence fields, including experiments and publication IDs.
saveRDS(target_evidence, evidence_rds_file)

write.csv(
  target_evidence,
  evidence_csv_file,
  row.names = FALSE,
  na = ""
)

# ------------------------------------------------------------------------------
# 6. Extract unique target-gene symbols
# ------------------------------------------------------------------------------

gene_symbols <- trimws(
  as.character(target_evidence$target_symbol)
)

valid_symbols <- (
  !is.na(gene_symbols) &
    nzchar(gene_symbols)
)

unique_genes <- sort(
  unique(gene_symbols[valid_symbols])
)

if (length(unique_genes) == 0L) {
  stop(
    "Evidence was saved, but no valid gene symbols were found.",
    call. = FALSE
  )
}

target_gene_table <- data.frame(
  miRNA = rep(selected_miRNA, length(unique_genes)),
  Target_symbol = unique_genes,
  Database = rep("miRTarBase", length(unique_genes)),
  Evidence_records = vapply(
    unique_genes,
    function(gene) {
      as.integer(sum(gene_symbols == gene, na.rm = TRUE))
    },
    integer(1)
  ),
  stringsAsFactors = FALSE
)

saveRDS(
  target_gene_table,
  targets_rds_file
)

write.csv(
  target_gene_table,
  targets_csv_file,
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 7. Completion report
# ------------------------------------------------------------------------------

cat("\n============================================\n")
cat("miRTarBase target extraction completed\n")
cat("============================================\n")

cat("miRNA:", selected_miRNA, "\n")
cat("Evidence records:", nrow(target_evidence), "\n")
cat("Unique target genes:", nrow(target_gene_table), "\n")

cat("\nFirst 20 target genes:\n")

print(
  utils::head(target_gene_table, 20),
  row.names = FALSE
)

cat("\nFiles created:\n")
cat(checkpoint_file, "\n")
cat(evidence_rds_file, "\n")
cat(targets_rds_file, "\n")
cat(evidence_csv_file, "\n")
cat(targets_csv_file, "\n")

cat("\nPrevious miRNA target files were not modified.\n")