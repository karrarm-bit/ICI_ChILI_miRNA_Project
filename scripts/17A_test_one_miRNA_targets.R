# ==============================================================================
# Script 17A: Install required packages and retrieve targets for ONE human miRNA
# miRNA: hsa-miR-509-3p
# Database: miRTarBase through multiMiR
# ==============================================================================

# Open the ICI_ChILI_miRNA_Project RStudio project before running this script.
# Installation requires an Internet connection.
# Previous expression data and analysis results are not modified.

options(
  timeout = 600,
  repos = c(CRAN = "https://cloud.r-project.org")
)

# ------------------------------------------------------------------------------
# 1. Install BiocManager if necessary
# ------------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  
  cat("Installing BiocManager...\n")
  
  tryCatch(
    install.packages("BiocManager"),
    error = function(e) {
      stop(
        "BiocManager installation failed:\n",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  stop(
    "BiocManager could not be loaded. Review the installation messages.",
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 2. Install multiMiR if necessary
# ------------------------------------------------------------------------------

if (!requireNamespace("multiMiR", quietly = TRUE)) {
  
  cat("\nInstalling multiMiR and its required dependencies...\n")
  cat("Package installation may take several minutes.\n")
  
  tryCatch(
    BiocManager::install(
      "multiMiR",
      ask = FALSE,
      update = FALSE
    ),
    error = function(e) {
      stop(
        "multiMiR installation failed:\n",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
}

# Confirm successful installation and namespace loading.
if (!requireNamespace("multiMiR", quietly = TRUE)) {
  stop(
    paste0(
      "multiMiR is still unavailable after the installation attempt.\n",
      "Do not continue to target extraction.\n",
      "Send the installation messages shown above this error."
    ),
    call. = FALSE
  )
}

cat(
  "\nmultiMiR is available. Package version:",
  as.character(utils::packageVersion("multiMiR")),
  "\n"
)

# ------------------------------------------------------------------------------
# 3. Check the project location
# ------------------------------------------------------------------------------

if (!dir.exists("data") || !dir.exists("results")) {
  stop(
    paste0(
      "The project folders 'data' and 'results' were not found.\n",
      "Open the correct RStudio project before running this script.\n",
      "Current working directory: ",
      getwd()
    ),
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 4. Define the selected miRNA and output files
# ------------------------------------------------------------------------------

selected_miRNA <- "hsa-miR-509-3p"

output_directory <- "results/R_objects/miRNA_targets"
csv_directory <- "data/processed/miRNA_targets"

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  csv_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

checkpoint_file <- file.path(
  output_directory,
  "hsa_miR_509_3p_miRTarBase_query.rds"
)

evidence_rds_file <- file.path(
  output_directory,
  "hsa_miR_509_3p_miRTarBase_evidence.rds"
)

targets_rds_file <- file.path(
  output_directory,
  "hsa_miR_509_3p_miRTarBase_unique_genes.rds"
)

evidence_csv_file <- file.path(
  csv_directory,
  "hsa_miR_509_3p_miRTarBase_evidence.csv"
)

targets_csv_file <- file.path(
  csv_directory,
  "hsa_miR_509_3p_miRTarBase_unique_genes.csv"
)

# ------------------------------------------------------------------------------
# 5. Download or load the saved query checkpoint
# ------------------------------------------------------------------------------

cat("\n============================================\n")
cat("Single-miRNA target extraction\n")
cat("============================================\n")
cat("Selected miRNA:", selected_miRNA, "\n")
cat("Database: miRTarBase only\n")

if (file.exists(checkpoint_file)) {
  
  cat("\nLoading the saved query checkpoint...\n")
  
  target_query <- readRDS(checkpoint_file)
  
} else {
  
  cat("\nStarting one online query...\n")
  cat("Do not run this script again while the query is active.\n")
  cat("If no progress occurs for five minutes, stop the query once.\n")
  
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
        paste0(
          "The online target query failed.\n",
          "Original message: ",
          conditionMessage(e),
          "\nNo expression data were modified."
        ),
        call. = FALSE
      )
    }
  )
  
  # Save the downloaded query immediately.
  saveRDS(target_query, checkpoint_file)
  
  cat("\nQuery checkpoint saved successfully.\n")
}

# ------------------------------------------------------------------------------
# 6. Extract and check the evidence table
# ------------------------------------------------------------------------------

if (
  !isS4(target_query) ||
  !"data" %in% methods::slotNames(target_query)
) {
  stop(
    "The query object has an unexpected structure.",
    call. = FALSE
  )
}

target_evidence <- as.data.frame(
  methods::slot(target_query, "data")
)

cat("Returned evidence records:", nrow(target_evidence), "\n")

if (nrow(target_evidence) == 0L) {
  stop(
    paste0(
      "No records were returned for this miRNA.\n",
      "This does not prove that the miRNA has no targets.\n",
      "Its identifier and database coverage must be checked."
    ),
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

# Keep records belonging to the selected database and exact miRNA.
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
    "No exact-match miRTarBase records remained after checking identifiers.",
    call. = FALSE
  )
}

target_evidence <- unique(target_evidence)
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
# 7. Extract unique target-gene symbols
# ------------------------------------------------------------------------------

gene_symbols <- trimws(
  as.character(target_evidence$target_symbol)
)

gene_symbols <- gene_symbols[
  !is.na(gene_symbols) & nzchar(gene_symbols)
]

gene_symbols <- sort(unique(gene_symbols))

if (length(gene_symbols) == 0L) {
  stop(
    "Evidence records were saved, but no valid target-gene symbols were found.",
    call. = FALSE
  )
}

unique_target_genes <- data.frame(
  miRNA = rep(selected_miRNA, length(gene_symbols)),
  Target_symbol = gene_symbols,
  Database = rep("miRTarBase", length(gene_symbols)),
  stringsAsFactors = FALSE
)

saveRDS(
  unique_target_genes,
  targets_rds_file
)

write.csv(
  unique_target_genes,
  targets_csv_file,
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 8. Print the completion report
# ------------------------------------------------------------------------------

cat("\n============================================\n")
cat("Target extraction completed successfully\n")
cat("============================================\n")

cat("miRNA:", selected_miRNA, "\n")
cat("Evidence records:", nrow(target_evidence), "\n")
cat("Unique target genes:", nrow(unique_target_genes), "\n")

cat("\nFirst 20 target genes:\n")

print(
  utils::head(unique_target_genes, 20),
  row.names = FALSE
)

cat("\nFiles created:\n")
cat(checkpoint_file, "\n")
cat(evidence_rds_file, "\n")
cat(targets_rds_file, "\n")
cat(evidence_csv_file, "\n")
cat(targets_csv_file, "\n")

cat("\nNext step:\n")
cat("Send this completion report before querying the second miRNA.\n")
# Inspect the saved query without modifying it.

checkpoint_object <- readRDS(
  "results/R_objects/miRNA_targets/hsa_miR_509_3p_miRTarBase_query.rds"
)

cat("\nOBJECT CLASS:\n")
print(class(checkpoint_object))

cat("\nIS S4 OBJECT:\n")
print(isS4(checkpoint_object))

cat("\nOBJECT STRUCTURE:\n")
str(checkpoint_object, max.level = 2)

if (isS4(checkpoint_object)) {
  cat("\nAVAILABLE SLOTS:\n")
  print(methods::slotNames(checkpoint_object))
} else {
  cat("\nAVAILABLE NAMES:\n")
  print(names(checkpoint_object))
}