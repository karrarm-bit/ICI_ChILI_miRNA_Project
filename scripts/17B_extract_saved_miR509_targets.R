# ==============================================================================
# Script 17B: Extract miRTarBase targets from the saved miR-509-3p query
# This script performs no online queries.
# ==============================================================================

# Load the package namespace to register the saved object's class.
if (!requireNamespace("multiMiR", quietly = TRUE)) {
  stop("The multiMiR package is required to read the saved query object.")
}

selected_miRNA <- "hsa-miR-509-3p"

checkpoint_file <- paste0(
  "results/R_objects/miRNA_targets/",
  "hsa_miR_509_3p_miRTarBase_query.rds"
)

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

# ------------------------------------------------------------------------------
# 1. Read the saved query
# ------------------------------------------------------------------------------

if (!file.exists(checkpoint_file)) {
  stop(
    "The checkpoint file was not found. Open the correct RStudio project."
  )
}

cat("Loading the saved miRTarBase query...\n")

saved_query <- readRDS(checkpoint_file)

# The diagnostic output confirmed that this object contains a data slot.
target_evidence <- as.data.frame(saved_query@data)

cat("Downloaded evidence records:", nrow(target_evidence), "\n")

required_columns <- c(
  "database",
  "mature_mirna_id",
  "target_symbol"
)

if (!all(required_columns %in% names(target_evidence))) {
  stop(
    "Required result columns are missing: ",
    paste(
      setdiff(required_columns, names(target_evidence)),
      collapse = ", "
    )
  )
}

# ------------------------------------------------------------------------------
# 2. Check database and miRNA identifiers
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
  stop("No exact-match records remained after checking identifiers.")
}

# Preserve all downloaded evidence fields.
rownames(target_evidence) <- NULL

# ------------------------------------------------------------------------------
# 3. Extract unique target genes
# ------------------------------------------------------------------------------

gene_symbols <- trimws(
  as.character(target_evidence$target_symbol)
)

valid_gene_symbols <- (
  !is.na(gene_symbols) &
    nzchar(gene_symbols)
)

unique_genes <- sort(unique(gene_symbols[valid_gene_symbols]))

if (length(unique_genes) == 0L) {
  stop("No valid target-gene symbols were found.")
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

# ------------------------------------------------------------------------------
# 4. Save the evidence and unique target-gene tables
# ------------------------------------------------------------------------------

evidence_rds_file <- file.path(
  rds_directory,
  "hsa_miR_509_3p_miRTarBase_evidence.rds"
)

targets_rds_file <- file.path(
  rds_directory,
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

saveRDS(target_evidence, evidence_rds_file)
saveRDS(target_gene_table, targets_rds_file)

write.csv(
  target_evidence,
  evidence_csv_file,
  row.names = FALSE,
  na = ""
)

write.csv(
  target_gene_table,
  targets_csv_file,
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 5. Completion report
# ------------------------------------------------------------------------------

cat("\n============================================\n")
cat("Saved target extraction completed\n")
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
cat(evidence_rds_file, "\n")
cat(targets_rds_file, "\n")
cat(evidence_csv_file, "\n")
cat(targets_csv_file, "\n")

cat("\nNo online query was performed.\n")