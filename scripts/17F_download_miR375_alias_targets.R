# ==============================================================================
# Script 17F: Query miRTarBase using an alternative miR-375 identifier
# The previous zero-result checkpoint is preserved.
# ==============================================================================

options(timeout = 600)

if (!requireNamespace("multiMiR", quietly = TRUE)) {
  stop("The multiMiR package is required.", call. = FALSE)
}

if (!dir.exists("data") || !dir.exists("results")) {
  stop("Open the correct RStudio project.", call. = FALSE)
}

# Keep the study identifier separate from the database query identifier.
study_miRNA <- "hsa-miR-375-3p"
query_identifier <- "hsa-miR-375"

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

# A different checkpoint prevents reuse of the previous empty query.
checkpoint_file <- file.path(
  rds_directory,
  "hsa_miR_375_alias_miRTarBase_query.rds"
)

cat("Study miRNA:", study_miRNA, "\n")
cat("Alternative query identifier:", query_identifier, "\n")

# ------------------------------------------------------------------------------
# 1. Query or load the alternative-identifier checkpoint
# ------------------------------------------------------------------------------

if (file.exists(checkpoint_file)) {
  
  cat("\nLoading the alternative-identifier checkpoint...\n")
  target_query <- readRDS(checkpoint_file)
  
} else {
  
  cat("\nQuerying miRTarBase using the alternative identifier...\n")
  
  target_query <- tryCatch(
    {
      multiMiR::get_multimir(
        org = "hsa",
        mirna = query_identifier,
        table = "mirtarbase",
        summary = FALSE,
        add.link = FALSE,
        use.tibble = FALSE,
        legacy.out = FALSE
      )
    },
    error = function(e) {
      stop(
        "Query failed:\n",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  saveRDS(target_query, checkpoint_file)
  cat("Alternative query checkpoint saved.\n")
}

# ------------------------------------------------------------------------------
# 2. Extract and inspect returned identifiers
# ------------------------------------------------------------------------------

target_evidence <- as.data.frame(target_query@data)

cat("Returned evidence records:", nrow(target_evidence), "\n")

if (nrow(target_evidence) == 0L) {
  stop(
    paste0(
      "The alternative identifier also returned no records.\n",
      "Send this output before trying another identifier or database."
    ),
    call. = FALSE
  )
}

required_columns <- c(
  "database",
  "mature_mirna_id",
  "target_symbol"
)

if (!all(required_columns %in% names(target_evidence))) {
  stop(
    "Missing columns: ",
    paste(
      setdiff(required_columns, names(target_evidence)),
      collapse = ", "
    ),
    call. = FALSE
  )
}

cat("\nReturned database miRNA names:\n")
print(unique(target_evidence$mature_mirna_id))

if ("mature_mirna_acc" %in% names(target_evidence)) {
  cat("\nReturned mature miRNA accession numbers:\n")
  print(unique(target_evidence$mature_mirna_acc))
}

# ------------------------------------------------------------------------------
# 3. Keep only exact miR-375 aliases from miRTarBase
# ------------------------------------------------------------------------------

database_names <- tolower(
  trimws(as.character(target_evidence$database))
)

miRNA_names <- tolower(
  trimws(as.character(target_evidence$mature_mirna_id))
)

accepted_names <- c(
  "hsa-mir-375",
  "hsa-mir-375-3p"
)

keep_records <- (
  !is.na(database_names) &
    !is.na(miRNA_names) &
    database_names == "mirtarbase" &
    miRNA_names %in% accepted_names
)

target_evidence <- target_evidence[
  keep_records,
  ,
  drop = FALSE
]

if (nrow(target_evidence) == 0L) {
  stop(
    "No exact accepted miR-375 names remained. Review returned identifiers.",
    call. = FALSE
  )
}

# Preserve the database name; do not replace the original evidence identifier.
target_evidence$Study_miRNA <- study_miRNA
target_evidence$Query_identifier <- query_identifier
target_evidence$Identifier_mapping_status <- "Pending_accession_verification"

rownames(target_evidence) <- NULL

# ------------------------------------------------------------------------------
# 4. Extract unique genes
# ------------------------------------------------------------------------------

gene_symbols <- trimws(
  as.character(target_evidence$target_symbol)
)

valid_symbols <- !is.na(gene_symbols) & nzchar(gene_symbols)

unique_genes <- sort(unique(gene_symbols[valid_symbols]))

if (length(unique_genes) == 0L) {
  stop("No valid gene symbols were found.", call. = FALSE)
}

database_accessions <- NA_character_

if ("mature_mirna_acc" %in% names(target_evidence)) {
  accession_values <- as.character(target_evidence$mature_mirna_acc)
  accession_values <- accession_values[
    !is.na(accession_values) & nzchar(accession_values)
  ]
  
  if (length(accession_values) > 0L) {
    database_accessions <- paste(
      sort(unique(accession_values)),
      collapse = "; "
    )
  }
}

target_gene_table <- data.frame(
  miRNA = rep(study_miRNA, length(unique_genes)),
  Target_symbol = unique_genes,
  Database = rep("miRTarBase", length(unique_genes)),
  Evidence_records = vapply(
    unique_genes,
    function(gene) {
      as.integer(sum(gene_symbols == gene, na.rm = TRUE))
    },
    integer(1)
  ),
  Database_miRNA = paste(
    sort(unique(as.character(target_evidence$mature_mirna_id))),
    collapse = "; "
  ),
  Database_accession = database_accessions,
  Identifier_mapping_status = "Pending_accession_verification",
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# 5. Save provisional alias results separately
# ------------------------------------------------------------------------------

evidence_rds_file <- file.path(
  rds_directory,
  "hsa_miR_375_3p_alias_miRTarBase_evidence.rds"
)

targets_rds_file <- file.path(
  rds_directory,
  "hsa_miR_375_3p_alias_miRTarBase_unique_genes.rds"
)

evidence_csv_file <- file.path(
  csv_directory,
  "hsa_miR_375_3p_alias_miRTarBase_evidence.csv"
)

targets_csv_file <- file.path(
  csv_directory,
  "hsa_miR_375_3p_alias_miRTarBase_unique_genes.csv"
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
  row.names = FALSE,
  na = ""
)

cat("\n============================================\n")
cat("Alternative-identifier target extraction completed\n")
cat("============================================\n")

cat("Study miRNA:", study_miRNA, "\n")
cat("Query identifier:", query_identifier, "\n")
cat("Evidence records:", nrow(target_evidence), "\n")
cat("Unique target genes:", nrow(target_gene_table), "\n")

cat("\nFirst 20 target genes:\n")
print(utils::head(target_gene_table, 20), row.names = FALSE)

cat("\nFiles created:\n")
cat(checkpoint_file, "\n")
cat(evidence_rds_file, "\n")
cat(targets_rds_file, "\n")
cat(evidence_csv_file, "\n")
cat(targets_csv_file, "\n")

cat("\nImportant:\n")
cat("Alias results are provisional until accession identity is verified.\n")
cat("Do not integrate these targets with mRNA results yet.\n")
mrna_results <- readRDS(
  "results/R_objects/GSE287540_DE_results.rds"
)

cat("\nOBJECT CLASS:\n")
print(class(mrna_results))

cat("\nCOLUMN NAMES:\n")
print(names(mrna_results))

if (is.data.frame(mrna_results)) {
  cat("\nFIRST THREE ROWS:\n")
  print(head(mrna_results, 3))
}