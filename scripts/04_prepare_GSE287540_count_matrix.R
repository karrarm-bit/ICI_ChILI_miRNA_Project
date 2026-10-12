# ============================================================
# Script 04: Prepare GSE287540 Raw Count Matrix
# Project: ICI-Induced Liver Injury and miRNA Analysis
# ============================================================

# Clear the R environment
rm(list = ls())

# Load required packages
suppressPackageStartupMessages({
  library(data.table)
  library(GEOquery)
  library(openxlsx)
})

# ============================================================
# Define project directories
# ============================================================

project_directory <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

raw_data_directory <- file.path(
  project_directory,
  "data",
  "raw"
)

processed_data_directory <- file.path(
  project_directory,
  "data",
  "processed"
)

tables_directory <- file.path(
  project_directory,
  "results",
  "tables"
)

objects_directory <- file.path(
  project_directory,
  "results",
  "R_objects"
)

# ============================================================
# Define input files
# ============================================================

counts_file <- file.path(
  raw_data_directory,
  "GSE287540_ChILI_bulk_readcount.tsv.gz"
)

geo_metadata_file <- file.path(
  raw_data_directory,
  "GSE287540_series_matrix.txt.gz"
)

clean_metadata_file <- file.path(
  processed_data_directory,
  "GSE287540_metadata_clean.csv"
)

# Check the required input files
required_files <- c(
  counts_file,
  geo_metadata_file,
  clean_metadata_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  cat("\nMissing files:\n")
  print(missing_files)
  
  stop(
    "Required input files are missing. ",
    "Run the previous scripts first."
  )
}

# ============================================================
# Read the clean metadata
# ============================================================

cat("\nReading the cleaned metadata...\n")

clean_metadata <- read.csv(
  clean_metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(
  "Metadata samples:",
  nrow(clean_metadata),
  "\n"
)

# ============================================================
# Extract library identifiers from the GEO metadata
# ============================================================

cat("\nExtracting sequencing library identifiers...\n")

gse_object <- GEOquery::getGEO(
  filename = geo_metadata_file
)

raw_geo_metadata <- Biobase::pData(
  gse_object
)

# Check whether the description column exists
if (!"description" %in% colnames(raw_geo_metadata)) {
  stop(
    "The description column was not found in the GEO metadata."
  )
}

library_id <- sub(
  pattern = "^Library name:[[:space:]]*",
  replacement = "",
  x = as.character(raw_geo_metadata$description),
  ignore.case = TRUE
)

library_map <- data.frame(
  Sample_ID = rownames(raw_geo_metadata),
  Library_ID = library_id,
  stringsAsFactors = FALSE
)

# Check for missing library identifiers
if (any(is.na(library_map$Library_ID))) {
  stop("Missing sequencing library identifiers were detected.")
}

if (any(library_map$Library_ID == "")) {
  stop("Empty sequencing library identifiers were detected.")
}

# Add library identifiers to the clean metadata
library_match <- match(
  clean_metadata$Sample_ID,
  library_map$Sample_ID
)

if (any(is.na(library_match))) {
  stop(
    "Some GEO sample IDs could not be matched ",
    "to sequencing library identifiers."
  )
}

clean_metadata$Library_ID <- library_map$Library_ID[
  library_match
]

# ============================================================
# Read the raw count matrix
# ============================================================

cat("\nReading the raw gene-count matrix...\n")

raw_counts <- data.table::fread(
  counts_file,
  header = TRUE,
  check.names = FALSE,
  data.table = FALSE
)

cat(
  "Raw count matrix dimensions:",
  nrow(raw_counts),
  "genes and",
  ncol(raw_counts) - 1,
  "samples\n"
)

# The first column contains gene identifiers
gene_ids <- as.character(
  raw_counts[[1]]
)

# Remove the gene identifier column
count_matrix <- as.matrix(
  raw_counts[, -1, drop = FALSE]
)

# Convert all expression values to numeric
storage.mode(count_matrix) <- "numeric"

# Assign gene identifiers as row names
rownames(count_matrix) <- gene_ids

# ============================================================
# Validate the raw count matrix
# ============================================================

if (any(is.na(count_matrix))) {
  stop("Missing values were detected in the count matrix.")
}

if (any(count_matrix < 0)) {
  stop("Negative values were detected in the count matrix.")
}

if (anyDuplicated(colnames(count_matrix)) > 0) {
  stop("Duplicated sequencing library identifiers were detected.")
}

if (anyDuplicated(rownames(count_matrix)) > 0) {
  
  cat(
    "\nDuplicated gene identifiers were detected.",
    "\nCollapsing duplicated genes by summing their counts...\n"
  )
  
  count_matrix <- rowsum(
    count_matrix,
    group = rownames(count_matrix),
    reorder = FALSE
  )
}

# ============================================================
# Match count columns with sample metadata
# ============================================================

cat("\nMatching count columns with sample metadata...\n")

count_column_match <- match(
  clean_metadata$Library_ID,
  colnames(count_matrix)
)

unmatched_metadata <- clean_metadata$Library_ID[
  is.na(count_column_match)
]

if (length(unmatched_metadata) > 0) {
  
  cat("\nUnmatched library identifiers:\n")
  print(unmatched_metadata)
  
  stop(
    "Some metadata samples could not be matched ",
    "with count-matrix columns."
  )
}

# Reorder count columns according to the metadata
count_matrix <- count_matrix[
  ,
  count_column_match,
  drop = FALSE
]

# Preserve the sequencing library names
original_library_ids <- colnames(count_matrix)

# Replace library names with GEO sample IDs
colnames(count_matrix) <- clean_metadata$Sample_ID

# Confirm that the order is identical
if (!identical(
  colnames(count_matrix),
  clean_metadata$Sample_ID
)) {
  stop(
    "The count matrix and metadata are not in the same order."
  )
}

# ============================================================
# Remove invalid and completely unexpressed genes
# ============================================================

valid_gene_id <- !is.na(rownames(count_matrix)) &
  rownames(count_matrix) != ""

count_matrix <- count_matrix[
  valid_gene_id,
  ,
  drop = FALSE
]

genes_before_zero_filter <- nrow(count_matrix)

nonzero_genes <- rowSums(count_matrix) > 0

count_matrix <- count_matrix[
  nonzero_genes,
  ,
  drop = FALSE
]

genes_removed_as_zero <- genes_before_zero_filter -
  nrow(count_matrix)

cat(
  "\nCompletely unexpressed genes removed:",
  genes_removed_as_zero,
  "\n"
)

# ============================================================
# Create the final sample map
# ============================================================

sample_map <- data.frame(
  Matrix_column = seq_len(ncol(count_matrix)),
  Sample_ID = clean_metadata$Sample_ID,
  Library_ID = original_library_ids,
  Sample_title = clean_metadata$Sample_title,
  Subject_ID = clean_metadata$Subject_ID,
  Group = clean_metadata$Group,
  Time_point = clean_metadata$Time_point,
  stringsAsFactors = FALSE
)

# ============================================================
# Calculate basic count statistics
# ============================================================

sample_statistics <- data.frame(
  Sample_ID = colnames(count_matrix),
  Library_ID = original_library_ids,
  Subject_ID = clean_metadata$Subject_ID,
  Group = clean_metadata$Group,
  Time_point = clean_metadata$Time_point,
  
  Total_reads = colSums(
    count_matrix
  ),
  
  Expressed_genes = colSums(
    count_matrix > 0
  ),
  
  stringsAsFactors = FALSE
)

sample_statistics$Library_size_million <- round(
  sample_statistics$Total_reads / 1000000,
  3
)

# Create a dataset summary
dataset_summary <- data.frame(
  Measurement = c(
    "Number of genes before zero-count filtering",
    "Number of zero-count genes removed",
    "Number of genes after zero-count filtering",
    "Number of samples",
    "Number of unique subjects",
    "Number of ChILI samples",
    "Number of control samples"
  ),
  
  Value = c(
    genes_before_zero_filter,
    genes_removed_as_zero,
    nrow(count_matrix),
    ncol(count_matrix),
    length(unique(clean_metadata$Subject_ID)),
    sum(clean_metadata$Group == "chili"),
    sum(clean_metadata$Group == "control")
  ),
  
  stringsAsFactors = FALSE
)

# ============================================================
# Save the prepared data
# ============================================================

saveRDS(
  count_matrix,
  file = file.path(
    objects_directory,
    "GSE287540_raw_count_matrix.rds"
  )
)

saveRDS(
  clean_metadata,
  file = file.path(
    objects_directory,
    "GSE287540_metadata_analysis_ready.rds"
  )
)

write.csv(
  clean_metadata,
  file = file.path(
    processed_data_directory,
    "GSE287540_metadata_analysis_ready.csv"
  ),
  row.names = FALSE
)

write.csv(
  sample_map,
  file = file.path(
    processed_data_directory,
    "GSE287540_sample_map.csv"
  ),
  row.names = FALSE
)

write.csv(
  sample_statistics,
  file = file.path(
    processed_data_directory,
    "GSE287540_sample_count_statistics.csv"
  ),
  row.names = FALSE
)

# ============================================================
# Create an Excel quality-control workbook
# ============================================================

count_workbook <- createWorkbook()

addWorksheet(
  count_workbook,
  "Dataset summary"
)

writeData(
  count_workbook,
  sheet = "Dataset summary",
  x = dataset_summary
)

addWorksheet(
  count_workbook,
  "Sample map"
)

writeData(
  count_workbook,
  sheet = "Sample map",
  x = sample_map
)

addWorksheet(
  count_workbook,
  "Sample statistics"
)

writeData(
  count_workbook,
  sheet = "Sample statistics",
  x = sample_statistics
)

for (sheet_name in names(count_workbook)) {
  
  setColWidths(
    count_workbook,
    sheet = sheet_name,
    cols = 1:10,
    widths = "auto"
  )
}

saveWorkbook(
  count_workbook,
  file = file.path(
    tables_directory,
    "GSE287540_count_matrix_preparation.xlsx"
  ),
  overwrite = TRUE
)

# ============================================================
# Final validation
# ============================================================

expected_sample_number <- 128

if (ncol(count_matrix) != expected_sample_number) {
  stop(
    "The final count matrix does not contain ",
    expected_sample_number,
    " samples."
  )
}

if (nrow(clean_metadata) != ncol(count_matrix)) {
  stop(
    "The number of metadata rows does not match ",
    "the number of count-matrix columns."
  )
}

if (!identical(
  clean_metadata$Sample_ID,
  colnames(count_matrix)
)) {
  stop(
    "The metadata and count matrix are not aligned."
  )
}

cat("\n============================================\n")
cat("GSE287540 count matrix preparation completed\n")
cat("============================================\n")

cat(
  "\nFinal count matrix:",
  nrow(count_matrix),
  "genes and",
  ncol(count_matrix),
  "samples\n"
)

cat("\nFiles created:\n")
cat("results/R_objects/GSE287540_raw_count_matrix.rds\n")
cat("results/R_objects/GSE287540_metadata_analysis_ready.rds\n")
cat("data/processed/GSE287540_metadata_analysis_ready.csv\n")
cat("data/processed/GSE287540_sample_map.csv\n")
cat("data/processed/GSE287540_sample_count_statistics.csv\n")
cat("results/tables/GSE287540_count_matrix_preparation.xlsx\n")