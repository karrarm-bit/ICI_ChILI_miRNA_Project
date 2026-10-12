# ============================================================
# Script 06: Review QC and Define the Primary Analysis Samples
# Project: ICI-Induced Liver Injury and miRNA Analysis
# ============================================================

# Clear the R environment
rm(list = ls())

# Load required packages
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(openxlsx)
  library(patchwork)
})

# ============================================================
# Define project directories
# ============================================================

project_directory <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

objects_directory <- file.path(
  project_directory,
  "results",
  "R_objects"
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

figures_directory <- file.path(
  project_directory,
  "results",
  "figures"
)

# ============================================================
# Define input files
# ============================================================

dge_file <- file.path(
  objects_directory,
  "GSE287540_DGEList_filtered_TMM.rds"
)

log_cpm_file <- file.path(
  objects_directory,
  "GSE287540_logCPM_normalized.rds"
)

metadata_file <- file.path(
  objects_directory,
  "GSE287540_metadata_analysis_ready.rds"
)

qc_file <- file.path(
  processed_data_directory,
  "GSE287540_sample_QC_table.csv"
)

# Check required files
required_files <- c(
  dge_file,
  log_cpm_file,
  metadata_file,
  qc_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  cat("\nMissing files:\n")
  print(missing_files)
  
  stop(
    "Required files are missing. ",
    "Run 05_GSE287540_quality_control.R first."
  )
}

# ============================================================
# Load data
# ============================================================

cat("\nLoading quality-control results...\n")

dge_object <- readRDS(
  dge_file
)

log_cpm <- readRDS(
  log_cpm_file
)

metadata <- readRDS(
  metadata_file
)

qc_table <- read.csv(
  qc_file,
  stringsAsFactors = FALSE
)

# Confirm sample alignment
if (!identical(
  colnames(dge_object),
  metadata$Sample_ID
)) {
  stop(
    "The DGEList and metadata samples are not aligned."
  )
}

if (!identical(
  colnames(log_cpm),
  metadata$Sample_ID
)) {
  stop(
    "The logCPM matrix and metadata samples are not aligned."
  )
}

# ============================================================
# Prepare metadata variables
# ============================================================

metadata$Group <- factor(
  metadata$Group,
  levels = c("control", "chili")
)

metadata$Time_point <- factor(
  metadata$Time_point,
  levels = c("BC", "AC", "TOL", "FU1", "FU2")
)

metadata$Batch <- factor(
  metadata$Batch
)

metadata$Sex <- factor(
  metadata$Sex
)

metadata$Treatment_class <- factor(
  metadata$Treatment_class
)

# Add QC status
qc_match <- match(
  metadata$Sample_ID,
  qc_table$Sample_ID
)

if (any(is.na(qc_match))) {
  stop(
    "Some samples could not be matched with the QC table."
  )
}

metadata$QC_status <- qc_table$QC_status[
  qc_match
]

metadata$Potential_QC_outlier <-
  qc_table$Potential_QC_outlier[
    qc_match
  ]

# ============================================================
# Review study structure
# ============================================================

group_time_table <- as.data.frame(
  table(
    metadata$Group,
    metadata$Time_point
  )
)

colnames(group_time_table) <- c(
  "Group",
  "Time_point",
  "Number_of_samples"
)

batch_group_table <- as.data.frame(
  table(
    metadata$Batch,
    metadata$Group
  )
)

colnames(batch_group_table) <- c(
  "Batch",
  "Group",
  "Number_of_samples"
)

batch_time_table <- as.data.frame(
  table(
    metadata$Batch,
    metadata$Time_point
  )
)

colnames(batch_time_table) <- c(
  "Batch",
  "Time_point",
  "Number_of_samples"
)

batch_group_time_table <- as.data.frame(
  table(
    metadata$Batch,
    metadata$Group,
    metadata$Time_point
  )
)

colnames(batch_group_time_table) <- c(
  "Batch",
  "Group",
  "Time_point",
  "Number_of_samples"
)

# Remove combinations containing zero samples
group_time_table <- subset(
  group_time_table,
  Number_of_samples > 0
)

batch_group_table <- subset(
  batch_group_table,
  Number_of_samples > 0
)

batch_time_table <- subset(
  batch_time_table,
  Number_of_samples > 0
)

batch_group_time_table <- subset(
  batch_group_time_table,
  Number_of_samples > 0
)

cat("\nSamples by group and time point:\n")
print(group_time_table)

cat("\nSamples by batch and group:\n")
print(batch_group_table)

cat("\nSamples by batch, group and time point:\n")
print(batch_group_time_table)

# ============================================================
# Define the primary comparison
# ============================================================

# Primary comparison:
# ChILI samples collected at the time of liver injury
# versus control samples collected after ICI treatment.

primary_sample <- (
  metadata$Group == "chili" &
    metadata$Time_point == "TOL"
) | (
  metadata$Group == "control" &
    metadata$Time_point == "AC"
)

primary_metadata <- droplevels(
  metadata[
    primary_sample,
    ,
    drop = FALSE
  ]
)

primary_metadata$Primary_group <- ifelse(
  primary_metadata$Group == "chili" &
    primary_metadata$Time_point == "TOL",
  "ChILI_TOL",
  "Control_AC"
)

primary_metadata$Primary_group <- factor(
  primary_metadata$Primary_group,
  levels = c(
    "Control_AC",
    "ChILI_TOL"
  )
)

rownames(primary_metadata) <-
  primary_metadata$Sample_ID

# Select the same samples from the expression objects
primary_dge <- dge_object[
  ,
  primary_metadata$Sample_ID,
  keep.lib.sizes = FALSE
]

primary_dge <- calcNormFactors(
  primary_dge,
  method = "TMM"
)

primary_log_cpm <- log_cpm[
  ,
  primary_metadata$Sample_ID,
  drop = FALSE
]

# ============================================================
# Summarize the primary comparison
# ============================================================

primary_group_summary <- as.data.frame(
  table(primary_metadata$Primary_group)
)

colnames(primary_group_summary) <- c(
  "Primary_group",
  "Number_of_samples"
)

primary_batch_summary <- as.data.frame(
  table(
    primary_metadata$Batch,
    primary_metadata$Primary_group
  )
)

colnames(primary_batch_summary) <- c(
  "Batch",
  "Primary_group",
  "Number_of_samples"
)

primary_batch_summary_nonzero <- subset(
  primary_batch_summary,
  Number_of_samples > 0
)

cat("\nPrimary comparison sample numbers:\n")
print(primary_group_summary)

cat("\nPrimary comparison by sequencing batch:\n")
print(primary_batch_summary)

# ============================================================
# Check whether the statistical design is estimable
# ============================================================

primary_design <- model.matrix(
  ~ Primary_group + Batch,
  data = primary_metadata
)

design_rank <- qr(primary_design)$rank

design_columns <- ncol(primary_design)

design_is_full_rank <- design_rank == design_columns

design_summary <- data.frame(
  Measurement = c(
    "Number of primary samples",
    "Number of Control_AC samples",
    "Number of ChILI_TOL samples",
    "Number of batches",
    "Number of design columns",
    "Design matrix rank",
    "Full-rank design"
  ),
  
  Value = c(
    nrow(primary_metadata),
    sum(
      primary_metadata$Primary_group ==
        "Control_AC"
    ),
    sum(
      primary_metadata$Primary_group ==
        "ChILI_TOL"
    ),
    nlevels(primary_metadata$Batch),
    design_columns,
    design_rank,
    design_is_full_rank
  ),
  
  stringsAsFactors = FALSE
)

cat("\nStatistical design summary:\n")
print(design_summary)

if (!design_is_full_rank) {
  stop(
    "The primary statistical design is not full rank. ",
    "Batch and study group may be completely confounded."
  )
}

# ============================================================
# Identify batches represented by only one primary group
# ============================================================

batch_representation <- aggregate(
  Primary_group ~ Batch,
  data = primary_metadata,
  FUN = function(x) {
    length(unique(x))
  }
)

colnames(batch_representation)[2] <-
  "Number_of_primary_groups"

batch_representation$Interpretation <- ifelse(
  batch_representation$Number_of_primary_groups == 2,
  "Both primary groups are represented",
  "Only one primary group is represented"
)

cat("\nBatch representation in the primary comparison:\n")
print(batch_representation)

# ============================================================
# PCA before and after batch correction
# ============================================================

gene_variances <- apply(
  primary_log_cpm,
  MARGIN = 1,
  FUN = var
)

number_of_pca_genes <- min(
  5000,
  length(gene_variances)
)

top_variable_genes <- names(
  sort(
    gene_variances,
    decreasing = TRUE
  )
)[seq_len(number_of_pca_genes)]

primary_pca_raw <- prcomp(
  t(
    primary_log_cpm[
      top_variable_genes,
      ,
      drop = FALSE
    ]
  ),
  center = TRUE,
  scale. = FALSE
)

raw_variance <- round(
  100 *
    primary_pca_raw$sdev^2 /
    sum(primary_pca_raw$sdev^2),
  1
)

raw_pca_table <- data.frame(
  Sample_ID = rownames(primary_pca_raw$x),
  PC1 = primary_pca_raw$x[, 1],
  PC2 = primary_pca_raw$x[, 2],
  primary_metadata,
  stringsAsFactors = FALSE
)

raw_pca_plot <- ggplot(
  raw_pca_table,
  aes(
    x = PC1,
    y = PC2,
    color = Primary_group,
    shape = Batch
  )
) +
  geom_point(
    size = 3,
    alpha = 0.85
  ) +
  scale_color_manual(
    values = c(
      "Control_AC" = "#2166AC",
      "ChILI_TOL" = "#B2182B"
    )
  ) +
  labs(
    title = "Primary Samples before Batch Correction",
    x = paste0(
      "PC1 (",
      raw_variance[1],
      "%)"
    ),
    y = paste0(
      "PC2 (",
      raw_variance[2],
      "%)"
    ),
    color = "Primary group",
    shape = "Batch"
  ) +
  theme_classic(
    base_size = 11
  )

# Batch correction is used for PCA visualization only.
# Differential expression will use the original count data
# with Batch included in the statistical model.

group_design_for_visualization <- model.matrix(
  ~ Primary_group,
  data = primary_metadata
)

batch_corrected_log_cpm <- removeBatchEffect(
  primary_log_cpm,
  batch = primary_metadata$Batch,
  design = group_design_for_visualization
)

primary_pca_corrected <- prcomp(
  t(
    batch_corrected_log_cpm[
      top_variable_genes,
      ,
      drop = FALSE
    ]
  ),
  center = TRUE,
  scale. = FALSE
)

corrected_variance <- round(
  100 *
    primary_pca_corrected$sdev^2 /
    sum(primary_pca_corrected$sdev^2),
  1
)

corrected_pca_table <- data.frame(
  Sample_ID = rownames(primary_pca_corrected$x),
  PC1 = primary_pca_corrected$x[, 1],
  PC2 = primary_pca_corrected$x[, 2],
  primary_metadata,
  stringsAsFactors = FALSE
)

corrected_pca_plot <- ggplot(
  corrected_pca_table,
  aes(
    x = PC1,
    y = PC2,
    color = Primary_group,
    shape = Batch
  )
) +
  geom_point(
    size = 3,
    alpha = 0.85
  ) +
  scale_color_manual(
    values = c(
      "Control_AC" = "#2166AC",
      "ChILI_TOL" = "#B2182B"
    )
  ) +
  labs(
    title = "Primary Samples after Batch Correction",
    subtitle = "Batch correction is used for visualization only",
    x = paste0(
      "PC1 (",
      corrected_variance[1],
      "%)"
    ),
    y = paste0(
      "PC2 (",
      corrected_variance[2],
      "%)"
    ),
    color = "Primary group",
    shape = "Batch"
  ) +
  theme_classic(
    base_size = 11
  )

combined_pca_plot <- raw_pca_plot +
  corrected_pca_plot +
  plot_layout(
    guides = "collect"
  )

ggsave(
  filename = file.path(
    figures_directory,
    "Figure_06_primary_PCA_before_after_batch_correction.pdf"
  ),
  plot = combined_pca_plot,
  width = 15,
  height = 6
)

ggsave(
  filename = file.path(
    figures_directory,
    "Figure_06_primary_PCA_before_after_batch_correction.png"
  ),
  plot = combined_pca_plot,
  width = 15,
  height = 6,
  dpi = 300
)

# ============================================================
# Save primary analysis objects
# ============================================================

saveRDS(
  primary_dge,
  file = file.path(
    objects_directory,
    "GSE287540_primary_DGEList.rds"
  )
)

saveRDS(
  primary_metadata,
  file = file.path(
    objects_directory,
    "GSE287540_primary_metadata.rds"
  )
)

saveRDS(
  primary_design,
  file = file.path(
    objects_directory,
    "GSE287540_primary_design_matrix.rds"
  )
)

write.csv(
  primary_metadata,
  file = file.path(
    processed_data_directory,
    "GSE287540_primary_analysis_samples.csv"
  ),
  row.names = FALSE
)

write.csv(
  batch_representation,
  file = file.path(
    processed_data_directory,
    "GSE287540_primary_batch_representation.csv"
  ),
  row.names = FALSE
)

# ============================================================
# Create an Excel workbook
# ============================================================

review_workbook <- createWorkbook()

addWorksheet(
  review_workbook,
  "Design summary"
)

writeData(
  review_workbook,
  sheet = "Design summary",
  x = design_summary
)

addWorksheet(
  review_workbook,
  "Primary samples"
)

writeData(
  review_workbook,
  sheet = "Primary samples",
  x = primary_metadata
)

addWorksheet(
  review_workbook,
  "Primary groups"
)

writeData(
  review_workbook,
  sheet = "Primary groups",
  x = primary_group_summary
)

addWorksheet(
  review_workbook,
  "Primary batches"
)

writeData(
  review_workbook,
  sheet = "Primary batches",
  x = primary_batch_summary
)

addWorksheet(
  review_workbook,
  "Batch representation"
)

writeData(
  review_workbook,
  sheet = "Batch representation",
  x = batch_representation
)

addWorksheet(
  review_workbook,
  "All group and time"
)

writeData(
  review_workbook,
  sheet = "All group and time",
  x = group_time_table
)

addWorksheet(
  review_workbook,
  "All batch group time"
)

writeData(
  review_workbook,
  sheet = "All batch group time",
  x = batch_group_time_table
)

for (sheet_name in names(review_workbook)) {
  
  setColWidths(
    review_workbook,
    sheet = sheet_name,
    cols = 1:20,
    widths = "auto"
  )
}

saveWorkbook(
  review_workbook,
  file = file.path(
    tables_directory,
    "GSE287540_QC_review_and_primary_design.xlsx"
  ),
  overwrite = TRUE
)

# ============================================================
# Final summary
# ============================================================

cat("\n============================================\n")
cat("QC review and sample definition completed\n")
cat("============================================\n")

cat(
  "\nPrimary samples:",
  nrow(primary_metadata),
  "\n"
)

cat(
  "Control_AC samples:",
  sum(
    primary_metadata$Primary_group ==
      "Control_AC"
  ),
  "\n"
)

cat(
  "ChILI_TOL samples:",
  sum(
    primary_metadata$Primary_group ==
      "ChILI_TOL"
  ),
  "\n"
)

cat(
  "Full-rank statistical design:",
  design_is_full_rank,
  "\n"
)

cat("\nImportant analytical decision:\n")
cat("No samples were removed.\n")
cat("Batch will be included in the differential-expression model.\n")
cat("Batch-corrected expression was used only for PCA visualization.\n")

cat("\nFiles created:\n")
cat("results/R_objects/GSE287540_primary_DGEList.rds\n")
cat("results/R_objects/GSE287540_primary_metadata.rds\n")
cat("results/R_objects/GSE287540_primary_design_matrix.rds\n")
cat("data/processed/GSE287540_primary_analysis_samples.csv\n")
cat("results/tables/GSE287540_QC_review_and_primary_design.xlsx\n")
cat("results/figures/Figure_06_primary_PCA_before_after_batch_correction.png\n")