# ============================================================
# Script 05: Quality Control of GSE287540 RNA-seq Data
# Project: ICI-Induced Liver Injury and miRNA Analysis
# ============================================================

# Clear the R environment
rm(list = ls())

# Load required packages
suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
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

objects_directory <- file.path(
  project_directory,
  "results",
  "R_objects"
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

processed_data_directory <- file.path(
  project_directory,
  "data",
  "processed"
)

# ============================================================
# Define input files
# ============================================================

count_matrix_file <- file.path(
  objects_directory,
  "GSE287540_raw_count_matrix.rds"
)

metadata_file <- file.path(
  objects_directory,
  "GSE287540_metadata_analysis_ready.rds"
)

# Check required input files
required_files <- c(
  count_matrix_file,
  metadata_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  cat("\nMissing input files:\n")
  print(missing_files)
  
  stop(
    "Required files are missing. ",
    "Run 04_prepare_GSE287540_count_matrix.R first."
  )
}

# ============================================================
# Load count matrix and metadata
# ============================================================

cat("\nLoading count matrix and metadata...\n")

count_matrix <- readRDS(
  count_matrix_file
)

metadata <- readRDS(
  metadata_file
)

cat(
  "Loaded count matrix:",
  nrow(count_matrix),
  "genes and",
  ncol(count_matrix),
  "samples\n"
)

# ============================================================
# Validate sample alignment
# ============================================================

if (ncol(count_matrix) != nrow(metadata)) {
  stop(
    "The number of count-matrix columns does not match ",
    "the number of metadata rows."
  )
}

if (!identical(
  colnames(count_matrix),
  metadata$Sample_ID
)) {
  stop(
    "The count matrix and metadata are not aligned."
  )
}

# Assign sample identifiers as metadata row names
rownames(metadata) <- metadata$Sample_ID

# Prepare factors for analysis and plotting
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

metadata$Analysis_group <- interaction(
  metadata$Group,
  metadata$Time_point,
  drop = TRUE,
  sep = "_"
)

# ============================================================
# Create an edgeR object
# ============================================================

dge_object <- DGEList(
  counts = count_matrix,
  samples = metadata
)

# Store the original library sizes
original_library_sizes <- colSums(
  count_matrix
)

# ============================================================
# Filter low-expression genes
# ============================================================

cat("\nFiltering low-expression genes...\n")

keep_gene <- filterByExpr(
  dge_object,
  group = metadata$Analysis_group,
  min.count = 10,
  min.total.count = 20
)

genes_before_filtering <- nrow(
  dge_object
)

dge_filtered <- dge_object[
  keep_gene,
  ,
  keep.lib.sizes = FALSE
]

genes_after_filtering <- nrow(
  dge_filtered
)

genes_removed <- genes_before_filtering -
  genes_after_filtering

cat(
  "Genes before expression filtering:",
  genes_before_filtering,
  "\n"
)

cat(
  "Genes removed as low expression:",
  genes_removed,
  "\n"
)

cat(
  "Genes retained for analysis:",
  genes_after_filtering,
  "\n"
)

# ============================================================
# TMM normalization
# ============================================================

dge_filtered <- calcNormFactors(
  dge_filtered,
  method = "TMM"
)

# Calculate normalized log2 counts per million
log_cpm <- cpm(
  dge_filtered,
  normalized.lib.sizes = TRUE,
  log = TRUE,
  prior.count = 2
)

# ============================================================
# Calculate sample-level QC measurements
# ============================================================

qc_table <- data.frame(
  Sample_ID = metadata$Sample_ID,
  Library_ID = metadata$Library_ID,
  Subject_ID = metadata$Subject_ID,
  Group = metadata$Group,
  Time_point = metadata$Time_point,
  Batch = metadata$Batch,
  Treatment = metadata$Treatment,
  
  Raw_library_size = original_library_sizes,
  
  Library_size_million = round(
    original_library_sizes / 1000000,
    3
  ),
  
  Detected_genes_raw = colSums(
    count_matrix > 0
  ),
  
  Detected_genes_filtered = colSums(
    dge_filtered$counts > 0
  ),
  
  TMM_normalization_factor = round(
    dge_filtered$samples$norm.factors,
    5
  ),
  
  stringsAsFactors = FALSE
)

# ============================================================
# Identify potential QC outliers
# ============================================================

log_library_size <- log10(
  qc_table$Raw_library_size + 1
)

library_median <- median(
  log_library_size,
  na.rm = TRUE
)

library_mad <- mad(
  log_library_size,
  na.rm = TRUE
)

detected_median <- median(
  qc_table$Detected_genes_filtered,
  na.rm = TRUE
)

detected_mad <- mad(
  qc_table$Detected_genes_filtered,
  na.rm = TRUE
)

if (library_mad == 0) {
  library_mad <- 1
}

if (detected_mad == 0) {
  detected_mad <- 1
}

qc_table$Library_size_outlier <- abs(
  log_library_size - library_median
) > 3 * library_mad

qc_table$Detected_gene_outlier <- abs(
  qc_table$Detected_genes_filtered -
    detected_median
) > 3 * detected_mad

qc_table$Potential_QC_outlier <- (
  qc_table$Library_size_outlier |
    qc_table$Detected_gene_outlier
)

qc_table$QC_status <- ifelse(
  qc_table$Potential_QC_outlier,
  "Review",
  "Accept"
)

number_of_potential_outliers <- sum(
  qc_table$Potential_QC_outlier
)

cat(
  "\nPotential QC outliers detected:",
  number_of_potential_outliers,
  "\n"
)

cat(
  "Important: no samples are removed by this script.\n"
)

# ============================================================
# Figure 1: Raw library sizes
# ============================================================

library_plot_data <- qc_table

library_plot_data$Sample_order <- seq_len(
  nrow(library_plot_data)
)

library_size_plot <- ggplot(
  library_plot_data,
  aes(
    x = Sample_order,
    y = Library_size_million,
    color = Group
  )
) +
  geom_point(
    size = 2.2,
    alpha = 0.8
  ) +
  geom_point(
    data = subset(
      library_plot_data,
      Potential_QC_outlier
    ),
    shape = 21,
    fill = "yellow",
    color = "red",
    size = 4,
    stroke = 1
  ) +
  scale_color_manual(
    values = c(
      "control" = "#2166AC",
      "chili" = "#B2182B"
    )
  ) +
  labs(
    title = "Raw Library Sizes",
    subtitle = "Yellow circles indicate samples requiring review",
    x = "Sample order",
    y = "Library size (million reads)",
    color = "Study group"
  ) +
  theme_classic(
    base_size = 12
  )

ggsave(
  filename = file.path(
    figures_directory,
    "Figure_01_raw_library_sizes.pdf"
  ),
  plot = library_size_plot,
  width = 10,
  height = 6
)

ggsave(
  filename = file.path(
    figures_directory,
    "Figure_01_raw_library_sizes.png"
  ),
  plot = library_size_plot,
  width = 10,
  height = 6,
  dpi = 300
)

# ============================================================
# Figure 2: Number of detected genes
# ============================================================

detected_gene_plot <- ggplot(
  qc_table,
  aes(
    x = Group,
    y = Detected_genes_filtered,
    fill = Group
  )
) +
  geom_boxplot(
    width = 0.6,
    alpha = 0.7,
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(
      color = Time_point
    ),
    width = 0.15,
    size = 2,
    alpha = 0.8
  ) +
  scale_fill_manual(
    values = c(
      "control" = "#67A9CF",
      "chili" = "#EF8A62"
    )
  ) +
  labs(
    title = "Detected Genes after Expression Filtering",
    x = "Study group",
    y = "Number of detected genes",
    fill = "Study group",
    color = "Time point"
  ) +
  theme_classic(
    base_size = 12
  )

ggsave(
  filename = file.path(
    figures_directory,
    "Figure_02_detected_genes.pdf"
  ),
  plot = detected_gene_plot,
  width = 8,
  height = 6
)

ggsave(
  filename = file.path(
    figures_directory,
    "Figure_02_detected_genes.png"
  ),
  plot = detected_gene_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# ============================================================
# Principal component analysis
# ============================================================

cat("\nPerforming principal component analysis...\n")

# Select the 5,000 most variable genes
gene_variances <- apply(
  log_cpm,
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

pca_result <- prcomp(
  t(
    log_cpm[
      top_variable_genes,
      ,
      drop = FALSE
    ]
  ),
  center = TRUE,
  scale. = FALSE
)

explained_variance <- round(
  100 *
    pca_result$sdev^2 /
    sum(pca_result$sdev^2),
  1
)

pca_table <- data.frame(
  Sample_ID = rownames(pca_result$x),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  metadata,
  stringsAsFactors = FALSE
)

# ============================================================
# Figure 3: PCA by study group and time point
# ============================================================

pca_group_plot <- ggplot(
  pca_table,
  aes(
    x = PC1,
    y = PC2,
    color = Group,
    shape = Time_point
  )
) +
  geom_point(
    size = 3,
    alpha = 0.85
  ) +
  scale_color_manual(
    values = c(
      "control" = "#2166AC",
      "chili" = "#B2182B"
    )
  ) +
  labs(
    title = "Principal Component Analysis",
    subtitle = "Based on the 5,000 most variable genes",
    x = paste0(
      "PC1 (",
      explained_variance[1],
      "%)"
    ),
    y = paste0(
      "PC2 (",
      explained_variance[2],
      "%)"
    ),
    color = "Study group",
    shape = "Time point"
  ) +
  theme_classic(
    base_size = 12
  )

ggsave(
  filename = file.path(
    figures_directory,
    "Figure_03_PCA_group_time.pdf"
  ),
  plot = pca_group_plot,
  width = 9,
  height = 7
)

ggsave(
  filename = file.path(
    figures_directory,
    "Figure_03_PCA_group_time.png"
  ),
  plot = pca_group_plot,
  width = 9,
  height = 7,
  dpi = 300
)

# ============================================================
# Figure 4: PCA by sequencing batch
# ============================================================

pca_batch_plot <- ggplot(
  pca_table,
  aes(
    x = PC1,
    y = PC2,
    color = Batch,
    shape = Group
  )
) +
  geom_point(
    size = 3,
    alpha = 0.85
  ) +
  labs(
    title = "PCA Assessment of Sequencing Batch",
    subtitle = "The batch variable will be considered in later analyses",
    x = paste0(
      "PC1 (",
      explained_variance[1],
      "%)"
    ),
    y = paste0(
      "PC2 (",
      explained_variance[2],
      "%)"
    ),
    color = "Batch",
    shape = "Study group"
  ) +
  theme_classic(
    base_size = 12
  )

ggsave(
  filename = file.path(
    figures_directory,
    "Figure_04_PCA_batch.pdf"
  ),
  plot = pca_batch_plot,
  width = 9,
  height = 7
)

ggsave(
  filename = file.path(
    figures_directory,
    "Figure_04_PCA_batch.png"
  ),
  plot = pca_batch_plot,
  width = 9,
  height = 7,
  dpi = 300
)

# ============================================================
# Figure 5: Sample correlation heatmap
# ============================================================

sample_correlation <- cor(
  log_cpm,
  method = "pearson"
)

heatmap_annotation <- data.frame(
  Group = metadata$Group,
  Time_point = metadata$Time_point,
  Batch = metadata$Batch
)

rownames(heatmap_annotation) <- metadata$Sample_ID

pdf(
  file = file.path(
    figures_directory,
    "Figure_05_sample_correlation_heatmap.pdf"
  ),
  width = 12,
  height = 11
)

pheatmap(
  sample_correlation,
  annotation_col = heatmap_annotation,
  annotation_row = heatmap_annotation,
  show_colnames = FALSE,
  show_rownames = FALSE,
  border_color = NA,
  main = "Sample Correlation Heatmap"
)

dev.off()

# ============================================================
# Save QC tables and R objects
# ============================================================

write.csv(
  qc_table,
  file = file.path(
    processed_data_directory,
    "GSE287540_sample_QC_table.csv"
  ),
  row.names = FALSE
)

write.csv(
  pca_table,
  file = file.path(
    processed_data_directory,
    "GSE287540_PCA_coordinates.csv"
  ),
  row.names = FALSE
)

saveRDS(
  dge_filtered,
  file = file.path(
    objects_directory,
    "GSE287540_DGEList_filtered_TMM.rds"
  )
)

saveRDS(
  log_cpm,
  file = file.path(
    objects_directory,
    "GSE287540_logCPM_normalized.rds"
  )
)

saveRDS(
  pca_result,
  file = file.path(
    objects_directory,
    "GSE287540_PCA_result.rds"
  )
)

# ============================================================
# Create the QC Excel workbook
# ============================================================

qc_summary <- data.frame(
  Measurement = c(
    "Genes before expression filtering",
    "Low-expression genes removed",
    "Genes retained for analysis",
    "Number of samples",
    "Potential QC outliers",
    "PCA genes"
  ),
  
  Value = c(
    genes_before_filtering,
    genes_removed,
    genes_after_filtering,
    ncol(dge_filtered),
    number_of_potential_outliers,
    number_of_pca_genes
  ),
  
  stringsAsFactors = FALSE
)

qc_workbook <- createWorkbook()

addWorksheet(
  qc_workbook,
  "QC summary"
)

writeData(
  qc_workbook,
  sheet = "QC summary",
  x = qc_summary
)

addWorksheet(
  qc_workbook,
  "Sample QC"
)

writeData(
  qc_workbook,
  sheet = "Sample QC",
  x = qc_table
)

addWorksheet(
  qc_workbook,
  "PCA coordinates"
)

writeData(
  qc_workbook,
  sheet = "PCA coordinates",
  x = pca_table
)

for (sheet_name in names(qc_workbook)) {
  
  setColWidths(
    qc_workbook,
    sheet = sheet_name,
    cols = 1:20,
    widths = "auto"
  )
}

saveWorkbook(
  qc_workbook,
  file = file.path(
    tables_directory,
    "GSE287540_quality_control.xlsx"
  ),
  overwrite = TRUE
)

# ============================================================
# Final summary
# ============================================================

cat("\n============================================\n")
cat("GSE287540 quality control completed\n")
cat("============================================\n")

cat(
  "\nGenes retained:",
  genes_after_filtering,
  "\n"
)

cat(
  "Samples retained:",
  ncol(dge_filtered),
  "\n"
)

cat(
  "Potential QC outliers:",
  number_of_potential_outliers,
  "\n"
)

cat("\nImportant:\n")
cat("No samples were removed by this script.\n")
cat("Potential outliers must be reviewed before exclusion.\n")

cat("\nMain files created:\n")
cat("results/R_objects/GSE287540_DGEList_filtered_TMM.rds\n")
cat("results/R_objects/GSE287540_logCPM_normalized.rds\n")
cat("data/processed/GSE287540_sample_QC_table.csv\n")
cat("results/tables/GSE287540_quality_control.xlsx\n")
cat("results/figures/Figure_01_raw_library_sizes.pdf\n")
cat("results/figures/Figure_02_detected_genes.pdf\n")
cat("results/figures/Figure_03_PCA_group_time.pdf\n")
cat("results/figures/Figure_04_PCA_batch.pdf\n")
cat("results/figures/Figure_05_sample_correlation_heatmap.pdf\n")
qc_table <- read.csv(
  "data/processed/GSE287540_sample_QC_table.csv"
)

dge_object <- readRDS(
  "results/R_objects/GSE287540_DGEList_filtered_TMM.rds"
)

cat("Genes retained:", nrow(dge_object), "\n")
cat("Samples retained:", ncol(dge_object), "\n")
cat(
  "Potential QC outliers:",
  sum(qc_table$Potential_QC_outlier),
  "\n"
)

cat("\nQC status:\n")
print(table(qc_table$QC_status))

cat("\nPotential outlier samples:\n")
print(
  qc_table[
    qc_table$Potential_QC_outlier,
    c(
      "Sample_ID",
      "Subject_ID",
      "Group",
      "Time_point",
      "Batch",
      "Library_size_million",
      "Detected_genes_filtered"
    )
  ]
)
dge_object <- readRDS(
  "results/R_objects/GSE287540_DGEList_filtered_TMM.rds"
)

cat("Genes retained:", nrow(dge_object), "\n")
cat("Samples retained:", ncol(dge_object), "\n")