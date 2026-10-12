# ============================================================
# Script 14: Review miRNA QC and Define Final Samples
# Project: Integrated mRNA-miRNA Analysis of ICI Response
# Dataset: GSE274664
#
# Purpose:
# 1. Apply the documented QC decision
# 2. exclude the technically inadequate NR5 sample
# 3. retain all other samples
# 4. repeat expression filtering after sample exclusion
# 5. repeat TMM normalization
# 6. generate final PCA and correlation figures
# 7. save the final objects for differential expression
#
# QC decision:
# Hg19_699_CD19 is excluded because it contains only
# 28 human-miRNA reads and 20 detected human miRNAs.
# ============================================================


# ------------------------------------------------------------
# 1. Clean the R environment
# ------------------------------------------------------------

rm(list = ls())

options(stringsAsFactors = FALSE)

set.seed(123)


# ------------------------------------------------------------
# 2. Check required packages
# ------------------------------------------------------------

required_packages <- c(
  "edgeR",
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "openxlsx"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing package(s): ",
      paste(missing_packages, collapse = ", "),
      "\nRun script 00_install_packages.R first."
    )
  )
}

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(openxlsx)
})


# ------------------------------------------------------------
# 3. Create output directories
# ------------------------------------------------------------

output_directories <- c(
  "data/processed",
  "results/R_objects",
  "results/tables",
  "results/figures"
)

for (directory_name in output_directories) {
  if (!dir.exists(directory_name)) {
    dir.create(
      directory_name,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
}


# ------------------------------------------------------------
# 4. Define input files
# ------------------------------------------------------------

count_matrix_file <- paste0(
  "results/R_objects/",
  "GSE274664_mature_miRNA_count_matrix.rds"
)

metadata_file <- paste0(
  "results/R_objects/",
  "GSE274664_miRNA_metadata.rds"
)

qc_table_file <- paste0(
  "results/R_objects/",
  "GSE274664_miRNA_QC_table.rds"
)

required_input_files <- c(
  count_matrix_file,
  metadata_file,
  qc_table_file
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0) {
  stop(
    paste0(
      "The following required files were not found:\n",
      paste(
        missing_input_files,
        collapse = "\n"
      ),
      "\nRun script 13 first."
    )
  )
}


# ------------------------------------------------------------
# 5. Load the corrected count matrix and metadata
# ------------------------------------------------------------

cat("Loading GSE274664 data and QC results...\n")

mirna_counts <- readRDS(
  count_matrix_file
)

mirna_metadata <- readRDS(
  metadata_file
)

previous_qc_table <- readRDS(
  qc_table_file
)

mirna_counts <- as.matrix(
  mirna_counts
)

mirna_metadata <- as.data.frame(
  mirna_metadata
)

previous_qc_table <- as.data.frame(
  previous_qc_table
)

storage.mode(mirna_counts) <- "numeric"

cat(
  "Initial matrix:",
  nrow(mirna_counts),
  "miRNAs and",
  ncol(mirna_counts),
  "samples\n"
)


# ------------------------------------------------------------
# 6. Validate sample alignment
# ------------------------------------------------------------

if (!all(
  colnames(mirna_counts) %in%
  mirna_metadata$Sample_ID
)) {
  stop(
    "Some count-matrix samples are missing from the metadata."
  )
}

mirna_metadata <- mirna_metadata[
  match(
    colnames(mirna_counts),
    mirna_metadata$Sample_ID
  ),
  ,
  drop = FALSE
]

rownames(mirna_metadata) <-
  mirna_metadata$Sample_ID

mirna_metadata$Response_group <- factor(
  as.character(
    mirna_metadata$Response_group
  ),
  levels = c(
    "Non_responder",
    "Responder"
  )
)

if (!identical(
  colnames(mirna_counts),
  mirna_metadata$Sample_ID
)) {
  stop(
    "The count matrix and metadata are not aligned."
  )
}


# ------------------------------------------------------------
# 7. Define the documented technical exclusion
# ------------------------------------------------------------

excluded_sample <- "Hg19_699_CD19"

exclusion_reason <- paste0(
  "Extremely low human-miRNA library size: ",
  "28 reads with only 20 detected human miRNAs"
)

if (!excluded_sample %in% colnames(mirna_counts)) {
  stop(
    paste0(
      "The sample selected for exclusion was not found: ",
      excluded_sample
    )
  )
}

excluded_sample_qc <- previous_qc_table[
  previous_qc_table$Sample_ID ==
    excluded_sample,
  ,
  drop = FALSE
]

if (nrow(excluded_sample_qc) != 1) {
  stop(
    "The excluded sample could not be uniquely identified in the QC table."
  )
}

qc_decision_table <- data.frame(
  Sample_ID = colnames(mirna_counts),
  Sample_label = as.character(
    mirna_metadata[
      colnames(mirna_counts),
      "Sample_label"
    ]
  ),
  Response_group = as.character(
    mirna_metadata[
      colnames(mirna_counts),
      "Response_group"
    ]
  ),
  QC_decision = ifelse(
    colnames(mirna_counts) ==
      excluded_sample,
    "Exclude",
    "Retain"
  ),
  QC_reason = ifelse(
    colnames(mirna_counts) ==
      excluded_sample,
    exclusion_reason,
    "Retained after review of library size, PCA, and correlation"
  ),
  stringsAsFactors = FALSE
)

cat("\nDocumented sample exclusion:\n")

print(
  qc_decision_table[
    qc_decision_table$QC_decision ==
      "Exclude",
    ,
    drop = FALSE
  ],
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. Remove the technically inadequate sample
# ------------------------------------------------------------

samples_to_retain <- setdiff(
  colnames(mirna_counts),
  excluded_sample
)

final_counts <- mirna_counts[
  ,
  samples_to_retain,
  drop = FALSE
]

final_metadata <- mirna_metadata[
  samples_to_retain,
  ,
  drop = FALSE
]

final_metadata$Response_group <- droplevels(
  final_metadata$Response_group
)

if (!identical(
  colnames(final_counts),
  final_metadata$Sample_ID
)) {
  stop(
    "The final count matrix and metadata are not aligned."
  )
}

cat(
  "\nSamples retained after QC:",
  ncol(final_counts),
  "\n"
)

cat("\nFinal response groups:\n")
print(table(final_metadata$Response_group))


# ------------------------------------------------------------
# 9. Validate the final group sizes
# ------------------------------------------------------------

responder_number <- sum(
  final_metadata$Response_group ==
    "Responder"
)

non_responder_number <- sum(
  final_metadata$Response_group ==
    "Non_responder"
)

if (responder_number != 8) {
  stop(
    paste0(
      "Eight responders were expected, but ",
      responder_number,
      " were retained."
    )
  )
}

if (non_responder_number != 4) {
  stop(
    paste0(
      "Four non-responders were expected, but ",
      non_responder_number,
      " were retained."
    )
  )
}


# ------------------------------------------------------------
# 10. Create the final edgeR object
# ------------------------------------------------------------

final_y <- DGEList(
  counts = final_counts,
  samples = data.frame(
    Sample_ID = final_metadata$Sample_ID,
    Sample_label = final_metadata$Sample_label,
    Response_group =
      final_metadata$Response_group,
    row.names = final_metadata$Sample_ID
  )
)

if (any(final_y$samples$lib.size <= 0)) {
  stop(
    "A non-positive library size remains after QC."
  )
}


# ------------------------------------------------------------
# 11. Repeat expression filtering
# ------------------------------------------------------------

cat(
  "\nRepeating expression filtering after sample exclusion...\n"
)

keep_mirnas <- filterByExpr(
  final_y,
  group = final_metadata$Response_group,
  min.count = 3,
  min.total.count = 10
)

number_before_filtering <- nrow(
  final_y
)

number_retained <- sum(
  keep_mirnas
)

number_removed <- sum(
  !keep_mirnas
)

if (number_retained < 10) {
  stop(
    paste0(
      "Only ",
      number_retained,
      " miRNAs passed the final expression filter."
    )
  )
}

final_y <- final_y[
  keep_mirnas,
  ,
  keep.lib.sizes = FALSE
]

if (any(colSums(final_y$counts) <= 0)) {
  stop(
    "A zero filtered library was generated."
  )
}

cat(
  "miRNAs before filtering:",
  number_before_filtering,
  "\n"
)

cat(
  "Low-expression miRNAs removed:",
  number_removed,
  "\n"
)

cat(
  "Final miRNAs retained:",
  number_retained,
  "\n"
)


# ------------------------------------------------------------
# 12. Apply final TMM normalization
# ------------------------------------------------------------

cat("\nApplying final TMM normalization...\n")

if (
  "normLibSizes" %in%
  getNamespaceExports("edgeR")
) {
  final_y <- normLibSizes(
    final_y,
    method = "TMM"
  )
} else {
  final_y <- calcNormFactors(
    final_y,
    method = "TMM"
  )
}

if (
  any(!is.finite(
    final_y$samples$norm.factors
  )) ||
  any(
    final_y$samples$norm.factors <= 0
  )
) {
  stop(
    "Invalid final TMM normalization factors were generated."
  )
}


# ------------------------------------------------------------
# 13. Calculate final normalized log2 CPM
# ------------------------------------------------------------

final_log_cpm <- cpm(
  final_y,
  log = TRUE,
  prior.count = 2
)

if (
  anyNA(final_log_cpm) ||
  any(!is.finite(final_log_cpm))
) {
  stop(
    "The final logCPM matrix contains invalid values."
  )
}


# ------------------------------------------------------------
# 14. Create the statistical design matrix
# ------------------------------------------------------------

final_metadata$Response_group <- relevel(
  final_metadata$Response_group,
  ref = "Non_responder"
)

design_matrix <- model.matrix(
  ~ Response_group,
  data = final_metadata
)

rownames(design_matrix) <-
  final_metadata$Sample_ID

design_is_full_rank <- (
  qr(design_matrix)$rank ==
    ncol(design_matrix)
)

if (!design_is_full_rank) {
  stop(
    "The final statistical design matrix is not full rank."
  )
}

cat(
  "\nFull-rank statistical design:",
  design_is_full_rank,
  "\n"
)


# ------------------------------------------------------------
# 15. Perform final PCA
# ------------------------------------------------------------

mirna_variances <- apply(
  final_log_cpm,
  1,
  var,
  na.rm = TRUE
)

mirna_variances[
  !is.finite(mirna_variances)
] <- 0

number_for_pca <- min(
  100,
  sum(mirna_variances > 0)
)

if (number_for_pca < 2) {
  stop(
    "Too few variable miRNAs are available for final PCA."
  )
}

most_variable_mirnas <- names(
  sort(
    mirna_variances,
    decreasing = TRUE
  )
)[seq_len(number_for_pca)]

pca_result <- prcomp(
  t(
    final_log_cpm[
      most_variable_mirnas,
      ,
      drop = FALSE
    ]
  ),
  center = TRUE,
  scale. = FALSE
)

explained_variance <- (
  pca_result$sdev^2 /
    sum(pca_result$sdev^2)
) * 100

pca_table <- data.frame(
  Sample_ID = rownames(pca_result$x),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  stringsAsFactors = FALSE
)

pca_table$Sample_label <- as.character(
  final_metadata[
    pca_table$Sample_ID,
    "Sample_label"
  ]
)

pca_table$Response_group <- as.character(
  final_metadata[
    pca_table$Sample_ID,
    "Response_group"
  ]
)

final_pca_plot <- ggplot(
  pca_table,
  aes(
    x = PC1,
    y = PC2,
    color = Response_group,
    label = Sample_label
  )
) +
  geom_point(
    size = 4,
    alpha = 0.85
  ) +
  geom_text_repel(
    size = 3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Non_responder" = "#D55E00",
      "Responder" = "#0072B2"
    )
  ) +
  labs(
    title = "GSE274664 Final miRNA PCA",
    subtitle = paste0(
      "After exclusion of the technically inadequate NR5 sample; ",
      number_for_pca,
      " miRNAs used"
    ),
    x = paste0(
      "PC1 (",
      round(explained_variance[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(explained_variance[2], 1),
      "%)"
    ),
    color = "Response group"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = paste0(
    "results/figures/",
    "Figure_21_final_miRNA_PCA_after_QC.png"
  ),
  plot = final_pca_plot,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  filename = paste0(
    "results/figures/",
    "Figure_21_final_miRNA_PCA_after_QC.pdf"
  ),
  plot = final_pca_plot,
  width = 10,
  height = 7
)


# ------------------------------------------------------------
# 16. Create the final sample-correlation heatmap
# ------------------------------------------------------------

final_sample_correlation <- cor(
  final_log_cpm,
  method = "spearman",
  use = "pairwise.complete.obs"
)

correlation_annotation <- data.frame(
  Response_group = as.character(
    final_metadata[
      colnames(final_sample_correlation),
      "Response_group"
    ]
  ),
  row.names = colnames(
    final_sample_correlation
  ),
  stringsAsFactors = FALSE
)

annotation_colors <- list(
  Response_group = c(
    "Non_responder" = "#D55E00",
    "Responder" = "#0072B2"
  )
)

png(
  filename = paste0(
    "results/figures/",
    "Figure_22_final_miRNA_sample_correlation.png"
  ),
  width = 2800,
  height = 2400,
  res = 300
)

pheatmap(
  final_sample_correlation,
  annotation_col = correlation_annotation,
  annotation_row = correlation_annotation,
  annotation_colors = annotation_colors,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = "GSE274664 Final miRNA Sample Correlation",
  border_color = NA,
  fontsize = 9
)

dev.off()

pdf(
  file = paste0(
    "results/figures/",
    "Figure_22_final_miRNA_sample_correlation.pdf"
  ),
  width = 9,
  height = 8
)

pheatmap(
  final_sample_correlation,
  annotation_col = correlation_annotation,
  annotation_row = correlation_annotation,
  annotation_colors = annotation_colors,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = "GSE274664 Final miRNA Sample Correlation",
  border_color = NA,
  fontsize = 9
)

dev.off()


# ------------------------------------------------------------
# 17. Create the final sample summary
# ------------------------------------------------------------

final_sample_summary <- data.frame(
  Sample_ID = colnames(final_y$counts),
  Sample_label = as.character(
    final_metadata[
      colnames(final_y$counts),
      "Sample_label"
    ]
  ),
  Response_group = as.character(
    final_metadata[
      colnames(final_y$counts),
      "Response_group"
    ]
  ),
  Filtered_library_size = as.numeric(
    colSums(final_y$counts)
  ),
  Detected_filtered_miRNAs = as.numeric(
    colSums(final_y$counts > 0)
  ),
  TMM_normalization_factor = as.numeric(
    final_y$samples$norm.factors
  ),
  Effective_library_size = as.numeric(
    final_y$samples$lib.size *
      final_y$samples$norm.factors
  ),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 18. Save final R objects
# ------------------------------------------------------------

saveRDS(
  final_y,
  file = paste0(
    "results/R_objects/",
    "GSE274664_final_miRNA_DGEList.rds"
  )
)

saveRDS(
  final_log_cpm,
  file = paste0(
    "results/R_objects/",
    "GSE274664_final_miRNA_logCPM.rds"
  )
)

saveRDS(
  final_metadata,
  file = paste0(
    "results/R_objects/",
    "GSE274664_final_miRNA_metadata.rds"
  )
)

saveRDS(
  design_matrix,
  file = paste0(
    "results/R_objects/",
    "GSE274664_final_miRNA_design_matrix.rds"
  )
)


# ------------------------------------------------------------
# 19. Save CSV files
# ------------------------------------------------------------

write.csv(
  qc_decision_table,
  file = paste0(
    "data/processed/",
    "GSE274664_miRNA_QC_decisions.csv"
  ),
  row.names = FALSE
)

write.csv(
  final_metadata,
  file = paste0(
    "data/processed/",
    "GSE274664_final_miRNA_metadata.csv"
  ),
  row.names = FALSE
)

write.csv(
  final_sample_summary,
  file = paste0(
    "data/processed/",
    "GSE274664_final_miRNA_sample_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  pca_table,
  file = paste0(
    "data/processed/",
    "GSE274664_final_miRNA_PCA_coordinates.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 20. Create the QC decision workbook
# ------------------------------------------------------------

qc_workbook <- createWorkbook()

addWorksheet(
  qc_workbook,
  "QC_decisions"
)

writeData(
  qc_workbook,
  sheet = "QC_decisions",
  x = qc_decision_table
)

addWorksheet(
  qc_workbook,
  "Excluded_sample"
)

writeData(
  qc_workbook,
  sheet = "Excluded_sample",
  x = excluded_sample_qc
)

addWorksheet(
  qc_workbook,
  "Final_samples"
)

writeData(
  qc_workbook,
  sheet = "Final_samples",
  x = final_sample_summary
)

addWorksheet(
  qc_workbook,
  "Final_design"
)

design_export <- data.frame(
  Sample_ID = rownames(design_matrix),
  design_matrix,
  check.names = FALSE,
  row.names = NULL
)

writeData(
  qc_workbook,
  sheet = "Final_design",
  x = design_export
)

addWorksheet(
  qc_workbook,
  "PCA_coordinates"
)

writeData(
  qc_workbook,
  sheet = "PCA_coordinates",
  x = pca_table
)

saveWorkbook(
  qc_workbook,
  file = paste0(
    "results/tables/",
    "GSE274664_QC_review_and_final_samples.xlsx"
  ),
  overwrite = TRUE
)


# ------------------------------------------------------------
# 21. Print final confirmation
# ------------------------------------------------------------

cat("\n")
cat("============================================\n")
cat("GSE274664 QC review completed\n")
cat("============================================\n")

cat(
  "\nExcluded sample:",
  excluded_sample,
  "\n"
)

cat(
  "Exclusion reason:",
  exclusion_reason,
  "\n"
)

cat(
  "\nFinal samples:",
  ncol(final_y),
  "\n"
)

cat(
  "Responders:",
  responder_number,
  "\n"
)

cat(
  "Non-responders:",
  non_responder_number,
  "\n"
)

cat(
  "Final miRNAs retained:",
  nrow(final_y),
  "\n"
)

cat(
  "Full-rank statistical design:",
  design_is_full_rank,
  "\n"
)

cat("\nImportant analytical decision:\n")
cat("NR5 was excluded because of an extremely low miRNA library.\n")
cat("R4, R7, and NR4 were retained after PCA and correlation review.\n")
cat("TMM normalization was repeated after sample exclusion.\n")

cat("\nMain files created:\n")
cat("results/R_objects/GSE274664_final_miRNA_DGEList.rds\n")
cat("results/R_objects/GSE274664_final_miRNA_logCPM.rds\n")
cat("results/R_objects/GSE274664_final_miRNA_metadata.rds\n")
cat("results/R_objects/GSE274664_final_miRNA_design_matrix.rds\n")
cat("data/processed/GSE274664_miRNA_QC_decisions.csv\n")
cat("data/processed/GSE274664_final_miRNA_metadata.csv\n")
cat("results/tables/GSE274664_QC_review_and_final_samples.xlsx\n")
cat("results/figures/Figure_21_final_miRNA_PCA_after_QC.png\n")
cat("results/figures/Figure_22_final_miRNA_sample_correlation.png\n")