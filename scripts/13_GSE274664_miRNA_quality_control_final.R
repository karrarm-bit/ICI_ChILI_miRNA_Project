# ============================================================
# Script 13: Final Quality Control of GSE274664 miRNA Data
# Project: Integrated mRNA-miRNA Analysis of ICI Response
#
# Input:
# - Corrected human miRNA count matrix from Script 11B
# - Corrected sample metadata from Script 11B
#
# Purpose:
# 1. Validate the corrected miRNA count matrix
# 2. filter low-expression miRNAs
# 3. apply TMM normalization
# 4. calculate log2 CPM expression
# 5. identify potential low-quality samples
# 6. generate library-size, detection, PCA, and correlation plots
#
# Important:
# This script does not automatically exclude any sample.
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

if (!file.exists(count_matrix_file)) {
  stop(
    paste0(
      "The corrected miRNA count matrix was not found:\n",
      count_matrix_file
    )
  )
}

if (!file.exists(metadata_file)) {
  stop(
    paste0(
      "The corrected miRNA metadata were not found:\n",
      metadata_file
    )
  )
}


# ------------------------------------------------------------
# 5. Load the corrected data
# ------------------------------------------------------------

cat("Loading the corrected GSE274664 miRNA data...\n")

mirna_counts <- readRDS(
  count_matrix_file
)

mirna_metadata <- readRDS(
  metadata_file
)

mirna_counts <- as.matrix(
  mirna_counts
)

mirna_metadata <- as.data.frame(
  mirna_metadata
)

storage.mode(mirna_counts) <- "numeric"

cat(
  "Loaded matrix:",
  nrow(mirna_counts),
  "miRNAs and",
  ncol(mirna_counts),
  "samples\n"
)


# ------------------------------------------------------------
# 6. Validate the count matrix
# ------------------------------------------------------------

if (nrow(mirna_counts) == 0) {
  stop("The miRNA count matrix contains no rows.")
}

if (ncol(mirna_counts) == 0) {
  stop("The miRNA count matrix contains no columns.")
}

if (anyNA(mirna_counts)) {
  stop("The miRNA count matrix contains missing values.")
}

if (any(!is.finite(mirna_counts))) {
  stop("The miRNA count matrix contains non-finite values.")
}

if (any(mirna_counts < 0)) {
  stop("The miRNA count matrix contains negative values.")
}

if (anyDuplicated(
  rownames(mirna_counts)
) > 0) {
  stop("Duplicated miRNA names were detected.")
}

if (anyDuplicated(
  colnames(mirna_counts)
) > 0) {
  stop("Duplicated sample names were detected.")
}

if (ncol(mirna_counts) != 13) {
  warning(
    paste0(
      "Thirteen samples were expected, but ",
      ncol(mirna_counts),
      " were detected."
    )
  )
}


# ------------------------------------------------------------
# 7. Validate and arrange the metadata
# ------------------------------------------------------------

required_metadata_columns <- c(
  "Sample_ID",
  "Response_group"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(mirna_metadata)
)

if (length(missing_metadata_columns) > 0) {
  stop(
    paste0(
      "Missing metadata column(s): ",
      paste(
        missing_metadata_columns,
        collapse = ", "
      )
    )
  )
}

mirna_metadata$Sample_ID <- as.character(
  mirna_metadata$Sample_ID
)

mirna_metadata$Response_group <- factor(
  as.character(
    mirna_metadata$Response_group
  ),
  levels = c(
    "Non_responder",
    "Responder"
  )
)

if (anyNA(mirna_metadata$Response_group)) {
  stop(
    "Missing or invalid response-group values were detected."
  )
}

missing_samples <- setdiff(
  colnames(mirna_counts),
  mirna_metadata$Sample_ID
)

if (length(missing_samples) > 0) {
  stop(
    paste0(
      "Metadata are missing for sample(s): ",
      paste(missing_samples, collapse = ", ")
    )
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

if (!identical(
  colnames(mirna_counts),
  mirna_metadata$Sample_ID
)) {
  stop(
    "The count matrix and metadata are not correctly aligned."
  )
}

cat("\nSample groups:\n")
print(table(mirna_metadata$Response_group))


# ------------------------------------------------------------
# 8. Calculate raw sample statistics
# ------------------------------------------------------------

raw_library_sizes <- colSums(
  mirna_counts
)

raw_detected_mirnas <- colSums(
  mirna_counts > 0
)

if (any(raw_library_sizes <= 0)) {
  stop(
    paste0(
      "Zero-library sample(s) detected: ",
      paste(
        names(
          raw_library_sizes[
            raw_library_sizes <= 0
          ]
        ),
        collapse = ", "
      )
    )
  )
}

cat("\nRaw library sizes:\n")

print(
  data.frame(
    Sample_ID = names(raw_library_sizes),
    Library_size = as.numeric(
      raw_library_sizes
    ),
    Detected_miRNAs = as.numeric(
      raw_detected_mirnas
    ),
    row.names = NULL
  )
)


# ------------------------------------------------------------
# 9. Create an edgeR DGEList
# ------------------------------------------------------------

y <- DGEList(
  counts = mirna_counts,
  samples = data.frame(
    Sample_ID = mirna_metadata$Sample_ID,
    Response_group =
      mirna_metadata$Response_group,
    row.names = mirna_metadata$Sample_ID
  )
)


# ------------------------------------------------------------
# 10. Filter low-expression miRNAs
# ------------------------------------------------------------

cat("\nFiltering low-expression miRNAs...\n")

keep_mirnas <- filterByExpr(
  y,
  group = mirna_metadata$Response_group,
  min.count = 3,
  min.total.count = 10
)

number_before_filtering <- nrow(y)
number_retained <- sum(keep_mirnas)
number_removed <- sum(!keep_mirnas)

if (number_retained < 10) {
  stop(
    paste0(
      "Only ",
      number_retained,
      " miRNAs passed the expression filter. ",
      "The dataset requires further inspection."
    )
  )
}

y <- y[
  keep_mirnas,
  ,
  keep.lib.sizes = FALSE
]

filtered_library_sizes <- colSums(
  y$counts
)

if (any(filtered_library_sizes <= 0)) {
  stop(
    paste0(
      "A sample became zero after expression filtering: ",
      paste(
        names(
          filtered_library_sizes[
            filtered_library_sizes <= 0
          ]
        ),
        collapse = ", "
      )
    )
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
  "miRNAs retained:",
  number_retained,
  "\n"
)


# ------------------------------------------------------------
# 11. Apply TMM normalization
# ------------------------------------------------------------

cat("\nApplying TMM normalization...\n")

y <- calcNormFactors(
  y,
  method = "TMM"
)

if (
  any(!is.finite(y$samples$norm.factors)) ||
  any(y$samples$norm.factors <= 0)
) {
  stop(
    "Invalid TMM normalization factors were generated."
  )
}


# ------------------------------------------------------------
# 12. Calculate normalized log2 CPM
# ------------------------------------------------------------

log_cpm <- cpm(
  y,
  log = TRUE,
  prior.count = 2
)

if (anyNA(log_cpm)) {
  stop(
    "The normalized logCPM matrix contains missing values."
  )
}

if (any(!is.finite(log_cpm))) {
  stop(
    "The normalized logCPM matrix contains non-finite values."
  )
}


# ------------------------------------------------------------
# 13. Define a robust MAD outlier function
# ------------------------------------------------------------

flag_mad_outlier <- function(x, threshold = 3) {
  x <- as.numeric(x)
  
  median_value <- median(
    x,
    na.rm = TRUE
  )
  
  mad_value <- mad(
    x,
    center = median_value,
    constant = 1.4826,
    na.rm = TRUE
  )
  
  if (
    !is.finite(mad_value) ||
    mad_value == 0
  ) {
    return(
      rep(FALSE, length(x))
    )
  }
  
  abs(x - median_value) >
    threshold * mad_value
}


# ------------------------------------------------------------
# 14. Build the sample QC table
# ------------------------------------------------------------

qc_table <- data.frame(
  Sample_ID = colnames(y$counts),
  Sample_label = as.character(
    mirna_metadata[
      colnames(y$counts),
      "Sample_label"
    ]
  ),
  Response_group = as.character(
    mirna_metadata[
      colnames(y$counts),
      "Response_group"
    ]
  ),
  Raw_library_size = as.numeric(
    raw_library_sizes[
      colnames(y$counts)
    ]
  ),
  Raw_library_size_thousand =
    as.numeric(
      raw_library_sizes[
        colnames(y$counts)
      ]
    ) / 1000,
  Raw_detected_miRNAs = as.numeric(
    raw_detected_mirnas[
      colnames(y$counts)
    ]
  ),
  Filtered_library_size = as.numeric(
    colSums(y$counts)
  ),
  Filtered_detected_miRNAs = as.numeric(
    colSums(y$counts > 0)
  ),
  TMM_normalization_factor = as.numeric(
    y$samples$norm.factors
  ),
  stringsAsFactors = FALSE
)

qc_table$Effective_library_size <-
  qc_table$Filtered_library_size *
  qc_table$TMM_normalization_factor

qc_table$Library_size_MAD_outlier <-
  flag_mad_outlier(
    log10(
      qc_table$Raw_library_size + 1
    )
  )

qc_table$Detected_miRNA_MAD_outlier <-
  flag_mad_outlier(
    qc_table$Raw_detected_miRNAs
  )

qc_table$Very_low_library <- (
  qc_table$Raw_library_size < 1000
)

qc_table$Very_few_detected_miRNAs <- (
  qc_table$Raw_detected_miRNAs < 25
)

qc_table$Potential_QC_outlier <-
  qc_table$Library_size_MAD_outlier |
  qc_table$Detected_miRNA_MAD_outlier |
  qc_table$Very_low_library |
  qc_table$Very_few_detected_miRNAs

number_of_potential_outliers <- sum(
  qc_table$Potential_QC_outlier
)


# ------------------------------------------------------------
# 15. Create the library-size figure
# ------------------------------------------------------------

library_size_plot <- ggplot(
  qc_table,
  aes(
    x = reorder(
      Sample_ID,
      Raw_library_size_thousand
    ),
    y = Raw_library_size_thousand,
    fill = Response_group
  )
) +
  geom_col(
    width = 0.75,
    color = "black",
    linewidth = 0.25
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Non_responder" = "#D55E00",
      "Responder" = "#0072B2"
    )
  ) +
  labs(
    title = "GSE274664 Human miRNA Library Sizes",
    subtitle = "Before expression filtering",
    x = "Sample",
    y = "Human miRNA library size (thousand reads)",
    fill = "Response group"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = paste0(
    "results/figures/",
    "Figure_17_miRNA_library_sizes.png"
  ),
  plot = library_size_plot,
  width = 9,
  height = 6,
  dpi = 300
)

ggsave(
  filename = paste0(
    "results/figures/",
    "Figure_17_miRNA_library_sizes.pdf"
  ),
  plot = library_size_plot,
  width = 9,
  height = 6
)


# ------------------------------------------------------------
# 16. Create the detected-miRNA figure
# ------------------------------------------------------------

detected_plot <- ggplot(
  qc_table,
  aes(
    x = reorder(
      Sample_ID,
      Raw_detected_miRNAs
    ),
    y = Raw_detected_miRNAs,
    fill = Response_group
  )
) +
  geom_col(
    width = 0.75,
    color = "black",
    linewidth = 0.25
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Non_responder" = "#D55E00",
      "Responder" = "#0072B2"
    )
  ) +
  labs(
    title = "Detected Human miRNAs",
    subtitle = "Number of human miRNAs with non-zero counts",
    x = "Sample",
    y = "Detected human miRNAs",
    fill = "Response group"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = paste0(
    "results/figures/",
    "Figure_18_detected_miRNAs.png"
  ),
  plot = detected_plot,
  width = 9,
  height = 6,
  dpi = 300
)

ggsave(
  filename = paste0(
    "results/figures/",
    "Figure_18_detected_miRNAs.pdf"
  ),
  plot = detected_plot,
  width = 9,
  height = 6
)


# ------------------------------------------------------------
# 17. Perform principal component analysis
# ------------------------------------------------------------

mirna_variances <- apply(
  log_cpm,
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
    "Too few variable miRNAs are available for PCA."
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
    log_cpm[
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

pca_table$Sample_label <- qc_table$Sample_label[
  match(
    pca_table$Sample_ID,
    qc_table$Sample_ID
  )
]

pca_table$Response_group <-
  qc_table$Response_group[
    match(
      pca_table$Sample_ID,
      qc_table$Sample_ID
    )
  ]

pca_table$Potential_QC_outlier <-
  qc_table$Potential_QC_outlier[
    match(
      pca_table$Sample_ID,
      qc_table$Sample_ID
    )
  ]

pca_table$Plot_label <- ifelse(
  pca_table$Potential_QC_outlier,
  paste0(
    pca_table$Sample_label,
    " (",
    pca_table$Sample_ID,
    ")"
  ),
  NA_character_
)

pca_plot <- ggplot(
  pca_table,
  aes(
    x = PC1,
    y = PC2,
    color = Response_group
  )
) +
  geom_point(
    size = 4,
    alpha = 0.85
  ) +
  geom_text_repel(
    data = subset(
      pca_table,
      Potential_QC_outlier
    ),
    aes(label = Plot_label),
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
    title = "GSE274664 miRNA Principal Component Analysis",
    subtitle = paste0(
      "Based on the ",
      number_for_pca,
      " most variable human miRNAs"
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
    "Figure_19_miRNA_PCA.png"
  ),
  plot = pca_plot,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  filename = paste0(
    "results/figures/",
    "Figure_19_miRNA_PCA.pdf"
  ),
  plot = pca_plot,
  width = 10,
  height = 7
)


# ------------------------------------------------------------
# 18. Create the sample-correlation heatmap
# ------------------------------------------------------------

sample_correlation <- cor(
  log_cpm,
  method = "spearman",
  use = "pairwise.complete.obs"
)

correlation_annotation <- data.frame(
  Response_group = qc_table$Response_group[
    match(
      colnames(sample_correlation),
      qc_table$Sample_ID
    )
  ],
  row.names = colnames(sample_correlation),
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
    "Figure_20_miRNA_sample_correlation.png"
  ),
  width = 2800,
  height = 2400,
  res = 300
)

pheatmap(
  sample_correlation,
  annotation_col = correlation_annotation,
  annotation_row = correlation_annotation,
  annotation_colors = annotation_colors,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = "GSE274664 miRNA Sample Correlation",
  border_color = NA,
  fontsize = 9
)

dev.off()

pdf(
  file = paste0(
    "results/figures/",
    "Figure_20_miRNA_sample_correlation.pdf"
  ),
  width = 9,
  height = 8
)

pheatmap(
  sample_correlation,
  annotation_col = correlation_annotation,
  annotation_row = correlation_annotation,
  annotation_colors = annotation_colors,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = "GSE274664 miRNA Sample Correlation",
  border_color = NA,
  fontsize = 9
)

dev.off()


# ------------------------------------------------------------
# 19. Save R objects
# ------------------------------------------------------------

saveRDS(
  y,
  file = paste0(
    "results/R_objects/",
    "GSE274664_miRNA_DGEList_filtered_TMM.rds"
  )
)

saveRDS(
  log_cpm,
  file = paste0(
    "results/R_objects/",
    "GSE274664_miRNA_logCPM_normalized.rds"
  )
)

saveRDS(
  mirna_metadata,
  file = paste0(
    "results/R_objects/",
    "GSE274664_miRNA_metadata_QC_ready.rds"
  )
)

saveRDS(
  qc_table,
  file = paste0(
    "results/R_objects/",
    "GSE274664_miRNA_QC_table.rds"
  )
)


# ------------------------------------------------------------
# 20. Save CSV files
# ------------------------------------------------------------

write.csv(
  qc_table,
  file = paste0(
    "data/processed/",
    "GSE274664_miRNA_QC_table.csv"
  ),
  row.names = FALSE
)

write.csv(
  pca_table,
  file = paste0(
    "data/processed/",
    "GSE274664_miRNA_PCA_coordinates.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 21. Create the Excel workbook
# ------------------------------------------------------------

qc_workbook <- createWorkbook()

addWorksheet(
  qc_workbook,
  "QC_summary"
)

writeData(
  qc_workbook,
  sheet = "QC_summary",
  x = qc_table
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

addWorksheet(
  qc_workbook,
  "Normalization"
)

normalization_table <- data.frame(
  Sample_ID = rownames(y$samples),
  Library_size = y$samples$lib.size,
  TMM_factor = y$samples$norm.factors,
  Effective_library_size =
    y$samples$lib.size *
    y$samples$norm.factors,
  stringsAsFactors = FALSE
)

writeData(
  qc_workbook,
  sheet = "Normalization",
  x = normalization_table
)

saveWorkbook(
  qc_workbook,
  file = paste0(
    "results/tables/",
    "GSE274664_miRNA_quality_control_final.xlsx"
  ),
  overwrite = TRUE
)


# ------------------------------------------------------------
# 22. Print the final QC summary
# ------------------------------------------------------------

cat("\n")
cat("============================================\n")
cat("GSE274664 final miRNA quality control completed\n")
cat("============================================\n")

cat(
  "\nSamples analyzed:",
  ncol(y),
  "\n"
)

cat(
  "miRNAs before filtering:",
  number_before_filtering,
  "\n"
)

cat(
  "miRNAs retained:",
  number_retained,
  "\n"
)

cat(
  "Potential QC outliers:",
  number_of_potential_outliers,
  "\n"
)

cat("\nSamples per response group:\n")
print(table(mirna_metadata$Response_group))

if (number_of_potential_outliers > 0) {
  cat("\nPotential QC outlier samples:\n")
  
  print(
    qc_table[
      qc_table$Potential_QC_outlier,
      c(
        "Sample_ID",
        "Sample_label",
        "Response_group",
        "Raw_library_size",
        "Raw_detected_miRNAs",
        "Very_low_library",
        "Very_few_detected_miRNAs",
        "Library_size_MAD_outlier",
        "Detected_miRNA_MAD_outlier"
      ),
      drop = FALSE
    ],
    row.names = FALSE
  )
}

cat("\nImportant analytical note:\n")
cat("No samples were removed by this script.\n")
cat("Potential outliers must be reviewed using PCA and correlation.\n")

cat("\nMain files created:\n")
cat("results/R_objects/GSE274664_miRNA_DGEList_filtered_TMM.rds\n")
cat("results/R_objects/GSE274664_miRNA_logCPM_normalized.rds\n")
cat("results/R_objects/GSE274664_miRNA_metadata_QC_ready.rds\n")
cat("results/R_objects/GSE274664_miRNA_QC_table.rds\n")
cat("data/processed/GSE274664_miRNA_QC_table.csv\n")
cat("results/tables/GSE274664_miRNA_quality_control_final.xlsx\n")
cat("results/figures/Figure_17_miRNA_library_sizes.png\n")
cat("results/figures/Figure_18_detected_miRNAs.png\n")
cat("results/figures/Figure_19_miRNA_PCA.png\n")
cat("results/figures/Figure_20_miRNA_sample_correlation.png\n")