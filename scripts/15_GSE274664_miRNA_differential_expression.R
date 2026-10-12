# ============================================================
# Script 15: GSE274664 miRNA Differential Expression
# Project: Integrated mRNA-miRNA Analysis of ICI Response
#
# Comparison:
# Responder versus Non_responder
#
# Interpretation:
# Positive logFC = higher expression in Responders
# Negative logFC = higher expression in Non-responders
#
# Primary threshold:
# FDR < 0.05 and absolute log2 fold change >= 1
#
# Exploratory threshold:
# P-value < 0.05 and absolute log2 fold change >= 1
# ============================================================


# ------------------------------------------------------------
# 1. Clean the R environment
# ------------------------------------------------------------

rm(list = ls())

options(stringsAsFactors = FALSE)

set.seed(123)


# ------------------------------------------------------------
# 2. Load required packages
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

dge_file <- paste0(
  "results/R_objects/",
  "GSE274664_final_miRNA_DGEList.rds"
)

metadata_file <- paste0(
  "results/R_objects/",
  "GSE274664_final_miRNA_metadata.rds"
)

log_cpm_file <- paste0(
  "results/R_objects/",
  "GSE274664_final_miRNA_logCPM.rds"
)

required_files <- c(
  dge_file,
  metadata_file,
  log_cpm_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    paste0(
      "Missing input file(s):\n",
      paste(missing_files, collapse = "\n"),
      "\nRun script 14 first."
    )
  )
}


# ------------------------------------------------------------
# 5. Load final analysis-ready data
# ------------------------------------------------------------

cat("Loading final GSE274664 miRNA data...\n")

y <- readRDS(dge_file)

mirna_metadata <- readRDS(
  metadata_file
)

log_cpm <- readRDS(
  log_cpm_file
)

mirna_metadata <- as.data.frame(
  mirna_metadata
)

log_cpm <- as.matrix(
  log_cpm
)

cat(
  "Loaded:",
  nrow(y),
  "miRNAs and",
  ncol(y),
  "samples\n"
)


# ------------------------------------------------------------
# 6. Validate and align metadata
# ------------------------------------------------------------

required_metadata_columns <- c(
  "Sample_ID",
  "Sample_label",
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

if (!all(
  colnames(y$counts) %in%
  mirna_metadata$Sample_ID
)) {
  stop(
    "Some DGEList samples are missing from the metadata."
  )
}

mirna_metadata <- mirna_metadata[
  match(
    colnames(y$counts),
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

if (anyNA(mirna_metadata$Response_group)) {
  stop(
    "Missing or invalid response-group values were detected."
  )
}

if (!identical(
  colnames(y$counts),
  mirna_metadata$Sample_ID
)) {
  stop(
    "The DGEList and metadata are not aligned."
  )
}

if (!identical(
  colnames(y$counts),
  colnames(log_cpm)
)) {
  stop(
    "The DGEList and logCPM samples are not aligned."
  )
}

if (!identical(
  rownames(y$counts),
  rownames(log_cpm)
)) {
  stop(
    "The DGEList and logCPM miRNAs are not aligned."
  )
}

cat("\nSamples per response group:\n")
print(table(mirna_metadata$Response_group))


# ------------------------------------------------------------
# 7. Construct the statistical design
# ------------------------------------------------------------

design_matrix <- model.matrix(
  ~ Response_group,
  data = mirna_metadata
)

rownames(design_matrix) <-
  mirna_metadata$Sample_ID

design_is_full_rank <- (
  qr(design_matrix)$rank ==
    ncol(design_matrix)
)

if (!design_is_full_rank) {
  stop(
    "The differential-expression design is not full rank."
  )
}

coefficient_name <- "Response_groupResponder"

if (!coefficient_name %in%
    colnames(design_matrix)) {
  stop(
    paste0(
      "The expected comparison coefficient was not found: ",
      coefficient_name
    )
  )
}

cat("\nStatistical design matrix:\n")
print(design_matrix)


# ------------------------------------------------------------
# 8. Estimate dispersion
# ------------------------------------------------------------

cat("\nEstimating miRNA biological dispersion...\n")

y <- estimateDisp(
  y,
  design = design_matrix,
  robust = TRUE
)

common_dispersion <- y$common.dispersion

biological_coefficient_variation <- sqrt(
  common_dispersion
)

cat(
  "Common dispersion:",
  round(common_dispersion, 5),
  "\n"
)

cat(
  "Biological coefficient of variation:",
  round(
    biological_coefficient_variation,
    4
  ),
  "\n"
)


# ------------------------------------------------------------
# 9. Fit the quasi-likelihood model
# ------------------------------------------------------------

cat("\nFitting the quasi-likelihood model...\n")

ql_fit <- glmQLFit(
  y,
  design = design_matrix,
  robust = TRUE
)

ql_test <- glmQLFTest(
  ql_fit,
  coef = coefficient_name
)


# ------------------------------------------------------------
# 10. Extract differential-expression results
# ------------------------------------------------------------

de_table <- topTags(
  ql_test,
  n = Inf,
  sort.by = "PValue"
)$table

de_results <- data.frame(
  miRNA = rownames(de_table),
  de_table,
  row.names = NULL,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_result_columns <- c(
  "miRNA",
  "logFC",
  "logCPM",
  "F",
  "PValue",
  "FDR"
)

if (!all(
  required_result_columns %in%
  colnames(de_results)
)) {
  stop(
    "The expected differential-expression columns were not generated."
  )
}


# ------------------------------------------------------------
# 11. Apply significance thresholds
# ------------------------------------------------------------

fdr_threshold <- 0.05
pvalue_threshold <- 0.05
logfc_threshold <- 1

de_results$Fold_change <- 2^(
  de_results$logFC
)

de_results$Primary_significant <- (
  de_results$FDR < fdr_threshold &
    abs(de_results$logFC) >=
    logfc_threshold
)

de_results$Nominal_significant <- (
  de_results$PValue < pvalue_threshold &
    abs(de_results$logFC) >=
    logfc_threshold
)

de_results$Direction <- ifelse(
  de_results$Primary_significant &
    de_results$logFC > 0,
  "Higher_in_Responder",
  ifelse(
    de_results$Primary_significant &
      de_results$logFC < 0,
    "Higher_in_Non_responder",
    ifelse(
      de_results$Nominal_significant &
        de_results$logFC > 0,
      "Nominal_higher_in_Responder",
      ifelse(
        de_results$Nominal_significant &
          de_results$logFC < 0,
        "Nominal_higher_in_Non_responder",
        "Not_significant"
      )
    )
  )
)

de_results$Comparison <-
  "Responder_vs_Non_responder"

de_results$Volcano_category <- factor(
  de_results$Direction,
  levels = c(
    "Higher_in_Responder",
    "Higher_in_Non_responder",
    "Nominal_higher_in_Responder",
    "Nominal_higher_in_Non_responder",
    "Not_significant"
  )
)


# ------------------------------------------------------------
# 12. Create result subsets
# ------------------------------------------------------------

primary_results <- de_results[
  de_results$Primary_significant,
  ,
  drop = FALSE
]

nominal_results <- de_results[
  de_results$Nominal_significant,
  ,
  drop = FALSE
]

primary_upregulated <- primary_results[
  primary_results$logFC >=
    logfc_threshold,
  ,
  drop = FALSE
]

primary_downregulated <- primary_results[
  primary_results$logFC <=
    -logfc_threshold,
  ,
  drop = FALSE
]

nominal_upregulated <- nominal_results[
  nominal_results$logFC >=
    logfc_threshold,
  ,
  drop = FALSE
]

nominal_downregulated <- nominal_results[
  nominal_results$logFC <=
    -logfc_threshold,
  ,
  drop = FALSE
]

number_primary <- nrow(
  primary_results
)

number_primary_up <- nrow(
  primary_upregulated
)

number_primary_down <- nrow(
  primary_downregulated
)

number_nominal <- nrow(
  nominal_results
)

number_nominal_up <- nrow(
  nominal_upregulated
)

number_nominal_down <- nrow(
  nominal_downregulated
)


# ------------------------------------------------------------
# 13. Select miRNAs for volcano labels
# ------------------------------------------------------------

label_candidates <- de_results[
  de_results$Nominal_significant,
  ,
  drop = FALSE
]

label_candidates <- label_candidates[
  order(
    label_candidates$FDR,
    label_candidates$PValue
  ),
  ,
  drop = FALSE
]

if (nrow(label_candidates) == 0) {
  label_candidates <- de_results[
    order(
      de_results$PValue
    ),
    ,
    drop = FALSE
  ]
}

number_of_labels <- min(
  10,
  nrow(label_candidates)
)

label_candidates <- head(
  label_candidates,
  number_of_labels
)

label_candidates$Volcano_y <- -log10(
  pmax(
    label_candidates$FDR,
    1e-300
  )
)


# ------------------------------------------------------------
# 14. Create the corrected volcano plot
# ------------------------------------------------------------

volcano_plot <- ggplot(
  de_results,
  aes(
    x = logFC,
    y = -log10(
      pmax(FDR, 1e-300)
    ),
    color = Volcano_category
  )
) +
  geom_point(
    size = 3,
    alpha = 0.80
  ) +
  geom_vline(
    xintercept = c(
      -logfc_threshold,
      logfc_threshold
    ),
    linetype = "dashed",
    color = "grey40"
  ) +
  geom_hline(
    yintercept = -log10(
      fdr_threshold
    ),
    linetype = "dashed",
    color = "grey40"
  ) +
  geom_text_repel(
    data = label_candidates,
    mapping = aes(
      x = logFC,
      y = Volcano_y,
      label = miRNA
    ),
    inherit.aes = FALSE,
    size = 3.2,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Higher_in_Responder" = "#D73027",
      "Higher_in_Non_responder" = "#2166AC",
      "Nominal_higher_in_Responder" = "#FC8D59",
      "Nominal_higher_in_Non_responder" = "#67A9CF",
      "Not_significant" = "grey70"
    ),
    drop = FALSE
  ) +
  labs(
    title = "GSE274664 Differentially Expressed miRNAs",
    subtitle = "Responder versus Non-responder",
    x = "log2 fold change",
    y = "-log10 FDR",
    color = "Classification"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(
  filename = paste0(
    "results/figures/",
    "Figure_23_miRNA_volcano_Responder_vs_Non_responder.png"
  ),
  plot = volcano_plot,
  width = 11,
  height = 7,
  dpi = 300
)

ggsave(
  filename = paste0(
    "results/figures/",
    "Figure_23_miRNA_volcano_Responder_vs_Non_responder.pdf"
  ),
  plot = volcano_plot,
  width = 11,
  height = 7
)


# ------------------------------------------------------------
# 15. Select top miRNAs for visualization
# ------------------------------------------------------------

ordered_results <- de_results[
  order(
    de_results$FDR,
    de_results$PValue
  ),
  ,
  drop = FALSE
]

number_for_heatmap <- min(
  20,
  nrow(ordered_results)
)

heatmap_mirnas <- head(
  ordered_results$miRNA,
  number_for_heatmap
)


# ------------------------------------------------------------
# 16. Create the top-miRNA heatmap
# ------------------------------------------------------------

heatmap_matrix <- log_cpm[
  heatmap_mirnas,
  ,
  drop = FALSE
]

heatmap_scaled <- t(
  scale(
    t(heatmap_matrix)
  )
)

heatmap_scaled[
  !is.finite(heatmap_scaled)
] <- 0

heatmap_annotation <- data.frame(
  Response_group = as.character(
    mirna_metadata[
      colnames(heatmap_scaled),
      "Response_group"
    ]
  ),
  row.names = colnames(
    heatmap_scaled
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
    "Figure_24_top20_miRNA_heatmap.png"
  ),
  width = 2800,
  height = 2500,
  res = 300
)

pheatmap(
  heatmap_scaled,
  annotation_col = heatmap_annotation,
  annotation_colors = annotation_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 8,
  border_color = NA,
  main = "Top Differentially Expressed miRNAs"
)

dev.off()

pdf(
  file = paste0(
    "results/figures/",
    "Figure_24_top20_miRNA_heatmap.pdf"
  ),
  width = 9,
  height = 8
)

pheatmap(
  heatmap_scaled,
  annotation_col = heatmap_annotation,
  annotation_colors = annotation_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 8,
  border_color = NA,
  main = "Top Differentially Expressed miRNAs"
)

dev.off()


# ------------------------------------------------------------
# 17. Prepare top-candidate expression data
# ------------------------------------------------------------

number_for_boxplot <- min(
  6,
  nrow(ordered_results)
)

boxplot_mirnas <- head(
  ordered_results$miRNA,
  number_for_boxplot
)

boxplot_list <- lapply(
  boxplot_mirnas,
  function(mirna_name) {
    data.frame(
      Sample_ID = colnames(log_cpm),
      Sample_label = as.character(
        mirna_metadata[
          colnames(log_cpm),
          "Sample_label"
        ]
      ),
      Response_group = as.character(
        mirna_metadata[
          colnames(log_cpm),
          "Response_group"
        ]
      ),
      miRNA = mirna_name,
      logCPM = as.numeric(
        log_cpm[
          mirna_name,
          colnames(log_cpm)
        ]
      ),
      stringsAsFactors = FALSE
    )
  }
)

boxplot_data <- do.call(
  rbind,
  boxplot_list
)

boxplot_data$Response_group <- factor(
  boxplot_data$Response_group,
  levels = c(
    "Non_responder",
    "Responder"
  )
)


# ------------------------------------------------------------
# 18. Create top-candidate expression plots
# ------------------------------------------------------------

candidate_boxplot <- ggplot(
  boxplot_data,
  aes(
    x = Response_group,
    y = logCPM,
    fill = Response_group
  )
) +
  geom_boxplot(
    width = 0.60,
    alpha = 0.65,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 2.3,
    alpha = 0.85
  ) +
  facet_wrap(
    ~ miRNA,
    scales = "free_y",
    ncol = 3
  ) +
  scale_fill_manual(
    values = c(
      "Non_responder" = "#D55E00",
      "Responder" = "#0072B2"
    )
  ) +
  labs(
    title = "Expression of Top Candidate miRNAs",
    subtitle = "TMM-normalized log2 CPM expression",
    x = "Response group",
    y = "log2 CPM"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    ),
    legend.position = "none"
  )

ggsave(
  filename = paste0(
    "results/figures/",
    "Figure_25_top_candidate_miRNA_expression.png"
  ),
  plot = candidate_boxplot,
  width = 12,
  height = 8,
  dpi = 300
)

ggsave(
  filename = paste0(
    "results/figures/",
    "Figure_25_top_candidate_miRNA_expression.pdf"
  ),
  plot = candidate_boxplot,
  width = 12,
  height = 8
)


# ------------------------------------------------------------
# 19. Save CSV result tables
# ------------------------------------------------------------

write.csv(
  de_results,
  file = paste0(
    "data/processed/",
    "GSE274664_all_miRNA_DE_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  primary_results,
  file = paste0(
    "data/processed/",
    "GSE274664_FDR_significant_miRNAs.csv"
  ),
  row.names = FALSE
)

write.csv(
  nominal_results,
  file = paste0(
    "data/processed/",
    "GSE274664_nominal_candidate_miRNAs.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 20. Create the Excel workbook
# ------------------------------------------------------------

analysis_summary <- data.frame(
  Measurement = c(
    "Comparison",
    "Samples analyzed",
    "Responder samples",
    "Non-responder samples",
    "miRNAs tested",
    "FDR-significant miRNAs",
    "FDR higher in Responders",
    "FDR higher in Non-responders",
    "Nominal candidate miRNAs",
    "Nominal higher in Responders",
    "Nominal higher in Non-responders",
    "Common dispersion",
    "Biological coefficient of variation",
    "FDR threshold",
    "Exploratory P-value threshold",
    "Absolute log2FC threshold"
  ),
  Value = c(
    "Responder versus Non-responder",
    ncol(y),
    sum(
      mirna_metadata$Response_group ==
        "Responder"
    ),
    sum(
      mirna_metadata$Response_group ==
        "Non_responder"
    ),
    nrow(de_results),
    number_primary,
    number_primary_up,
    number_primary_down,
    number_nominal,
    number_nominal_up,
    number_nominal_down,
    common_dispersion,
    biological_coefficient_variation,
    fdr_threshold,
    pvalue_threshold,
    logfc_threshold
  ),
  stringsAsFactors = FALSE
)

de_workbook <- createWorkbook()

addWorksheet(
  de_workbook,
  "Analysis_summary"
)

writeData(
  de_workbook,
  sheet = "Analysis_summary",
  x = analysis_summary
)

addWorksheet(
  de_workbook,
  "All_miRNAs"
)

writeData(
  de_workbook,
  sheet = "All_miRNAs",
  x = de_results
)

addWorksheet(
  de_workbook,
  "FDR_significant"
)

if (nrow(primary_results) > 0) {
  writeData(
    de_workbook,
    sheet = "FDR_significant",
    x = primary_results
  )
} else {
  writeData(
    de_workbook,
    sheet = "FDR_significant",
    x = data.frame(
      Message = paste0(
        "No miRNAs met FDR < 0.05 and ",
        "absolute log2FC >= 1."
      )
    )
  )
}

addWorksheet(
  de_workbook,
  "Nominal_candidates"
)

if (nrow(nominal_results) > 0) {
  writeData(
    de_workbook,
    sheet = "Nominal_candidates",
    x = nominal_results
  )
} else {
  writeData(
    de_workbook,
    sheet = "Nominal_candidates",
    x = data.frame(
      Message = paste0(
        "No miRNAs met P-value < 0.05 and ",
        "absolute log2FC >= 1."
      )
    )
  )
}

addWorksheet(
  de_workbook,
  "Top20"
)

writeData(
  de_workbook,
  sheet = "Top20",
  x = head(
    ordered_results,
    20
  )
)

saveWorkbook(
  de_workbook,
  file = paste0(
    "results/tables/",
    "GSE274664_miRNA_differential_expression.xlsx"
  ),
  overwrite = TRUE
)


# ------------------------------------------------------------
# 21. Save analysis R objects
# ------------------------------------------------------------

saveRDS(
  y,
  file = paste0(
    "results/R_objects/",
    "GSE274664_miRNA_DGEList_with_dispersion.rds"
  )
)

saveRDS(
  ql_fit,
  file = paste0(
    "results/R_objects/",
    "GSE274664_miRNA_QL_model.rds"
  )
)

saveRDS(
  ql_test,
  file = paste0(
    "results/R_objects/",
    "GSE274664_miRNA_QL_test.rds"
  )
)

saveRDS(
  de_results,
  file = paste0(
    "results/R_objects/",
    "GSE274664_miRNA_DE_results.rds"
  )
)


# ------------------------------------------------------------
# 22. Print final confirmation
# ------------------------------------------------------------

cat("\n")
cat("============================================\n")
cat("GSE274664 miRNA differential expression completed\n")
cat("============================================\n")

cat(
  "\nComparison: Responder versus Non-responder\n"
)

cat(
  "Samples analyzed:",
  ncol(y),
  "\n"
)

cat(
  "Responder samples:",
  sum(
    mirna_metadata$Response_group ==
      "Responder"
  ),
  "\n"
)

cat(
  "Non-responder samples:",
  sum(
    mirna_metadata$Response_group ==
      "Non_responder"
  ),
  "\n"
)

cat(
  "miRNAs tested:",
  nrow(de_results),
  "\n"
)

cat(
  "FDR-significant miRNAs:",
  number_primary,
  "\n"
)

cat(
  "FDR higher in Responders:",
  number_primary_up,
  "\n"
)

cat(
  "FDR higher in Non-responders:",
  number_primary_down,
  "\n"
)

cat(
  "Nominal candidate miRNAs:",
  number_nominal,
  "\n"
)

cat(
  "Nominal higher in Responders:",
  number_nominal_up,
  "\n"
)

cat(
  "Nominal higher in Non-responders:",
  number_nominal_down,
  "\n"
)

cat("\nThresholds:\n")
cat("Primary: FDR < 0.05 and |log2FC| >= 1\n")
cat("Exploratory: P-value < 0.05 and |log2FC| >= 1\n")

cat("\nInterpretation:\n")
cat("Positive logFC indicates higher expression in Responders.\n")
cat("Negative logFC indicates higher expression in Non-responders.\n")
cat("Nominal candidates require integration or validation.\n")

cat("\nMain files created:\n")
cat("results/tables/GSE274664_miRNA_differential_expression.xlsx\n")
cat("results/R_objects/GSE274664_miRNA_DE_results.rds\n")
cat("results/figures/Figure_23_miRNA_volcano_Responder_vs_Non_responder.png\n")
cat("results/figures/Figure_24_top20_miRNA_heatmap.png\n")
cat("results/figures/Figure_25_top_candidate_miRNA_expression.png\n")
de_results <- readRDS(
  "results/R_objects/GSE274664_miRNA_DE_results.rds"
)

cat("\nFDR-significant miRNAs:\n")

print(
  de_results[
    de_results$Primary_significant,
    c(
      "miRNA",
      "logFC",
      "Fold_change",
      "logCPM",
      "PValue",
      "FDR",
      "Direction"
    )
  ],
  row.names = FALSE
)

cat("\nTop 20 miRNAs:\n")

print(
  head(
    de_results[
      order(
        de_results$PValue
      ),
      c(
        "miRNA",
        "logFC",
        "Fold_change",
        "PValue",
        "FDR",
        "Direction"
      )
    ],
    20
  ),
  row.names = FALSE
)