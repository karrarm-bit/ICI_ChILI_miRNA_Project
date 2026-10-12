# ============================================================
# Script 07: Differential Expression Analysis of GSE287540
# Comparison: ChILI_TOL versus Control_AC
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
  library(pheatmap)
  library(openxlsx)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
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

primary_dge_file <- file.path(
  objects_directory,
  "GSE287540_primary_DGEList.rds"
)

primary_metadata_file <- file.path(
  objects_directory,
  "GSE287540_primary_metadata.rds"
)

# Check input files
required_files <- c(
  primary_dge_file,
  primary_metadata_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  cat("\nMissing files:\n")
  print(missing_files)
  
  stop(
    "Required primary-analysis files are missing. ",
    "Run 06_review_QC_and_define_samples.R first."
  )
}

# ============================================================
# Load primary analysis data
# ============================================================

cat("\nLoading primary analysis data...\n")

primary_dge <- readRDS(
  primary_dge_file
)

primary_metadata <- readRDS(
  primary_metadata_file
)

# Confirm sample alignment
if (!identical(
  colnames(primary_dge),
  primary_metadata$Sample_ID
)) {
  stop(
    "The expression matrix and metadata are not aligned."
  )
}

rownames(primary_metadata) <-
  primary_metadata$Sample_ID

# Prepare analysis factors
primary_metadata$Primary_group <- factor(
  primary_metadata$Primary_group,
  levels = c(
    "Control_AC",
    "ChILI_TOL"
  )
)

primary_metadata$Batch <- factor(
  primary_metadata$Batch
)

# ============================================================
# Refilter genes for the 59 primary samples
# ============================================================

cat("\nFiltering genes for the primary comparison...\n")

genes_before_primary_filter <- nrow(
  primary_dge
)

keep_primary_gene <- filterByExpr(
  primary_dge,
  group = primary_metadata$Primary_group,
  min.count = 10,
  min.total.count = 20
)

primary_dge <- primary_dge[
  keep_primary_gene,
  ,
  keep.lib.sizes = FALSE
]

primary_dge <- calcNormFactors(
  primary_dge,
  method = "TMM"
)

genes_after_primary_filter <- nrow(
  primary_dge
)

genes_removed_primary <- genes_before_primary_filter -
  genes_after_primary_filter

cat(
  "Genes before primary filtering:",
  genes_before_primary_filter,
  "\n"
)

cat(
  "Genes removed:",
  genes_removed_primary,
  "\n"
)

cat(
  "Genes retained:",
  genes_after_primary_filter,
  "\n"
)

# ============================================================
# Create the statistical design
# ============================================================

# Batch is included to control the strong sequencing batch effect.
design_matrix <- model.matrix(
  ~ Batch + Primary_group,
  data = primary_metadata
)

cat("\nDesign matrix columns:\n")
print(colnames(design_matrix))

# Confirm that the design matrix is full rank
design_rank <- qr(design_matrix)$rank

if (design_rank != ncol(design_matrix)) {
  stop(
    "The differential-expression design matrix is not full rank."
  )
}

# ============================================================
# Estimate dispersion and fit the edgeR model
# ============================================================

cat("\nEstimating gene-wise dispersion...\n")

primary_dge <- estimateDisp(
  primary_dge,
  design = design_matrix,
  robust = TRUE
)

cat("\nFitting the quasi-likelihood model...\n")

ql_fit <- glmQLFit(
  primary_dge,
  design = design_matrix,
  robust = TRUE
)

# Identify the coefficient representing ChILI_TOL versus Control_AC
comparison_coefficient <- grep(
  pattern = "^Primary_groupChILI_TOL$",
  x = colnames(design_matrix),
  value = TRUE
)

if (length(comparison_coefficient) != 1) {
  stop(
    "The ChILI_TOL comparison coefficient could not be identified."
  )
}

ql_test <- glmQLFTest(
  ql_fit,
  coef = comparison_coefficient
)

# ============================================================
# Extract differential-expression results
# ============================================================

de_results <- topTags(
  ql_test,
  n = Inf,
  sort.by = "PValue"
)$table

de_results$Gene_ID <- rownames(
  de_results
)

# Remove Ensembl version numbers if present
de_results$Ensembl_ID <- sub(
  pattern = "\\..*$",
  replacement = "",
  x = de_results$Gene_ID
)

# ============================================================
# Add gene symbols and gene names
# ============================================================

cat("\nAdding human gene annotations...\n")

valid_ensembl_keys <- keys(
  org.Hs.eg.db,
  keytype = "ENSEMBL"
)

is_valid_ensembl <- de_results$Ensembl_ID %in%
  valid_ensembl_keys

de_results$Gene_symbol <- NA_character_
de_results$Gene_name <- NA_character_
de_results$Entrez_ID <- NA_character_

if (any(is_valid_ensembl)) {
  
  valid_ids <- unique(
    de_results$Ensembl_ID[
      is_valid_ensembl
    ]
  )
  
  gene_annotations <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = valid_ids,
    keytype = "ENSEMBL",
    columns = c(
      "SYMBOL",
      "GENENAME",
      "ENTREZID"
    )
  )
  
  gene_annotations <- gene_annotations[
    !duplicated(gene_annotations$ENSEMBL),
    ,
    drop = FALSE
  ]
  
  annotation_match <- match(
    de_results$Ensembl_ID,
    gene_annotations$ENSEMBL
  )
  
  de_results$Gene_symbol <-
    gene_annotations$SYMBOL[
      annotation_match
    ]
  
  de_results$Gene_name <-
    gene_annotations$GENENAME[
      annotation_match
    ]
  
  de_results$Entrez_ID <-
    gene_annotations$ENTREZID[
      annotation_match
    ]
}

# ============================================================
# Classify differentially expressed genes
# ============================================================

fdr_threshold <- 0.05
logfc_threshold <- 1

de_results$Regulation <- "Not significant"

de_results$Regulation[
  de_results$FDR < fdr_threshold &
    de_results$logFC >= logfc_threshold
] <- "Upregulated"

de_results$Regulation[
  de_results$FDR < fdr_threshold &
    de_results$logFC <= -logfc_threshold
] <- "Downregulated"

significant_results <- de_results[
  de_results$Regulation != "Not significant",
  ,
  drop = FALSE
]

upregulated_results <- de_results[
  de_results$Regulation == "Upregulated",
  ,
  drop = FALSE
]

downregulated_results <- de_results[
  de_results$Regulation == "Downregulated",
  ,
  drop = FALSE
]

number_upregulated <- nrow(
  upregulated_results
)

number_downregulated <- nrow(
  downregulated_results
)

number_significant <- nrow(
  significant_results
)

cat("\nDifferential-expression summary:\n")

cat(
  "Upregulated genes:",
  number_upregulated,
  "\n"
)

cat(
  "Downregulated genes:",
  number_downregulated,
  "\n"
)

cat(
  "Total significant genes:",
  number_significant,
  "\n"
)

# ============================================================
# Create the volcano plot
# ============================================================

volcano_data <- de_results

volcano_data$Negative_log10_FDR <- -log10(
  pmax(
    volcano_data$FDR,
    1e-300
  )
)

# Label the 15 most statistically significant genes
label_candidates <- volcano_data[
  volcano_data$Regulation != "Not significant" &
    !is.na(volcano_data$Gene_symbol),
  ,
  drop = FALSE
]

label_candidates <- label_candidates[
  order(label_candidates$FDR),
  ,
  drop = FALSE
]

label_genes <- head(
  label_candidates$Gene_symbol,
  15
)

volcano_data$Plot_label <- ifelse(
  volcano_data$Gene_symbol %in% label_genes,
  volcano_data$Gene_symbol,
  NA
)

volcano_plot <- ggplot(
  volcano_data,
  aes(
    x = logFC,
    y = Negative_log10_FDR,
    color = Regulation
  )
) +
  geom_point(
    alpha = 0.65,
    size = 1.7
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
    yintercept = -log10(fdr_threshold),
    linetype = "dashed",
    color = "grey40"
  ) +
  geom_text_repel(
    aes(label = Plot_label),
    size = 3.2,
    max.overlaps = 20,
    box.padding = 0.4,
    point.padding = 0.3
  ) +
  scale_color_manual(
    values = c(
      "Downregulated" = "#2166AC",
      "Not significant" = "grey75",
      "Upregulated" = "#B2182B"
    )
  ) +
  labs(
    title = "Differential Expression in ICI-Induced Liver Injury",
    subtitle = "ChILI_TOL versus Control_AC, adjusted for sequencing batch",
    x = "log2 fold change",
    y = "-log10 adjusted P-value",
    color = "Regulation"
  ) +
  theme_classic(
    base_size = 12
  )

ggsave(
  filename = file.path(
    figures_directory,
    "Figure_07_volcano_ChILI_TOL_vs_Control_AC.pdf"
  ),
  plot = volcano_plot,
  width = 9,
  height = 7
)

ggsave(
  filename = file.path(
    figures_directory,
    "Figure_07_volcano_ChILI_TOL_vs_Control_AC.png"
  ),
  plot = volcano_plot,
  width = 9,
  height = 7,
  dpi = 300
)

# ============================================================
# Create an MD plot
# ============================================================

pdf(
  file = file.path(
    figures_directory,
    "Figure_08_MD_plot_ChILI_TOL_vs_Control_AC.pdf"
  ),
  width = 8,
  height = 6
)

plotMD(
  ql_test,
  status = de_results$Regulation,
  values = c(
    "Downregulated",
    "Upregulated"
  ),
  col = c(
    "#2166AC",
    "#B2182B"
  ),
  main = "ChILI_TOL versus Control_AC"
)

abline(
  h = c(
    -logfc_threshold,
    logfc_threshold
  ),
  lty = 2,
  col = "grey40"
)

dev.off()

# ============================================================
# Create the significant-gene heatmap
# ============================================================

normalized_log_cpm <- cpm(
  primary_dge,
  normalized.lib.sizes = TRUE,
  log = TRUE,
  prior.count = 2
)

# Select the 40 most significant annotated genes
heatmap_candidates <- de_results[
  !is.na(de_results$Gene_symbol),
  ,
  drop = FALSE
]

heatmap_candidates <- heatmap_candidates[
  order(heatmap_candidates$FDR),
  ,
  drop = FALSE
]

number_of_heatmap_genes <- min(
  40,
  nrow(heatmap_candidates)
)

heatmap_gene_ids <- head(
  heatmap_candidates$Gene_ID,
  number_of_heatmap_genes
)

heatmap_matrix <- normalized_log_cpm[
  heatmap_gene_ids,
  ,
  drop = FALSE
]

heatmap_gene_symbols <- heatmap_candidates$Gene_symbol[
  match(
    heatmap_gene_ids,
    heatmap_candidates$Gene_ID
  )
]

rownames(heatmap_matrix) <- make.unique(
  heatmap_gene_symbols
)

# Remove the batch effect for heatmap visualization only
heatmap_group_design <- model.matrix(
  ~ Primary_group,
  data = primary_metadata
)

heatmap_matrix_corrected <- removeBatchEffect(
  heatmap_matrix,
  batch = primary_metadata$Batch,
  design = heatmap_group_design
)

# Convert each gene to a row Z-score
heatmap_z_scores <- t(
  scale(
    t(heatmap_matrix_corrected)
  )
)

heatmap_annotation <- data.frame(
  Group = primary_metadata$Primary_group,
  Batch = primary_metadata$Batch
)

rownames(heatmap_annotation) <-
  primary_metadata$Sample_ID

heatmap_colors <- colorRampPalette(
  c(
    "#2166AC",
    "white",
    "#B2182B"
  )
)(100)

pdf(
  file = file.path(
    figures_directory,
    "Figure_09_top40_DEG_heatmap.pdf"
  ),
  width = 11,
  height = 10
)

pheatmap(
  heatmap_z_scores,
  annotation_col = heatmap_annotation,
  show_colnames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  color = heatmap_colors,
  main = "Top Differentially Expressed Genes"
)

dev.off()

png(
  filename = file.path(
    figures_directory,
    "Figure_09_top40_DEG_heatmap.png"
  ),
  width = 3300,
  height = 3000,
  res = 300
)

pheatmap(
  heatmap_z_scores,
  annotation_col = heatmap_annotation,
  show_colnames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  color = heatmap_colors,
  main = "Top Differentially Expressed Genes"
)

dev.off()

# ============================================================
# Save differential-expression tables
# ============================================================

write.csv(
  de_results,
  file = file.path(
    processed_data_directory,
    "GSE287540_all_differential_expression_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  significant_results,
  file = file.path(
    processed_data_directory,
    "GSE287540_significant_DEGs.csv"
  ),
  row.names = FALSE
)

write.csv(
  upregulated_results,
  file = file.path(
    processed_data_directory,
    "GSE287540_upregulated_DEGs.csv"
  ),
  row.names = FALSE
)

write.csv(
  downregulated_results,
  file = file.path(
    processed_data_directory,
    "GSE287540_downregulated_DEGs.csv"
  ),
  row.names = FALSE
)

# ============================================================
# Create the Excel results workbook
# ============================================================

analysis_summary <- data.frame(
  Measurement = c(
    "Comparison",
    "Control samples",
    "ChILI samples",
    "Genes tested",
    "FDR threshold",
    "Absolute log2 fold-change threshold",
    "Upregulated genes",
    "Downregulated genes",
    "Total significant genes",
    "Batch adjustment"
  ),
  
  Value = c(
    "ChILI_TOL versus Control_AC",
    sum(
      primary_metadata$Primary_group ==
        "Control_AC"
    ),
    sum(
      primary_metadata$Primary_group ==
        "ChILI_TOL"
    ),
    nrow(de_results),
    fdr_threshold,
    logfc_threshold,
    number_upregulated,
    number_downregulated,
    number_significant,
    "Included in edgeR design"
  ),
  
  stringsAsFactors = FALSE
)

de_workbook <- createWorkbook()

addWorksheet(
  de_workbook,
  "Analysis summary"
)

writeData(
  de_workbook,
  sheet = "Analysis summary",
  x = analysis_summary
)

addWorksheet(
  de_workbook,
  "All tested genes"
)

writeData(
  de_workbook,
  sheet = "All tested genes",
  x = de_results
)

addWorksheet(
  de_workbook,
  "Significant DEGs"
)

writeData(
  de_workbook,
  sheet = "Significant DEGs",
  x = significant_results
)

addWorksheet(
  de_workbook,
  "Upregulated"
)

writeData(
  de_workbook,
  sheet = "Upregulated",
  x = upregulated_results
)

addWorksheet(
  de_workbook,
  "Downregulated"
)

writeData(
  de_workbook,
  sheet = "Downregulated",
  x = downregulated_results
)

for (sheet_name in names(de_workbook)) {
  
  setColWidths(
    de_workbook,
    sheet = sheet_name,
    cols = 1:12,
    widths = "auto"
  )
}

saveWorkbook(
  de_workbook,
  file = file.path(
    tables_directory,
    "GSE287540_differential_expression_results.xlsx"
  ),
  overwrite = TRUE
)

# ============================================================
# Save statistical R objects
# ============================================================

saveRDS(
  primary_dge,
  file = file.path(
    objects_directory,
    "GSE287540_primary_DGEList_refiltered.rds"
  )
)

saveRDS(
  design_matrix,
  file = file.path(
    objects_directory,
    "GSE287540_DE_design_matrix.rds"
  )
)

saveRDS(
  ql_fit,
  file = file.path(
    objects_directory,
    "GSE287540_edgeR_QL_fit.rds"
  )
)

saveRDS(
  ql_test,
  file = file.path(
    objects_directory,
    "GSE287540_edgeR_QL_test.rds"
  )
)

saveRDS(
  de_results,
  file = file.path(
    objects_directory,
    "GSE287540_DE_results.rds"
  )
)

# ============================================================
# Final summary
# ============================================================

cat("\n============================================\n")
cat("Differential-expression analysis completed\n")
cat("============================================\n")

cat("\nComparison: ChILI_TOL versus Control_AC\n")

cat(
  "Genes tested:",
  nrow(de_results),
  "\n"
)

cat(
  "Upregulated genes:",
  number_upregulated,
  "\n"
)

cat(
  "Downregulated genes:",
  number_downregulated,
  "\n"
)

cat(
  "Total significant genes:",
  number_significant,
  "\n"
)

cat("\nThresholds:\n")
cat("FDR < 0.05\n")
cat("|log2 fold change| >= 1\n")

cat("\nMain files created:\n")
cat("results/tables/GSE287540_differential_expression_results.xlsx\n")
cat("results/figures/Figure_07_volcano_ChILI_TOL_vs_Control_AC.png\n")
cat("results/figures/Figure_09_top40_DEG_heatmap.png\n")
cat("results/R_objects/GSE287540_DE_results.rds\n")