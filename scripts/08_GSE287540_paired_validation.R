# ============================================================
# Script 08: Paired Validation Analysis of GSE287540
# Comparison: ChILI_TOL versus ChILI_BC
# Project: ICI-Induced Liver Injury and miRNA Analysis
# ============================================================

# Clear the R environment
rm(list = ls())

# Load required packages
suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
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

full_dge_file <- file.path(
  objects_directory,
  "GSE287540_DGEList_filtered_TMM.rds"
)

metadata_file <- file.path(
  objects_directory,
  "GSE287540_metadata_analysis_ready.rds"
)

primary_results_file <- file.path(
  objects_directory,
  "GSE287540_DE_results.rds"
)

# Check required files
required_files <- c(
  full_dge_file,
  metadata_file,
  primary_results_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  cat("\nMissing input files:\n")
  print(missing_files)
  
  stop(
    "Required files are missing. ",
    "Run the previous scripts first."
  )
}

# ============================================================
# Load data
# ============================================================

cat("\nLoading expression data and metadata...\n")

full_dge <- readRDS(
  full_dge_file
)

metadata <- readRDS(
  metadata_file
)

primary_results <- readRDS(
  primary_results_file
)

# Validate sample alignment
if (!identical(
  colnames(full_dge),
  metadata$Sample_ID
)) {
  stop(
    "The expression matrix and metadata are not aligned."
  )
}

# ============================================================
# Identify paired ChILI subjects
# ============================================================

chili_candidate <- metadata[
  metadata$Group == "chili" &
    metadata$Time_point %in% c("BC", "TOL"),
  ,
  drop = FALSE
]

subject_time_table <- table(
  chili_candidate$Subject_ID,
  chili_candidate$Time_point
)

# Retain subjects with both BC and TOL samples
paired_subjects <- rownames(
  subject_time_table
)[
  subject_time_table[, "BC"] >= 1 &
    subject_time_table[, "TOL"] >= 1
]

number_of_paired_subjects <- length(
  paired_subjects
)

cat(
  "\nNumber of paired ChILI subjects:",
  number_of_paired_subjects,
  "\n"
)

if (number_of_paired_subjects < 3) {
  stop(
    "Too few paired subjects were detected for paired analysis."
  )
}

# Select paired samples
paired_metadata <- metadata[
  metadata$Subject_ID %in% paired_subjects &
    metadata$Group == "chili" &
    metadata$Time_point %in% c("BC", "TOL"),
  ,
  drop = FALSE
]

# Ensure one BC and one TOL sample per subject
paired_metadata <- paired_metadata[
  order(
    paired_metadata$Subject_ID,
    paired_metadata$Time_point
  ),
  ,
  drop = FALSE
]

# Check duplicated subject-time combinations
subject_time_key <- paste(
  paired_metadata$Subject_ID,
  paired_metadata$Time_point,
  sep = "_"
)

if (anyDuplicated(subject_time_key) > 0) {
  stop(
    "Duplicated subject-time combinations were detected."
  )
}

# Prepare factors
paired_metadata$Subject_ID <- factor(
  paired_metadata$Subject_ID
)

paired_metadata$Paired_time <- factor(
  paired_metadata$Time_point,
  levels = c("BC", "TOL")
)

rownames(paired_metadata) <-
  paired_metadata$Sample_ID

cat("\nPaired samples:\n")

print(
  paired_metadata[
    ,
    c(
      "Sample_ID",
      "Subject_ID",
      "Time_point",
      "Batch",
      "Treatment"
    )
  ]
)

# ============================================================
# Select paired expression samples
# ============================================================

paired_dge <- full_dge[
  ,
  paired_metadata$Sample_ID,
  keep.lib.sizes = FALSE
]

# Confirm sample alignment
if (!identical(
  colnames(paired_dge),
  paired_metadata$Sample_ID
)) {
  stop(
    "The paired expression matrix and metadata are not aligned."
  )
}

# ============================================================
# Filter genes for the paired analysis
# ============================================================

genes_before_filtering <- nrow(
  paired_dge
)

keep_paired_gene <- filterByExpr(
  paired_dge,
  group = paired_metadata$Paired_time,
  min.count = 10,
  min.total.count = 15
)

paired_dge <- paired_dge[
  keep_paired_gene,
  ,
  keep.lib.sizes = FALSE
]

paired_dge <- calcNormFactors(
  paired_dge,
  method = "TMM"
)

genes_after_filtering <- nrow(
  paired_dge
)

cat(
  "\nGenes before paired filtering:",
  genes_before_filtering,
  "\n"
)

cat(
  "Genes retained for paired analysis:",
  genes_after_filtering,
  "\n"
)

# ============================================================
# Create the paired statistical design
# ============================================================

# Subject_ID controls individual patient differences.
# Batch is not added because paired samples from the same
# patient are controlled through the subject fixed effect.

paired_design <- model.matrix(
  ~ Subject_ID + Paired_time,
  data = paired_metadata
)

cat("\nPaired design columns:\n")
print(colnames(paired_design))

paired_design_rank <- qr(
  paired_design
)$rank

if (paired_design_rank != ncol(paired_design)) {
  stop(
    "The paired statistical design is not full rank."
  )
}

# ============================================================
# Fit the paired edgeR model
# ============================================================

cat("\nEstimating dispersion for paired analysis...\n")

paired_dge <- estimateDisp(
  paired_dge,
  design = paired_design,
  robust = TRUE
)

cat("\nFitting the paired quasi-likelihood model...\n")

paired_fit <- glmQLFit(
  paired_dge,
  design = paired_design,
  robust = TRUE
)

paired_coefficient <- grep(
  pattern = "^Paired_timeTOL$",
  x = colnames(paired_design),
  value = TRUE
)

if (length(paired_coefficient) != 1) {
  stop(
    "The TOL versus BC coefficient could not be identified."
  )
}

paired_test <- glmQLFTest(
  paired_fit,
  coef = paired_coefficient
)

paired_results <- topTags(
  paired_test,
  n = Inf,
  sort.by = "PValue"
)$table

paired_results$Gene_ID <- rownames(
  paired_results
)

paired_results$Ensembl_ID <- sub(
  pattern = "\\..*$",
  replacement = "",
  x = paired_results$Gene_ID
)

# ============================================================
# Add human gene annotations
# ============================================================

cat("\nAdding gene annotations...\n")

valid_ensembl_keys <- keys(
  org.Hs.eg.db,
  keytype = "ENSEMBL"
)

valid_gene <- paired_results$Ensembl_ID %in%
  valid_ensembl_keys

paired_results$Gene_symbol <- NA_character_
paired_results$Gene_name <- NA_character_
paired_results$Entrez_ID <- NA_character_

if (any(valid_gene)) {
  
  valid_ids <- unique(
    paired_results$Ensembl_ID[
      valid_gene
    ]
  )
  
  annotations <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = valid_ids,
    keytype = "ENSEMBL",
    columns = c(
      "SYMBOL",
      "GENENAME",
      "ENTREZID"
    )
  )
  
  annotations <- annotations[
    !duplicated(annotations$ENSEMBL),
    ,
    drop = FALSE
  ]
  
  annotation_match <- match(
    paired_results$Ensembl_ID,
    annotations$ENSEMBL
  )
  
  paired_results$Gene_symbol <-
    annotations$SYMBOL[
      annotation_match
    ]
  
  paired_results$Gene_name <-
    annotations$GENENAME[
      annotation_match
    ]
  
  paired_results$Entrez_ID <-
    annotations$ENTREZID[
      annotation_match
    ]
}

# ============================================================
# Classify paired differential-expression results
# ============================================================

fdr_threshold <- 0.05
logfc_threshold <- 1

paired_results$Regulation <- "Not significant"

paired_results$Regulation[
  paired_results$FDR < fdr_threshold &
    paired_results$logFC >= logfc_threshold
] <- "Upregulated"

paired_results$Regulation[
  paired_results$FDR < fdr_threshold &
    paired_results$logFC <= -logfc_threshold
] <- "Downregulated"

paired_significant <- paired_results[
  paired_results$Regulation != "Not significant",
  ,
  drop = FALSE
]

paired_upregulated <- paired_results[
  paired_results$Regulation == "Upregulated",
  ,
  drop = FALSE
]

paired_downregulated <- paired_results[
  paired_results$Regulation == "Downregulated",
  ,
  drop = FALSE
]

# ============================================================
# Compare paired results with the primary results
# ============================================================

primary_validation_columns <- primary_results[
  ,
  c(
    "Gene_ID",
    "Gene_symbol",
    "logFC",
    "PValue",
    "FDR",
    "Regulation"
  )
]

colnames(primary_validation_columns) <- c(
  "Gene_ID",
  "Gene_symbol_primary",
  "Primary_logFC",
  "Primary_PValue",
  "Primary_FDR",
  "Primary_regulation"
)

paired_validation_columns <- paired_results[
  ,
  c(
    "Gene_ID",
    "Gene_symbol",
    "logFC",
    "PValue",
    "FDR",
    "Regulation"
  )
]

colnames(paired_validation_columns) <- c(
  "Gene_ID",
  "Gene_symbol_paired",
  "Paired_logFC",
  "Paired_PValue",
  "Paired_FDR",
  "Paired_regulation"
)

validation_table <- merge(
  primary_validation_columns,
  paired_validation_columns,
  by = "Gene_ID",
  all = FALSE
)

validation_table$Gene_symbol <- ifelse(
  !is.na(validation_table$Gene_symbol_primary),
  validation_table$Gene_symbol_primary,
  validation_table$Gene_symbol_paired
)

validation_table$Direction_concordant <- sign(
  validation_table$Primary_logFC
) == sign(
  validation_table$Paired_logFC
)

# Primary significant genes with supporting paired evidence
validation_table$Validation_status <- "Not validated"

validation_table$Validation_status[
  validation_table$Primary_regulation !=
    "Not significant" &
    validation_table$Direction_concordant &
    validation_table$Paired_PValue < 0.05
] <- "Direction and nominal P-value validated"

validation_table$Validation_status[
  validation_table$Primary_regulation !=
    "Not significant" &
    validation_table$Direction_concordant &
    validation_table$Paired_FDR < 0.05
] <- "FDR validated"

validated_primary_genes <- validation_table[
  validation_table$Validation_status !=
    "Not validated",
  ,
  drop = FALSE
]

validated_primary_genes <- validated_primary_genes[
  order(
    validated_primary_genes$Paired_PValue
  ),
  ,
  drop = FALSE
]

# ============================================================
# Calculate direction concordance
# ============================================================

primary_significant_validation <- validation_table[
  validation_table$Primary_regulation !=
    "Not significant",
  ,
  drop = FALSE
]

direction_concordance_percent <- round(
  100 * mean(
    primary_significant_validation$Direction_concordant,
    na.rm = TRUE
  ),
  1
)

number_nominally_validated <- sum(
  validation_table$Validation_status ==
    "Direction and nominal P-value validated"
)

number_fdr_validated <- sum(
  validation_table$Validation_status ==
    "FDR validated"
)

# ============================================================
# Volcano plot for paired analysis
# ============================================================

paired_volcano <- paired_results

paired_volcano$Negative_log10_FDR <- -log10(
  pmax(
    paired_volcano$FDR,
    1e-300
  )
)

paired_label_candidates <- paired_volcano[
  paired_volcano$Regulation != "Not significant" &
    !is.na(paired_volcano$Gene_symbol),
  ,
  drop = FALSE
]

paired_label_candidates <- paired_label_candidates[
  order(paired_label_candidates$FDR),
  ,
  drop = FALSE
]

paired_label_genes <- head(
  paired_label_candidates$Gene_symbol,
  15
)

paired_volcano$Plot_label <- ifelse(
  paired_volcano$Gene_symbol %in%
    paired_label_genes,
  paired_volcano$Gene_symbol,
  NA
)

paired_volcano_plot <- ggplot(
  paired_volcano,
  aes(
    x = logFC,
    y = Negative_log10_FDR,
    color = Regulation
  )
) +
  geom_point(
    size = 1.7,
    alpha = 0.65
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
    max.overlaps = 20
  ) +
  scale_color_manual(
    values = c(
      "Downregulated" = "#2166AC",
      "Not significant" = "grey75",
      "Upregulated" = "#B2182B"
    )
  ) +
  labs(
    title = "Paired Differential-Expression Analysis",
    subtitle = "ChILI_TOL versus ChILI_BC",
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
    "Figure_10_paired_volcano_TOL_vs_BC.pdf"
  ),
  plot = paired_volcano_plot,
  width = 9,
  height = 7
)

ggsave(
  filename = file.path(
    figures_directory,
    "Figure_10_paired_volcano_TOL_vs_BC.png"
  ),
  plot = paired_volcano_plot,
  width = 9,
  height = 7,
  dpi = 300
)

# ============================================================
# Fold-change concordance plot
# ============================================================

concordance_plot_data <-
  primary_significant_validation[
    !is.na(
      primary_significant_validation$Gene_symbol
    ),
    ,
    drop = FALSE
  ]

fold_change_correlation <- suppressWarnings(
  cor(
    concordance_plot_data$Primary_logFC,
    concordance_plot_data$Paired_logFC,
    method = "spearman",
    use = "complete.obs"
  )
)

concordance_plot <- ggplot(
  concordance_plot_data,
  aes(
    x = Primary_logFC,
    y = Paired_logFC,
    color = Direction_concordant
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey60"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey60"
  ) +
  geom_point(
    size = 2.5,
    alpha = 0.8
  ) +
  scale_color_manual(
    values = c(
      "FALSE" = "grey65",
      "TRUE" = "#1B9E77"
    )
  ) +
  labs(
    title = "Fold-Change Concordance",
    subtitle = paste0(
      "Primary DEGs; Spearman rho = ",
      round(fold_change_correlation, 3)
    ),
    x = "Primary log2 fold change",
    y = "Paired log2 fold change",
    color = "Same direction"
  ) +
  theme_classic(
    base_size = 12
  )

ggsave(
  filename = file.path(
    figures_directory,
    "Figure_11_primary_paired_concordance.pdf"
  ),
  plot = concordance_plot,
  width = 8,
  height = 7
)

ggsave(
  filename = file.path(
    figures_directory,
    "Figure_11_primary_paired_concordance.png"
  ),
  plot = concordance_plot,
  width = 8,
  height = 7,
  dpi = 300
)

# ============================================================
# Save output tables
# ============================================================

write.csv(
  paired_results,
  file = file.path(
    processed_data_directory,
    "GSE287540_paired_all_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  paired_significant,
  file = file.path(
    processed_data_directory,
    "GSE287540_paired_significant_DEGs.csv"
  ),
  row.names = FALSE
)

write.csv(
  validation_table,
  file = file.path(
    processed_data_directory,
    "GSE287540_primary_paired_validation_table.csv"
  ),
  row.names = FALSE
)

write.csv(
  validated_primary_genes,
  file = file.path(
    processed_data_directory,
    "GSE287540_validated_primary_genes.csv"
  ),
  row.names = FALSE
)

# ============================================================
# Create Excel workbook
# ============================================================

validation_summary <- data.frame(
  Measurement = c(
    "Paired subjects",
    "Paired samples",
    "Genes tested",
    "Paired upregulated genes",
    "Paired downregulated genes",
    "Paired significant genes",
    "Primary DEG direction concordance (%)",
    "Nominally validated primary genes",
    "FDR-validated primary genes",
    "Fold-change Spearman correlation"
  ),
  
  Value = c(
    number_of_paired_subjects,
    nrow(paired_metadata),
    nrow(paired_results),
    nrow(paired_upregulated),
    nrow(paired_downregulated),
    nrow(paired_significant),
    direction_concordance_percent,
    number_nominally_validated,
    number_fdr_validated,
    round(fold_change_correlation, 3)
  ),
  
  stringsAsFactors = FALSE
)

validation_workbook <- createWorkbook()

addWorksheet(
  validation_workbook,
  "Validation summary"
)

writeData(
  validation_workbook,
  sheet = "Validation summary",
  x = validation_summary
)

addWorksheet(
  validation_workbook,
  "Paired samples"
)

writeData(
  validation_workbook,
  sheet = "Paired samples",
  x = paired_metadata
)

addWorksheet(
  validation_workbook,
  "All paired results"
)

writeData(
  validation_workbook,
  sheet = "All paired results",
  x = paired_results
)

addWorksheet(
  validation_workbook,
  "Paired significant DEGs"
)

writeData(
  validation_workbook,
  sheet = "Paired significant DEGs",
  x = paired_significant
)

addWorksheet(
  validation_workbook,
  "Primary validation"
)

writeData(
  validation_workbook,
  sheet = "Primary validation",
  x = validation_table
)

addWorksheet(
  validation_workbook,
  "Validated primary genes"
)

writeData(
  validation_workbook,
  sheet = "Validated primary genes",
  x = validated_primary_genes
)

for (sheet_name in names(validation_workbook)) {
  
  setColWidths(
    validation_workbook,
    sheet = sheet_name,
    cols = 1:15,
    widths = "auto"
  )
}

saveWorkbook(
  validation_workbook,
  file = file.path(
    tables_directory,
    "GSE287540_paired_validation_results.xlsx"
  ),
  overwrite = TRUE
)

# ============================================================
# Save R objects
# ============================================================

saveRDS(
  paired_metadata,
  file = file.path(
    objects_directory,
    "GSE287540_paired_metadata.rds"
  )
)

saveRDS(
  paired_results,
  file = file.path(
    objects_directory,
    "GSE287540_paired_DE_results.rds"
  )
)

saveRDS(
  validation_table,
  file = file.path(
    objects_directory,
    "GSE287540_primary_paired_validation.rds"
  )
)

# ============================================================
# Final summary
# ============================================================

cat("\n============================================\n")
cat("Paired validation analysis completed\n")
cat("============================================\n")

cat(
  "\nPaired subjects:",
  number_of_paired_subjects,
  "\n"
)

cat(
  "Paired samples:",
  nrow(paired_metadata),
  "\n"
)

cat(
  "Genes tested:",
  nrow(paired_results),
  "\n"
)

cat(
  "Paired upregulated genes:",
  nrow(paired_upregulated),
  "\n"
)

cat(
  "Paired downregulated genes:",
  nrow(paired_downregulated),
  "\n"
)

cat(
  "Paired significant genes:",
  nrow(paired_significant),
  "\n"
)

cat(
  "Primary DEG direction concordance:",
  direction_concordance_percent,
  "%\n"
)

cat(
  "Nominally validated primary genes:",
  number_nominally_validated,
  "\n"
)

cat(
  "FDR-validated primary genes:",
  number_fdr_validated,
  "\n"
)

cat(
  "Fold-change Spearman correlation:",
  round(fold_change_correlation, 3),
  "\n"
)

cat("\nMain files created:\n")
cat("results/tables/GSE287540_paired_validation_results.xlsx\n")
cat("results/figures/Figure_10_paired_volcano_TOL_vs_BC.png\n")
cat("results/figures/Figure_11_primary_paired_concordance.png\n")