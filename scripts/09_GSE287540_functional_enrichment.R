# ============================================================
# Script 09: Functional Enrichment Analysis of GSE287540
# Project: ICI-Induced Liver Injury and miRNA Analysis
# ============================================================

# Clear the R environment
rm(list = ls())

options(
  stringsAsFactors = FALSE,
  timeout = 1800
)

# ============================================================
# Install missing packages
# ============================================================

cran_packages <- c(
  "ggplot2",
  "openxlsx",
  "msigdbr"
)

bioconductor_packages <- c(
  "clusterProfiler",
  "enrichplot",
  "DOSE",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "fgsea"
)

missing_cran <- cran_packages[
  !vapply(
    cran_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_cran) > 0) {
  
  install.packages(
    missing_cran,
    repos = "https://cloud.r-project.org"
  )
}

if (!requireNamespace(
  "BiocManager",
  quietly = TRUE
)) {
  
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}

missing_bioconductor <- bioconductor_packages[
  !vapply(
    bioconductor_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_bioconductor) > 0) {
  
  BiocManager::install(
    missing_bioconductor,
    ask = FALSE,
    update = FALSE
  )
}

# ============================================================
# Load required packages
# ============================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
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

primary_results_file <- file.path(
  objects_directory,
  "GSE287540_DE_results.rds"
)

validation_results_file <- file.path(
  objects_directory,
  "GSE287540_primary_paired_validation.rds"
)

required_files <- c(
  primary_results_file,
  validation_results_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  cat("\nMissing files:\n")
  print(missing_files)
  
  stop(
    "Required differential-expression files are missing."
  )
}

# ============================================================
# Load differential-expression results
# ============================================================

cat("\nLoading differential-expression results...\n")

primary_results <- readRDS(
  primary_results_file
)

validation_results <- readRDS(
  validation_results_file
)

required_primary_columns <- c(
  "Gene_ID",
  "Gene_symbol",
  "Entrez_ID",
  "logFC",
  "PValue",
  "FDR",
  "Regulation"
)

missing_primary_columns <- setdiff(
  required_primary_columns,
  colnames(primary_results)
)

if (length(missing_primary_columns) > 0) {
  
  stop(
    "Missing columns in primary results: ",
    paste(
      missing_primary_columns,
      collapse = ", "
    )
  )
}

# ============================================================
# Define gene sets
# ============================================================

primary_significant <- primary_results[
  primary_results$Regulation != "Not significant" &
    !is.na(primary_results$Entrez_ID),
  ,
  drop = FALSE
]

primary_upregulated <- primary_results[
  primary_results$Regulation == "Upregulated" &
    !is.na(primary_results$Entrez_ID),
  ,
  drop = FALSE
]

primary_downregulated <- primary_results[
  primary_results$Regulation == "Downregulated" &
    !is.na(primary_results$Entrez_ID),
  ,
  drop = FALSE
]

validated_gene_ids <- validation_results$Gene_ID[
  validation_results$Validation_status != "Not validated"
]

fdr_validated_gene_ids <- validation_results$Gene_ID[
  validation_results$Validation_status == "FDR validated"
]

validated_genes <- primary_results[
  primary_results$Gene_ID %in% validated_gene_ids &
    !is.na(primary_results$Entrez_ID),
  ,
  drop = FALSE
]

fdr_validated_genes <- primary_results[
  primary_results$Gene_ID %in% fdr_validated_gene_ids &
    !is.na(primary_results$Entrez_ID),
  ,
  drop = FALSE
]

# Convert Entrez identifiers to unique character vectors
all_tested_entrez <- unique(
  as.character(
    primary_results$Entrez_ID[
      !is.na(primary_results$Entrez_ID)
    ]
  )
)

significant_entrez <- unique(
  as.character(
    primary_significant$Entrez_ID
  )
)

upregulated_entrez <- unique(
  as.character(
    primary_upregulated$Entrez_ID
  )
)

downregulated_entrez <- unique(
  as.character(
    primary_downregulated$Entrez_ID
  )
)

validated_entrez <- unique(
  as.character(
    validated_genes$Entrez_ID
  )
)

fdr_validated_entrez <- unique(
  as.character(
    fdr_validated_genes$Entrez_ID
  )
)

cat("\nGene sets prepared:\n")
cat("Primary significant genes:", length(significant_entrez), "\n")
cat("Upregulated genes:", length(upregulated_entrez), "\n")
cat("Downregulated genes:", length(downregulated_entrez), "\n")
cat("Nominally validated genes:", length(validated_entrez), "\n")
cat("FDR-validated genes:", length(fdr_validated_entrez), "\n")

# ============================================================
# Over-representation analysis function
# ============================================================

run_enrichment <- function(
    gene_list,
    gene_universe,
    analysis_name) {
  
  empty_result <- list(
    GO = data.frame(),
    KEGG = data.frame()
  )
  
  if (length(gene_list) < 5) {
    
    warning(
      paste(
        "Too few genes for enrichment:",
        analysis_name
      )
    )
    
    return(empty_result)
  }
  
  cat(
    "\nRunning enrichment for:",
    analysis_name,
    "\n"
  )
  
  go_object <- tryCatch(
    
    enrichGO(
      gene = gene_list,
      universe = gene_universe,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.20,
      readable = TRUE
    ),
    
    error = function(error_message) {
      
      warning(
        paste(
          "GO analysis failed for",
          analysis_name,
          ":",
          error_message$message
        )
      )
      
      return(NULL)
    }
  )
  
  kegg_object <- tryCatch(
    
    enrichKEGG(
      gene = gene_list,
      universe = gene_universe,
      organism = "hsa",
      keyType = "kegg",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.20
    ),
    
    error = function(error_message) {
      
      warning(
        paste(
          "KEGG analysis failed for",
          analysis_name,
          ":",
          error_message$message
        )
      )
      
      return(NULL)
    }
  )
  
  if (!is.null(go_object)) {
    
    go_table <- as.data.frame(
      go_object
    )
    
  } else {
    
    go_table <- data.frame()
  }
  
  if (!is.null(kegg_object)) {
    
    kegg_table <- as.data.frame(
      kegg_object
    )
    
  } else {
    
    kegg_table <- data.frame()
  }
  
  if (nrow(go_table) > 0) {
    go_table$Analysis <- analysis_name
  }
  
  if (nrow(kegg_table) > 0) {
    kegg_table$Analysis <- analysis_name
  }
  
  return(
    list(
      GO = go_table,
      KEGG = kegg_table
    )
  )
}

# ============================================================
# Run GO and KEGG analyses
# ============================================================

enrichment_all <- run_enrichment(
  gene_list = significant_entrez,
  gene_universe = all_tested_entrez,
  analysis_name = "All significant DEGs"
)

enrichment_up <- run_enrichment(
  gene_list = upregulated_entrez,
  gene_universe = all_tested_entrez,
  analysis_name = "Upregulated DEGs"
)

enrichment_down <- run_enrichment(
  gene_list = downregulated_entrez,
  gene_universe = all_tested_entrez,
  analysis_name = "Downregulated DEGs"
)

enrichment_validated <- run_enrichment(
  gene_list = validated_entrez,
  gene_universe = all_tested_entrez,
  analysis_name = "Paired-validated DEGs"
)

enrichment_fdr_validated <- run_enrichment(
  gene_list = fdr_validated_entrez,
  gene_universe = all_tested_entrez,
  analysis_name = "FDR-validated DEGs"
)

# ============================================================
# Combine non-empty tables
# ============================================================

combine_nonempty_tables <- function(table_list) {
  
  table_sizes <- vapply(
    table_list,
    nrow,
    integer(1)
  )
  
  nonempty_tables <- table_list[
    table_sizes > 0
  ]
  
  if (length(nonempty_tables) == 0) {
    return(data.frame())
  }
  
  combined_table <- do.call(
    rbind,
    nonempty_tables
  )
  
  rownames(combined_table) <- NULL
  
  return(combined_table)
}

combined_go_results <- combine_nonempty_tables(
  list(
    enrichment_all$GO,
    enrichment_up$GO,
    enrichment_down$GO,
    enrichment_validated$GO,
    enrichment_fdr_validated$GO
  )
)

combined_kegg_results <- combine_nonempty_tables(
  list(
    enrichment_all$KEGG,
    enrichment_up$KEGG,
    enrichment_down$KEGG,
    enrichment_validated$KEGG,
    enrichment_fdr_validated$KEGG
  )
)

# ============================================================
# Prepare ranked genes for Hallmark GSEA
# ============================================================

cat("\nPreparing ranked genes for Hallmark GSEA...\n")

ranking_table <- primary_results[
  !is.na(primary_results$Gene_symbol) &
    primary_results$Gene_symbol != "",
  ,
  drop = FALSE
]

# Keep one result per gene symbol
ranking_table <- ranking_table[
  order(ranking_table$PValue),
  ,
  drop = FALSE
]

ranking_table <- ranking_table[
  !duplicated(ranking_table$Gene_symbol),
  ,
  drop = FALSE
]

ranking_score <- sign(
  ranking_table$logFC
) * -log10(
  pmax(
    ranking_table$PValue,
    1e-300
  )
)

names(ranking_score) <-
  ranking_table$Gene_symbol

ranking_score <- ranking_score[
  is.finite(ranking_score)
]

ranking_score <- sort(
  ranking_score,
  decreasing = TRUE
)

# ============================================================
# Retrieve Hallmark gene sets
# ============================================================

hallmark_database <- tryCatch(
  
  msigdbr(
    species = "Homo sapiens",
    collection = "H"
  ),
  
  error = function(error_message) {
    
    msigdbr(
      species = "Homo sapiens",
      category = "H"
    )
  }
)

required_hallmark_columns <- c(
  "gs_name",
  "gene_symbol"
)

if (!all(
  required_hallmark_columns %in%
  colnames(hallmark_database)
)) {
  
  stop(
    "Required Hallmark database columns were not detected."
  )
}

hallmark_pathways <- split(
  hallmark_database$gene_symbol,
  hallmark_database$gs_name
)

# ============================================================
# Run Hallmark GSEA
# ============================================================

cat("\nRunning Hallmark GSEA...\n")

hallmark_gsea <- fgseaMultilevel(
  pathways = hallmark_pathways,
  stats = ranking_score,
  minSize = 15,
  maxSize = 500
)

hallmark_gsea <- as.data.frame(
  hallmark_gsea
)

hallmark_gsea <- hallmark_gsea[
  order(hallmark_gsea$padj),
  ,
  drop = FALSE
]

hallmark_gsea$Leading_edge_genes <- vapply(
  hallmark_gsea$leadingEdge,
  FUN = paste,
  collapse = "; ",
  FUN.VALUE = character(1)
)

hallmark_gsea$leadingEdge <- NULL

# ============================================================
# Save CSV results
# ============================================================

write.csv(
  combined_go_results,
  file = file.path(
    processed_data_directory,
    "GSE287540_GO_enrichment_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  combined_kegg_results,
  file = file.path(
    processed_data_directory,
    "GSE287540_KEGG_enrichment_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  hallmark_gsea,
  file = file.path(
    processed_data_directory,
    "GSE287540_Hallmark_GSEA_results.csv"
  ),
  row.names = FALSE
)

# ============================================================
# Enrichment dot-plot function
# ============================================================

create_enrichment_plot <- function(
    enrichment_table,
    analysis_name,
    plot_title,
    output_name,
    maximum_terms = 15) {
  
  if (
    is.null(enrichment_table) ||
    nrow(enrichment_table) == 0
  ) {
    
    cat(
      "\nNo enrichment results for:",
      plot_title,
      "\n"
    )
    
    return(NULL)
  }
  
  selected_results <- enrichment_table[
    enrichment_table$Analysis == analysis_name,
    ,
    drop = FALSE
  ]
  
  if (nrow(selected_results) == 0) {
    
    cat(
      "\nNo significant terms for:",
      plot_title,
      "\n"
    )
    
    return(NULL)
  }
  
  selected_results <- selected_results[
    order(selected_results$p.adjust),
    ,
    drop = FALSE
  ]
  
  selected_results <- head(
    selected_results,
    maximum_terms
  )
  
  selected_results$Description_short <- ifelse(
    nchar(selected_results$Description) > 65,
    paste0(
      substr(
        selected_results$Description,
        1,
        62
      ),
      "..."
    ),
    selected_results$Description
  )
  
  selected_results$Negative_log10_FDR <- -log10(
    pmax(
      selected_results$p.adjust,
      1e-300
    )
  )
  
  enrichment_plot <- ggplot(
    selected_results,
    aes(
      x = Negative_log10_FDR,
      y = reorder(
        Description_short,
        Negative_log10_FDR
      ),
      size = Count,
      color = p.adjust
    )
  ) +
    geom_point(
      alpha = 0.85
    ) +
    scale_color_gradient(
      low = "#B2182B",
      high = "#2166AC"
    ) +
    labs(
      title = plot_title,
      x = "-log10 adjusted P-value",
      y = NULL,
      size = "Gene count",
      color = "Adjusted P-value"
    ) +
    theme_classic(
      base_size = 11
    )
  
  ggsave(
    filename = file.path(
      figures_directory,
      paste0(output_name, ".pdf")
    ),
    plot = enrichment_plot,
    width = 10,
    height = 7
  )
  
  ggsave(
    filename = file.path(
      figures_directory,
      paste0(output_name, ".png")
    ),
    plot = enrichment_plot,
    width = 10,
    height = 7,
    dpi = 300
  )
  
  return(enrichment_plot)
}

# ============================================================
# Create GO enrichment figures
# ============================================================

create_enrichment_plot(
  enrichment_table = combined_go_results,
  analysis_name = "All significant DEGs",
  plot_title = "GO Biological Processes: All Significant DEGs",
  output_name = "Figure_12_GO_all_DEGs"
)

create_enrichment_plot(
  enrichment_table = combined_go_results,
  analysis_name = "Upregulated DEGs",
  plot_title = "GO Biological Processes: Upregulated DEGs",
  output_name = "Figure_13_GO_upregulated_DEGs"
)

create_enrichment_plot(
  enrichment_table = combined_go_results,
  analysis_name = "Paired-validated DEGs",
  plot_title = "GO Biological Processes: Paired-Validated DEGs",
  output_name = "Figure_14_GO_validated_DEGs"
)

# ============================================================
# Create KEGG enrichment figure
# ============================================================

create_enrichment_plot(
  enrichment_table = combined_kegg_results,
  analysis_name = "All significant DEGs",
  plot_title = "KEGG Pathways: All Significant DEGs",
  output_name = "Figure_15_KEGG_all_DEGs"
)

# ============================================================
# Create Hallmark GSEA figure
# ============================================================

hallmark_plot_data <- hallmark_gsea[
  !is.na(hallmark_gsea$padj),
  ,
  drop = FALSE
]

hallmark_plot_data <- hallmark_plot_data[
  order(hallmark_plot_data$padj),
  ,
  drop = FALSE
]

hallmark_plot_data <- head(
  hallmark_plot_data,
  20
)

if (nrow(hallmark_plot_data) > 0) {
  
  hallmark_plot_data$Significance <- ifelse(
    hallmark_plot_data$padj < 0.05,
    "FDR < 0.05",
    "Not significant"
  )
  
  hallmark_plot <- ggplot(
    hallmark_plot_data,
    aes(
      x = NES,
      y = reorder(
        pathway,
        NES
      ),
      fill = Significance
    )
  ) +
    geom_col() +
    scale_fill_manual(
      values = c(
        "FDR < 0.05" = "#B2182B",
        "Not significant" = "grey70"
      )
    ) +
    labs(
      title = "Hallmark Gene-Set Enrichment Analysis",
      subtitle = "ChILI_TOL versus Control_AC",
      x = "Normalized enrichment score",
      y = NULL,
      fill = "Significance"
    ) +
    theme_classic(
      base_size = 11
    )
  
  ggsave(
    filename = file.path(
      figures_directory,
      "Figure_16_Hallmark_GSEA.pdf"
    ),
    plot = hallmark_plot,
    width = 10,
    height = 8
  )
  
  ggsave(
    filename = file.path(
      figures_directory,
      "Figure_16_Hallmark_GSEA.png"
    ),
    plot = hallmark_plot,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# ============================================================
# Create summary table
# ============================================================

enrichment_summary <- data.frame(
  Gene_set = c(
    "All primary significant DEGs",
    "Upregulated DEGs",
    "Downregulated DEGs",
    "Nominally paired-validated DEGs",
    "FDR paired-validated DEGs"
  ),
  
  Number_of_genes = c(
    length(significant_entrez),
    length(upregulated_entrez),
    length(downregulated_entrez),
    length(validated_entrez),
    length(fdr_validated_entrez)
  ),
  
  stringsAsFactors = FALSE
)

# ============================================================
# Create Excel workbook
# ============================================================

enrichment_workbook <- createWorkbook()

addWorksheet(
  enrichment_workbook,
  "Enrichment summary"
)

writeData(
  enrichment_workbook,
  sheet = "Enrichment summary",
  x = enrichment_summary
)

addWorksheet(
  enrichment_workbook,
  "GO results"
)

writeData(
  enrichment_workbook,
  sheet = "GO results",
  x = combined_go_results
)

addWorksheet(
  enrichment_workbook,
  "KEGG results"
)

writeData(
  enrichment_workbook,
  sheet = "KEGG results",
  x = combined_kegg_results
)

addWorksheet(
  enrichment_workbook,
  "Hallmark GSEA"
)

writeData(
  enrichment_workbook,
  sheet = "Hallmark GSEA",
  x = hallmark_gsea
)

for (sheet_name in names(enrichment_workbook)) {
  
  setColWidths(
    enrichment_workbook,
    sheet = sheet_name,
    cols = 1:15,
    widths = "auto"
  )
}

saveWorkbook(
  enrichment_workbook,
  file = file.path(
    tables_directory,
    "GSE287540_functional_enrichment_results.xlsx"
  ),
  overwrite = TRUE
)

# ============================================================
# Save R objects
# ============================================================

saveRDS(
  combined_go_results,
  file = file.path(
    objects_directory,
    "GSE287540_GO_results.rds"
  )
)

saveRDS(
  combined_kegg_results,
  file = file.path(
    objects_directory,
    "GSE287540_KEGG_results.rds"
  )
)

saveRDS(
  hallmark_gsea,
  file = file.path(
    objects_directory,
    "GSE287540_Hallmark_GSEA_results.rds"
  )
)

# ============================================================
# Final summary
# ============================================================

significant_hallmark_number <- sum(
  hallmark_gsea$padj < 0.05,
  na.rm = TRUE
)

cat("\n============================================\n")
cat("Functional enrichment analysis completed\n")
cat("============================================\n")

cat(
  "\nGO terms identified:",
  nrow(combined_go_results),
  "\n"
)

cat(
  "KEGG pathways identified:",
  nrow(combined_kegg_results),
  "\n"
)

cat(
  "Significant Hallmark pathways:",
  significant_hallmark_number,
  "\n"
)

cat("\nMain files created:\n")
cat("results/tables/GSE287540_functional_enrichment_results.xlsx\n")
cat("results/figures/Figure_12_GO_all_DEGs.png\n")
cat("results/figures/Figure_13_GO_upregulated_DEGs.png\n")
cat("results/figures/Figure_14_GO_validated_DEGs.png\n")
cat("results/figures/Figure_15_KEGG_all_DEGs.png\n")
cat("results/figures/Figure_16_Hallmark_GSEA.png\n")
go_results <- readRDS(
  "results/R_objects/GSE287540_GO_results.rds"
)

kegg_results <- readRDS(
  "results/R_objects/GSE287540_KEGG_results.rds"
)

hallmark_results <- readRDS(
  "results/R_objects/GSE287540_Hallmark_GSEA_results.rds"
)

cat("GO terms:", nrow(go_results), "\n")
cat("KEGG pathways:", nrow(kegg_results), "\n")
cat(
  "Significant Hallmark pathways:",
  sum(hallmark_results$padj < 0.05, na.rm = TRUE),
  "\n"
)