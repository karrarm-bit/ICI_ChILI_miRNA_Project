# ============================================================
# Script 16: GSE274664 miRNA Leave-One-Out Sensitivity Analysis
# Project: Integrated mRNA-miRNA Analysis of ICI Response
#
# Purpose:
# 1. Test whether the four FDR-significant miRNAs are robust
# 2. omit one sample at a time
# 3. refit the edgeR quasi-likelihood model
# 4. evaluate fold-change direction consistency
# 5. evaluate dependence on individual samples
# 6. examine the effect of omitting NR4 specifically
#
# Important:
# No sample is permanently removed by this script.
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

de_results_file <- paste0(
  "results/R_objects/",
  "GSE274664_miRNA_DE_results.rds"
)

log_cpm_file <- paste0(
  "results/R_objects/",
  "GSE274664_final_miRNA_logCPM.rds"
)

required_files <- c(
  dge_file,
  metadata_file,
  de_results_file,
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
      "\nRun scripts 14 and 15 first."
    )
  )
}


# ------------------------------------------------------------
# 5. Load input objects
# ------------------------------------------------------------

cat("Loading GSE274664 differential-expression data...\n")

final_y <- readRDS(
  dge_file
)

final_metadata <- readRDS(
  metadata_file
)

baseline_results <- readRDS(
  de_results_file
)

final_log_cpm <- readRDS(
  log_cpm_file
)

final_metadata <- as.data.frame(
  final_metadata
)

baseline_results <- as.data.frame(
  baseline_results
)

final_log_cpm <- as.matrix(
  final_log_cpm
)

cat(
  "Loaded:",
  nrow(final_y),
  "miRNAs and",
  ncol(final_y),
  "samples\n"
)


# ------------------------------------------------------------
# 6. Validate and align the data
# ------------------------------------------------------------

final_metadata <- final_metadata[
  match(
    colnames(final_y$counts),
    final_metadata$Sample_ID
  ),
  ,
  drop = FALSE
]

rownames(final_metadata) <-
  final_metadata$Sample_ID

final_metadata$Response_group <- factor(
  as.character(
    final_metadata$Response_group
  ),
  levels = c(
    "Non_responder",
    "Responder"
  )
)

if (anyNA(final_metadata$Sample_ID)) {
  stop(
    "Some samples could not be matched to the metadata."
  )
}

if (!identical(
  colnames(final_y$counts),
  final_metadata$Sample_ID
)) {
  stop(
    "The count matrix and metadata are not aligned."
  )
}

if (!identical(
  colnames(final_y$counts),
  colnames(final_log_cpm)
)) {
  stop(
    "The count matrix and logCPM matrix are not aligned."
  )
}


# ------------------------------------------------------------
# 7. Identify the baseline FDR-significant miRNAs
# ------------------------------------------------------------

primary_candidates <- baseline_results[
  baseline_results$Primary_significant,
  ,
  drop = FALSE
]

if (nrow(primary_candidates) == 0) {
  stop(
    "No FDR-significant baseline miRNAs were available."
  )
}

candidate_mirnas <- primary_candidates$miRNA

missing_candidate_mirnas <- setdiff(
  candidate_mirnas,
  rownames(final_y$counts)
)

if (length(missing_candidate_mirnas) > 0) {
  stop(
    paste0(
      "Candidate miRNAs missing from the count matrix: ",
      paste(
        missing_candidate_mirnas,
        collapse = ", "
      )
    )
  )
}

cat("\nBaseline FDR-significant miRNAs:\n")

print(
  primary_candidates[
    ,
    c(
      "miRNA",
      "logFC",
      "PValue",
      "FDR",
      "Direction"
    )
  ],
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. Define a TMM normalization function
# ------------------------------------------------------------

apply_tmm_normalization <- function(dge_object) {
  if (
    "normLibSizes" %in%
    getNamespaceExports("edgeR")
  ) {
    dge_object <- normLibSizes(
      dge_object,
      method = "TMM"
    )
  } else {
    dge_object <- calcNormFactors(
      dge_object,
      method = "TMM"
    )
  }
  
  return(dge_object)
}


# ------------------------------------------------------------
# 9. Run leave-one-out models
# ------------------------------------------------------------

cat("\nRunning leave-one-out sensitivity models...\n")

sample_ids <- colnames(
  final_y$counts
)

loo_all_results <- vector(
  mode = "list",
  length = length(sample_ids)
)

loo_candidate_results <- vector(
  mode = "list",
  length = length(sample_ids)
)

for (sample_index in seq_along(sample_ids)) {
  omitted_sample <- sample_ids[
    sample_index
  ]
  
  retained_samples <- setdiff(
    sample_ids,
    omitted_sample
  )
  
  loo_metadata <- final_metadata[
    retained_samples,
    ,
    drop = FALSE
  ]
  
  loo_metadata$Response_group <- droplevels(
    loo_metadata$Response_group
  )
  
  if (nlevels(loo_metadata$Response_group) != 2) {
    stop(
      paste0(
        "Both response groups were not retained after omitting ",
        omitted_sample,
        "."
      )
    )
  }
  
  loo_counts <- final_y$counts[
    ,
    retained_samples,
    drop = FALSE
  ]
  
  loo_y <- DGEList(
    counts = loo_counts,
    samples = data.frame(
      Sample_ID = loo_metadata$Sample_ID,
      Response_group =
        loo_metadata$Response_group,
      row.names = loo_metadata$Sample_ID
    )
  )
  
  loo_y <- apply_tmm_normalization(
    loo_y
  )
  
  loo_design <- model.matrix(
    ~ Response_group,
    data = loo_metadata
  )
  
  rownames(loo_design) <-
    loo_metadata$Sample_ID
  
  if (
    qr(loo_design)$rank !=
    ncol(loo_design)
  ) {
    stop(
      paste0(
        "The design became rank-deficient after omitting ",
        omitted_sample,
        "."
      )
    )
  }
  
  loo_y <- estimateDisp(
    loo_y,
    design = loo_design,
    robust = TRUE
  )
  
  loo_fit <- glmQLFit(
    loo_y,
    design = loo_design,
    robust = TRUE
  )
  
  loo_test <- glmQLFTest(
    loo_fit,
    coef = "Response_groupResponder"
  )
  
  loo_table <- topTags(
    loo_test,
    n = Inf,
    sort.by = "PValue"
  )$table
  
  loo_table <- data.frame(
    miRNA = rownames(loo_table),
    loo_table,
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  omitted_label <- as.character(
    final_metadata[
      omitted_sample,
      "Sample_label"
    ]
  )
  
  omitted_group <- as.character(
    final_metadata[
      omitted_sample,
      "Response_group"
    ]
  )
  
  loo_table$Omitted_sample <-
    omitted_sample
  
  loo_table$Omitted_label <-
    omitted_label
  
  loo_table$Omitted_group <-
    omitted_group
  
  loo_all_results[[
    sample_index
  ]] <- loo_table
  
  candidate_table <- loo_table[
    match(
      candidate_mirnas,
      loo_table$miRNA
    ),
    ,
    drop = FALSE
  ]
  
  loo_candidate_results[[
    sample_index
  ]] <- candidate_table
  
  cat(
    "Completed model",
    sample_index,
    "of",
    length(sample_ids),
    "- omitted:",
    omitted_label,
    "\n"
  )
}

loo_all_results <- do.call(
  rbind,
  loo_all_results
)

loo_candidate_results <- do.call(
  rbind,
  loo_candidate_results
)

rownames(loo_all_results) <- NULL
rownames(loo_candidate_results) <- NULL


# ------------------------------------------------------------
# 10. Add baseline directions and sensitivity indicators
# ------------------------------------------------------------

baseline_lookup <- primary_candidates[
  match(
    loo_candidate_results$miRNA,
    primary_candidates$miRNA
  ),
  c(
    "miRNA",
    "logFC",
    "PValue",
    "FDR"
  ),
  drop = FALSE
]

loo_candidate_results$Baseline_logFC <-
  baseline_lookup$logFC

loo_candidate_results$Baseline_PValue <-
  baseline_lookup$PValue

loo_candidate_results$Baseline_FDR <-
  baseline_lookup$FDR

loo_candidate_results$Direction_concordant <- (
  sign(loo_candidate_results$logFC) ==
    sign(loo_candidate_results$Baseline_logFC)
)

loo_candidate_results$Nominal_with_effect <- (
  loo_candidate_results$PValue < 0.05 &
    abs(loo_candidate_results$logFC) >= 1
)

loo_candidate_results$FDR_with_effect <- (
  loo_candidate_results$FDR < 0.05 &
    abs(loo_candidate_results$logFC) >= 1
)


# ------------------------------------------------------------
# 11. Summarize robustness for each candidate
# ------------------------------------------------------------

candidate_summary_list <- lapply(
  candidate_mirnas,
  function(mirna_name) {
    candidate_data <- loo_candidate_results[
      loo_candidate_results$miRNA ==
        mirna_name,
      ,
      drop = FALSE
    ]
    
    baseline_row <- primary_candidates[
      primary_candidates$miRNA ==
        mirna_name,
      ,
      drop = FALSE
    ]
    
    nr4_removed_row <- candidate_data[
      candidate_data$Omitted_label ==
        "NR4",
      ,
      drop = FALSE
    ]
    
    data.frame(
      miRNA = mirna_name,
      Baseline_logFC =
        baseline_row$logFC[1],
      Baseline_PValue =
        baseline_row$PValue[1],
      Baseline_FDR =
        baseline_row$FDR[1],
      LOO_models = nrow(
        candidate_data
      ),
      Direction_concordant_models = sum(
        candidate_data$Direction_concordant,
        na.rm = TRUE
      ),
      Direction_concordance_percent = round(
        100 * mean(
          candidate_data$Direction_concordant,
          na.rm = TRUE
        ),
        1
      ),
      Median_LOO_logFC = median(
        candidate_data$logFC,
        na.rm = TRUE
      ),
      Minimum_LOO_logFC = min(
        candidate_data$logFC,
        na.rm = TRUE
      ),
      Maximum_LOO_logFC = max(
        candidate_data$logFC,
        na.rm = TRUE
      ),
      Nominal_supported_models = sum(
        candidate_data$Nominal_with_effect,
        na.rm = TRUE
      ),
      FDR_supported_models = sum(
        candidate_data$FDR_with_effect,
        na.rm = TRUE
      ),
      NR4_removed_logFC = ifelse(
        nrow(nr4_removed_row) == 1,
        nr4_removed_row$logFC,
        NA_real_
      ),
      NR4_removed_PValue = ifelse(
        nrow(nr4_removed_row) == 1,
        nr4_removed_row$PValue,
        NA_real_
      ),
      NR4_removed_FDR = ifelse(
        nrow(nr4_removed_row) == 1,
        nr4_removed_row$FDR,
        NA_real_
      ),
      stringsAsFactors = FALSE
    )
  }
)

candidate_sensitivity_summary <- do.call(
  rbind,
  candidate_summary_list
)

rownames(candidate_sensitivity_summary) <-
  NULL


# ------------------------------------------------------------
# 12. Assign transparent stability classifications
# ------------------------------------------------------------

candidate_sensitivity_summary$Stability_class <- ifelse(
  candidate_sensitivity_summary$
    Direction_concordance_percent >= 90 &
    candidate_sensitivity_summary$
    Nominal_supported_models >= 8,
  "Strong",
  ifelse(
    candidate_sensitivity_summary$
      Direction_concordance_percent >= 80 &
      candidate_sensitivity_summary$
      Nominal_supported_models >= 4,
    "Moderate",
    ifelse(
      candidate_sensitivity_summary$
        Direction_concordance_percent >= 80,
      "Direction_only",
      "Unstable"
    )
  )
)

candidate_sensitivity_summary <-
  candidate_sensitivity_summary[
    order(
      factor(
        candidate_sensitivity_summary$
          Stability_class,
        levels = c(
          "Strong",
          "Moderate",
          "Direction_only",
          "Unstable"
        )
      ),
      candidate_sensitivity_summary$
        Baseline_FDR
    ),
    ,
    drop = FALSE
  ]


# ------------------------------------------------------------
# 13. Calculate expression and detection summaries
# ------------------------------------------------------------

candidate_expression_list <- lapply(
  candidate_mirnas,
  function(mirna_name) {
    raw_candidate_counts <- final_y$counts[
      mirna_name,
      ,
      drop = TRUE
    ]
    
    normalized_candidate_expression <-
      final_log_cpm[
        mirna_name,
        ,
        drop = TRUE
      ]
    
    data.frame(
      miRNA = mirna_name,
      Mean_logCPM_Non_responder = mean(
        normalized_candidate_expression[
          final_metadata$Response_group ==
            "Non_responder"
        ]
      ),
      Mean_logCPM_Responder = mean(
        normalized_candidate_expression[
          final_metadata$Response_group ==
            "Responder"
        ]
      ),
      Median_logCPM_Non_responder = median(
        normalized_candidate_expression[
          final_metadata$Response_group ==
            "Non_responder"
        ]
      ),
      Median_logCPM_Responder = median(
        normalized_candidate_expression[
          final_metadata$Response_group ==
            "Responder"
        ]
      ),
      Expressing_Non_responders = sum(
        raw_candidate_counts[
          final_metadata$Response_group ==
            "Non_responder"
        ] > 0
      ),
      Total_Non_responders = sum(
        final_metadata$Response_group ==
          "Non_responder"
      ),
      Expressing_Responders = sum(
        raw_candidate_counts[
          final_metadata$Response_group ==
            "Responder"
        ] > 0
      ),
      Total_Responders = sum(
        final_metadata$Response_group ==
          "Responder"
      ),
      stringsAsFactors = FALSE
    )
  }
)

candidate_expression_summary <- do.call(
  rbind,
  candidate_expression_list
)

rownames(candidate_expression_summary) <-
  NULL


# ------------------------------------------------------------
# 14. Create the leave-one-out fold-change figure
# ------------------------------------------------------------

plot_data <- loo_candidate_results

plot_data$Omitted_label <- factor(
  plot_data$Omitted_label,
  levels = final_metadata$Sample_label
)

baseline_plot_data <- primary_candidates[
  primary_candidates$miRNA %in%
    candidate_mirnas,
  c(
    "miRNA",
    "logFC"
  ),
  drop = FALSE
]

loo_plot <- ggplot(
  plot_data,
  aes(
    x = Omitted_label,
    y = logFC,
    color = Omitted_group
  )
) +
  geom_hline(
    data = baseline_plot_data,
    mapping = aes(
      yintercept = logFC
    ),
    inherit.aes = FALSE,
    color = "black",
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_hline(
    yintercept = 0,
    color = "grey50",
    linewidth = 0.4
  ) +
  geom_point(
    size = 3,
    alpha = 0.85
  ) +
  facet_wrap(
    ~ miRNA,
    scales = "free_y",
    ncol = 2
  ) +
  scale_color_manual(
    values = c(
      "Non_responder" = "#D55E00",
      "Responder" = "#0072B2"
    )
  ) +
  labs(
    title = "Leave-One-Out Stability of Significant miRNAs",
    subtitle = paste0(
      "Each point represents the model after omitting one sample; ",
      "dashed line represents the complete-data logFC"
    ),
    x = "Omitted sample",
    y = "Responder versus Non-responder log2 fold change",
    color = "Omitted group"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "bottom"
  )

ggsave(
  filename = paste0(
    "results/figures/",
    "Figure_26_miRNA_leave_one_out_logFC.png"
  ),
  plot = loo_plot,
  width = 13,
  height = 9,
  dpi = 300
)

ggsave(
  filename = paste0(
    "results/figures/",
    "Figure_26_miRNA_leave_one_out_logFC.pdf"
  ),
  plot = loo_plot,
  width = 13,
  height = 9
)


# ------------------------------------------------------------
# 15. Create candidate expression data
# ------------------------------------------------------------

candidate_expression_long <- do.call(
  rbind,
  lapply(
    candidate_mirnas,
    function(mirna_name) {
      data.frame(
        Sample_ID = colnames(
          final_log_cpm
        ),
        Sample_label = as.character(
          final_metadata[
            colnames(final_log_cpm),
            "Sample_label"
          ]
        ),
        Response_group = as.character(
          final_metadata[
            colnames(final_log_cpm),
            "Response_group"
          ]
        ),
        miRNA = mirna_name,
        logCPM = as.numeric(
          final_log_cpm[
            mirna_name,
            colnames(final_log_cpm)
          ]
        ),
        Raw_count = as.numeric(
          final_y$counts[
            mirna_name,
            colnames(final_log_cpm)
          ]
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

candidate_expression_long$Response_group <-
  factor(
    candidate_expression_long$Response_group,
    levels = c(
      "Non_responder",
      "Responder"
    )
  )


# ------------------------------------------------------------
# 16. Create labeled candidate-expression figure
# ------------------------------------------------------------

candidate_expression_plot <- ggplot(
  candidate_expression_long,
  aes(
    x = Response_group,
    y = logCPM,
    fill = Response_group
  )
) +
  geom_boxplot(
    width = 0.6,
    alpha = 0.6,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.10,
    size = 2.5,
    alpha = 0.85
  ) +
  facet_wrap(
    ~ miRNA,
    scales = "free_y",
    ncol = 2
  ) +
  scale_fill_manual(
    values = c(
      "Non_responder" = "#D55E00",
      "Responder" = "#0072B2"
    )
  ) +
  labs(
    title = "Expression of FDR-Significant miRNAs",
    subtitle = "Final TMM-normalized log2 CPM values",
    x = "Response group",
    y = "log2 CPM"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    ),
    legend.position = "none"
  )

ggsave(
  filename = paste0(
    "results/figures/",
    "Figure_27_significant_miRNA_expression.png"
  ),
  plot = candidate_expression_plot,
  width = 11,
  height = 8,
  dpi = 300
)

ggsave(
  filename = paste0(
    "results/figures/",
    "Figure_27_significant_miRNA_expression.pdf"
  ),
  plot = candidate_expression_plot,
  width = 11,
  height = 8
)


# ------------------------------------------------------------
# 17. Save CSV files
# ------------------------------------------------------------

write.csv(
  loo_candidate_results,
  file = paste0(
    "data/processed/",
    "GSE274664_miRNA_leave_one_out_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  candidate_sensitivity_summary,
  file = paste0(
    "data/processed/",
    "GSE274664_miRNA_sensitivity_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  candidate_expression_summary,
  file = paste0(
    "data/processed/",
    "GSE274664_significant_miRNA_expression_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 18. Create the Excel workbook
# ------------------------------------------------------------

sensitivity_workbook <- createWorkbook()

addWorksheet(
  sensitivity_workbook,
  "Sensitivity_summary"
)

writeData(
  sensitivity_workbook,
  sheet = "Sensitivity_summary",
  x = candidate_sensitivity_summary
)

addWorksheet(
  sensitivity_workbook,
  "LOO_candidate_results"
)

writeData(
  sensitivity_workbook,
  sheet = "LOO_candidate_results",
  x = loo_candidate_results
)

addWorksheet(
  sensitivity_workbook,
  "Expression_summary"
)

writeData(
  sensitivity_workbook,
  sheet = "Expression_summary",
  x = candidate_expression_summary
)

addWorksheet(
  sensitivity_workbook,
  "Candidate_expression"
)

writeData(
  sensitivity_workbook,
  sheet = "Candidate_expression",
  x = candidate_expression_long
)

saveWorkbook(
  sensitivity_workbook,
  file = paste0(
    "results/tables/",
    "GSE274664_miRNA_sensitivity_analysis.xlsx"
  ),
  overwrite = TRUE
)


# ------------------------------------------------------------
# 19. Save R objects
# ------------------------------------------------------------

saveRDS(
  loo_candidate_results,
  file = paste0(
    "results/R_objects/",
    "GSE274664_miRNA_leave_one_out_results.rds"
  )
)

saveRDS(
  candidate_sensitivity_summary,
  file = paste0(
    "results/R_objects/",
    "GSE274664_miRNA_sensitivity_summary.rds"
  )
)

saveRDS(
  candidate_expression_summary,
  file = paste0(
    "results/R_objects/",
    "GSE274664_significant_miRNA_expression_summary.rds"
  )
)


# ------------------------------------------------------------
# 20. Print final confirmation
# ------------------------------------------------------------

cat("\n")
cat("============================================\n")
cat("GSE274664 miRNA sensitivity analysis completed\n")
cat("============================================\n")

cat(
  "\nLeave-one-out models:",
  length(sample_ids),
  "\n"
)

cat(
  "Baseline candidates:",
  length(candidate_mirnas),
  "\n"
)

cat("\nCandidate sensitivity summary:\n")

print(
  candidate_sensitivity_summary[
    ,
    c(
      "miRNA",
      "Baseline_logFC",
      "Direction_concordance_percent",
      "Nominal_supported_models",
      "FDR_supported_models",
      "NR4_removed_logFC",
      "NR4_removed_PValue",
      "NR4_removed_FDR",
      "Stability_class"
    )
  ],
  row.names = FALSE,
  digits = 4
)

cat("\nCandidate detection summary:\n")

print(
  candidate_expression_summary[
    ,
    c(
      "miRNA",
      "Expressing_Non_responders",
      "Total_Non_responders",
      "Expressing_Responders",
      "Total_Responders"
    )
  ],
  row.names = FALSE
)

cat("\nInterpretation guide:\n")
cat("Strong: consistent direction and broad nominal support.\n")
cat("Moderate: generally consistent with partial statistical support.\n")
cat("Direction_only: direction is stable but P-values are sensitive.\n")
cat("Unstable: direction changes across leave-one-out models.\n")

cat("\nMain files created:\n")
cat("results/tables/GSE274664_miRNA_sensitivity_analysis.xlsx\n")
cat("results/R_objects/GSE274664_miRNA_sensitivity_summary.rds\n")
cat("results/figures/Figure_26_miRNA_leave_one_out_logFC.png\n")
cat("results/figures/Figure_27_significant_miRNA_expression.png\n")