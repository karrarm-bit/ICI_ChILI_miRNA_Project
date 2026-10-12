# ============================================================
# Script 11: Prepare GSE274664 Human microRNA Count Data
# Project: Integrated mRNA-miRNA Analysis of ICI Response
# Dataset: GSE274664
#
# Purpose:
# 1. Read the original feature-count file
# 2. Remove the non-numeric title row
# 3. Convert all counts to numeric values
# 4. select human miRNA-associated features
# 5. combine matching _mature and _hairpin rows
# 6. prepare responder/non-responder metadata
# 7. validate and save the final count matrix
#
# Important:
# In this dataset, rows ending in both _hairpin and _mature
# contain mature-arm miRNA identifiers such as miR-21-5p.
# Matching rows are therefore combined by miRNA base name.
# ============================================================


# ------------------------------------------------------------
# 1. Clean the R environment
# ------------------------------------------------------------

rm(list = ls())

options(stringsAsFactors = FALSE)


# ------------------------------------------------------------
# 2. Check the project directory
# ------------------------------------------------------------

required_project_name <- "ICI_ChILI_miRNA_Project"

if (!grepl(required_project_name, getwd(), fixed = TRUE)) {
  stop(
    paste0(
      "The current working directory is not the expected project folder.\n",
      "Current directory: ",
      getwd()
    )
  )
}


# ------------------------------------------------------------
# 3. Check required packages
# ------------------------------------------------------------

required_packages <- c(
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
  library(openxlsx)
})


# ------------------------------------------------------------
# 4. Create output directories
# ------------------------------------------------------------

output_directories <- c(
  "data/processed",
  "results/R_objects",
  "results/tables"
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
# 5. Define and validate the input file
# ------------------------------------------------------------

raw_count_file <- paste0(
  "data/raw/",
  "GSE274664_CD19_EV_miRNA_readcounts.txt.gz"
)

if (!file.exists(raw_count_file)) {
  stop(
    paste0(
      "The raw count file was not found:\n",
      raw_count_file,
      "\nRun script 10_download_GSE274664_miRNA.R first."
    )
  )
}


# ------------------------------------------------------------
# 6. Read the original count file
# ------------------------------------------------------------

cat("Reading the original GSE274664 count file...\n")

raw_data <- read.delim(
  raw_count_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat(
  "Raw file:",
  nrow(raw_data),
  "rows and",
  ncol(raw_data),
  "columns\n"
)

if (!"sample" %in% colnames(raw_data)) {
  stop(
    "The required feature identifier column named 'sample' was not found."
  )
}

sample_columns <- setdiff(
  colnames(raw_data),
  "sample"
)

if (length(sample_columns) != 13) {
  stop(
    paste0(
      "Thirteen biological sample columns were expected, but ",
      length(sample_columns),
      " were found."
    )
  )
}


# ------------------------------------------------------------
# 7. Extract the sample labels
# ------------------------------------------------------------

title_row_index <- which(
  tolower(trimws(raw_data$sample)) == "title"
)

if (length(title_row_index) != 1) {
  stop(
    "Exactly one non-numeric title row was expected."
  )
}

sample_labels <- as.character(
  unlist(
    raw_data[
      title_row_index,
      sample_columns,
      drop = FALSE
    ],
    use.names = FALSE
  )
)

sample_label_lookup <- setNames(
  sample_labels,
  sample_columns
)

cat("\nSample labels detected:\n")

print(
  data.frame(
    Sample_ID = sample_columns,
    Sample_label = sample_labels,
    stringsAsFactors = FALSE
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. Remove the title row before numeric conversion
# ------------------------------------------------------------

raw_data <- raw_data[
  -title_row_index,
  ,
  drop = FALSE
]

rownames(raw_data) <- NULL

feature_ids <- trimws(
  as.character(raw_data$sample)
)

if (any(is.na(feature_ids) | feature_ids == "")) {
  stop(
    "Missing feature identifiers were detected."
  )
}

if (anyDuplicated(feature_ids) > 0) {
  stop(
    "Duplicated original feature identifiers were detected."
  )
}

cat(
  "\nFeature rows after removing the title row:",
  nrow(raw_data),
  "\n"
)


# ------------------------------------------------------------
# 9. Convert sample counts to numeric values
# ------------------------------------------------------------

cat("\nConverting sample counts to numeric values...\n")

numeric_count_list <- lapply(
  raw_data[
    ,
    sample_columns,
    drop = FALSE
  ],
  function(column_values) {
    suppressWarnings(
      as.numeric(
        gsub(
          ",",
          "",
          trimws(as.character(column_values))
        )
      )
    )
  }
)

count_matrix_all <- do.call(
  cbind,
  numeric_count_list
)

count_matrix_all <- as.matrix(
  count_matrix_all
)

rownames(count_matrix_all) <- feature_ids
colnames(count_matrix_all) <- sample_columns

storage.mode(count_matrix_all) <- "numeric"

if (anyNA(count_matrix_all)) {
  stop(
    paste0(
      "Numeric conversion generated ",
      sum(is.na(count_matrix_all)),
      " missing count values."
    )
  )
}

if (any(!is.finite(count_matrix_all))) {
  stop(
    "Non-finite count values were detected."
  )
}

if (any(count_matrix_all < 0)) {
  stop(
    "Negative count values were detected."
  )
}

if (any(
  abs(
    count_matrix_all -
    round(count_matrix_all)
  ) > 1e-6
)) {
  warning(
    "Non-integer values were detected and rounded."
  )
  
  count_matrix_all <- round(
    count_matrix_all
  )
}

cat("Numeric conversion completed successfully.\n")


# ------------------------------------------------------------
# 10. Classify the original features
# ------------------------------------------------------------

feature_type <- ifelse(
  grepl(
    "_mature$",
    feature_ids
  ),
  "Mature microRNA",
  ifelse(
    grepl(
      "_hairpin$",
      feature_ids
    ),
    "Hairpin-associated microRNA",
    "Other feature"
  )
)

feature_summary <- as.data.frame(
  table(feature_type),
  stringsAsFactors = FALSE
)

colnames(feature_summary) <- c(
  "Feature_type",
  "Number_of_features"
)

cat("\nOriginal feature types:\n")

print(
  feature_summary,
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Select human miRNA-associated rows
# ------------------------------------------------------------

human_mirna_rows <- grepl(
  "^hsa-",
  feature_ids
) &
  grepl(
    "_(mature|hairpin)$",
    feature_ids
  )

number_of_human_feature_rows <- sum(
  human_mirna_rows
)

if (number_of_human_feature_rows == 0) {
  stop(
    "No human microRNA-associated rows were detected."
  )
}

human_feature_matrix <- count_matrix_all[
  human_mirna_rows,
  ,
  drop = FALSE
]

human_feature_ids <- rownames(
  human_feature_matrix
)

human_feature_source <- ifelse(
  grepl(
    "_mature$",
    human_feature_ids
  ),
  "mature",
  "hairpin"
)

human_mirna_base_names <- sub(
  "_(mature|hairpin)$",
  "",
  human_feature_ids
)

cat(
  "\nHuman miRNA-associated rows:",
  nrow(human_feature_matrix),
  "\n"
)

cat(
  "Human mature rows:",
  sum(human_feature_source == "mature"),
  "\n"
)

cat(
  "Human hairpin-associated rows:",
  sum(human_feature_source == "hairpin"),
  "\n"
)

cat(
  "Unique human miRNA names before zero filtering:",
  length(unique(human_mirna_base_names)),
  "\n"
)


# ------------------------------------------------------------
# 12. Combine mature and hairpin-associated rows
# ------------------------------------------------------------

cat(
  "\nCombining matching mature and hairpin-associated rows...\n"
)

combined_mirna_counts <- rowsum(
  human_feature_matrix,
  group = human_mirna_base_names,
  reorder = FALSE
)

combined_mirna_counts <- as.matrix(
  combined_mirna_counts
)

storage.mode(combined_mirna_counts) <- "numeric"

colnames(combined_mirna_counts) <-
  sample_columns

if (anyDuplicated(
  rownames(combined_mirna_counts)
) > 0) {
  stop(
    "Duplicated miRNA names remain after row aggregation."
  )
}

if (anyNA(combined_mirna_counts)) {
  stop(
    "Missing values were generated during row aggregation."
  )
}

if (any(combined_mirna_counts < 0)) {
  stop(
    "Negative values were generated during row aggregation."
  )
}


# ------------------------------------------------------------
# 13. Remove completely unexpressed miRNAs
# ------------------------------------------------------------

total_counts_per_mirna <- rowSums(
  combined_mirna_counts
)

keep_expressed_mirnas <-
  total_counts_per_mirna > 0

number_before_zero_filtering <- nrow(
  combined_mirna_counts
)

number_zero_mirnas_removed <- sum(
  !keep_expressed_mirnas
)

combined_mirna_counts <- combined_mirna_counts[
  keep_expressed_mirnas,
  ,
  drop = FALSE
]

cat(
  "Human microRNAs before zero filtering:",
  number_before_zero_filtering,
  "\n"
)

cat(
  "Completely unexpressed microRNAs removed:",
  number_zero_mirnas_removed,
  "\n"
)

cat(
  "Human microRNAs retained:",
  nrow(combined_mirna_counts),
  "\n"
)


# ------------------------------------------------------------
# 14. Define response-group metadata
# ------------------------------------------------------------

responder_samples <- c(
  "Hg16_1263_CD19",
  "Hg16_884_CD19",
  "Hg17_2590_CD19",
  "Hg18_2387_CD19",
  "Hg19_1_CD19",
  "Hg19_4576_CD19",
  "Hg20_452_CD19",
  "Hg20_6327_CD19"
)

non_responder_samples <- c(
  "Hg17_3442_CD19",
  "Hg18_3693_CD19",
  "HG18.1201_CD19_NTLS",
  "HG18.924_TLS_3_CD19",
  "Hg19_699_CD19"
)

expected_samples <- c(
  responder_samples,
  non_responder_samples
)

if (anyDuplicated(expected_samples) > 0) {
  stop(
    "Duplicated sample identifiers were detected in the metadata."
  )
}

missing_samples <- setdiff(
  expected_samples,
  colnames(combined_mirna_counts)
)

unexpected_samples <- setdiff(
  colnames(combined_mirna_counts),
  expected_samples
)

if (length(missing_samples) > 0) {
  stop(
    paste0(
      "Expected samples missing from the matrix:\n",
      paste(missing_samples, collapse = ", ")
    )
  )
}

if (length(unexpected_samples) > 0) {
  stop(
    paste0(
      "Unexpected samples detected in the matrix:\n",
      paste(unexpected_samples, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------
# 15. Arrange the count-matrix columns
# ------------------------------------------------------------

combined_mirna_counts <- combined_mirna_counts[
  ,
  expected_samples,
  drop = FALSE
]

mirna_metadata <- data.frame(
  Matrix_column = seq_along(
    expected_samples
  ),
  Sample_ID = expected_samples,
  Sample_label = unname(
    sample_label_lookup[
      expected_samples
    ]
  ),
  Response_group = c(
    rep(
      "Responder",
      length(responder_samples)
    ),
    rep(
      "Non_responder",
      length(non_responder_samples)
    )
  ),
  stringsAsFactors = FALSE
)

mirna_metadata$Response_group <- factor(
  mirna_metadata$Response_group,
  levels = c(
    "Non_responder",
    "Responder"
  )
)

rownames(mirna_metadata) <-
  mirna_metadata$Sample_ID

if (!identical(
  colnames(combined_mirna_counts),
  mirna_metadata$Sample_ID
)) {
  stop(
    "The count matrix and metadata are not aligned."
  )
}


# ------------------------------------------------------------
# 16. Create the miRNA source annotation
# ------------------------------------------------------------

unique_mirna_names <- rownames(
  combined_mirna_counts
)

source_annotation_list <- lapply(
  unique_mirna_names,
  function(mirna_name) {
    matching_rows <- human_mirna_base_names ==
      mirna_name
    
    sources <- sort(
      unique(
        human_feature_source[
          matching_rows
        ]
      )
    )
    
    data.frame(
      Mature_miRNA = mirna_name,
      Original_feature_rows = sum(
        matching_rows
      ),
      Source_type = ifelse(
        length(sources) > 1,
        "mature_and_hairpin",
        sources
      ),
      Mature_row_available = any(
        human_feature_source[
          matching_rows
        ] == "mature"
      ),
      Hairpin_row_available = any(
        human_feature_source[
          matching_rows
        ] == "hairpin"
      ),
      Total_count = sum(
        combined_mirna_counts[
          mirna_name,
          ,
          drop = TRUE
        ]
      ),
      Number_of_expressing_samples = sum(
        combined_mirna_counts[
          mirna_name,
          ,
          drop = TRUE
        ] > 0
      ),
      stringsAsFactors = FALSE
    )
  }
)

mirna_annotation <- do.call(
  rbind,
  source_annotation_list
)

mirna_annotation <- mirna_annotation[
  order(
    mirna_annotation$Total_count,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

rownames(mirna_annotation) <- NULL


# ------------------------------------------------------------
# 17. Calculate sample-level statistics
# ------------------------------------------------------------

sample_library_sizes <- colSums(
  combined_mirna_counts
)

sample_detected_mirnas <- colSums(
  combined_mirna_counts > 0
)

sample_statistics <- data.frame(
  Sample_ID = colnames(
    combined_mirna_counts
  ),
  Sample_label = unname(
    sample_label_lookup[
      colnames(combined_mirna_counts)
    ]
  ),
  Response_group = as.character(
    mirna_metadata[
      colnames(combined_mirna_counts),
      "Response_group"
    ]
  ),
  Human_miRNA_library_size = as.numeric(
    sample_library_sizes
  ),
  Human_miRNA_library_size_thousand =
    as.numeric(sample_library_sizes) / 1000,
  Detected_human_miRNAs = as.numeric(
    sample_detected_mirnas
  ),
  stringsAsFactors = FALSE
)

zero_library_samples <- sample_statistics$Sample_ID[
  sample_statistics$Human_miRNA_library_size <= 0
]

if (length(zero_library_samples) > 0) {
  stop(
    paste0(
      "Zero-library samples remain:\n",
      paste(
        zero_library_samples,
        collapse = ", "
      )
    )
  )
}

cat("\nSample groups:\n")
print(table(mirna_metadata$Response_group))

cat("\nFinal sample-level miRNA statistics:\n")

print(
  sample_statistics,
  row.names = FALSE,
  digits = 5
)


# ------------------------------------------------------------
# 18. Validate the final matrix
# ------------------------------------------------------------

if (ncol(combined_mirna_counts) != 13) {
  stop(
    paste0(
      "The final matrix contains ",
      ncol(combined_mirna_counts),
      " samples instead of 13."
    )
  )
}

if (sum(
  mirna_metadata$Response_group ==
  "Responder"
) != 8) {
  stop(
    "The final metadata must contain eight responders."
  )
}

if (sum(
  mirna_metadata$Response_group ==
  "Non_responder"
) != 5) {
  stop(
    "The final metadata must contain five non-responders."
  )
}


# ------------------------------------------------------------
# 19. Save the corrected R objects
# ------------------------------------------------------------

saveRDS(
  combined_mirna_counts,
  file = paste0(
    "results/R_objects/",
    "GSE274664_mature_miRNA_count_matrix.rds"
  )
)

saveRDS(
  mirna_metadata,
  file = paste0(
    "results/R_objects/",
    "GSE274664_miRNA_metadata.rds"
  )
)

saveRDS(
  mirna_annotation,
  file = paste0(
    "results/R_objects/",
    "GSE274664_miRNA_annotation.rds"
  )
)


# ------------------------------------------------------------
# 20. Save processed CSV files
# ------------------------------------------------------------

write.csv(
  mirna_metadata,
  file = paste0(
    "data/processed/",
    "GSE274664_miRNA_metadata.csv"
  ),
  row.names = FALSE
)

write.csv(
  sample_statistics,
  file = paste0(
    "data/processed/",
    "GSE274664_miRNA_sample_statistics.csv"
  ),
  row.names = FALSE
)

write.csv(
  mirna_annotation,
  file = paste0(
    "data/processed/",
    "GSE274664_miRNA_annotation.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 21. Create the preparation summary
# ------------------------------------------------------------

preparation_summary <- data.frame(
  Measurement = c(
    "Original feature rows",
    "Biological samples",
    "Human mature rows",
    "Human hairpin-associated rows",
    "Unique human miRNAs before zero filtering",
    "Completely unexpressed miRNAs removed",
    "Final human miRNAs",
    "Responder samples",
    "Non-responder samples",
    "Zero-library samples"
  ),
  Value = c(
    nrow(raw_data),
    length(sample_columns),
    sum(human_feature_source == "mature"),
    sum(human_feature_source == "hairpin"),
    number_before_zero_filtering,
    number_zero_mirnas_removed,
    nrow(combined_mirna_counts),
    sum(
      mirna_metadata$Response_group ==
        "Responder"
    ),
    sum(
      mirna_metadata$Response_group ==
        "Non_responder"
    ),
    length(zero_library_samples)
  ),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 22. Create the Excel workbook
# ------------------------------------------------------------

preparation_workbook <- createWorkbook()

addWorksheet(
  preparation_workbook,
  "Preparation_summary"
)

writeData(
  preparation_workbook,
  sheet = "Preparation_summary",
  x = preparation_summary
)

addWorksheet(
  preparation_workbook,
  "Sample_metadata"
)

writeData(
  preparation_workbook,
  sheet = "Sample_metadata",
  x = mirna_metadata,
  rowNames = FALSE
)

addWorksheet(
  preparation_workbook,
  "Sample_statistics"
)

writeData(
  preparation_workbook,
  sheet = "Sample_statistics",
  x = sample_statistics
)

addWorksheet(
  preparation_workbook,
  "Feature_types"
)

writeData(
  preparation_workbook,
  sheet = "Feature_types",
  x = feature_summary
)

addWorksheet(
  preparation_workbook,
  "miRNA_annotation"
)

writeData(
  preparation_workbook,
  sheet = "miRNA_annotation",
  x = mirna_annotation
)

header_style <- createStyle(
  fontColour = "#FFFFFF",
  fgFill = "#1F4E78",
  textDecoration = "bold",
  halign = "center"
)

for (sheet_name in names(preparation_workbook)) {
  sheet_data <- readWorkbook(
    preparation_workbook,
    sheet = sheet_name
  )
  
  if (ncol(sheet_data) > 0) {
    addStyle(
      preparation_workbook,
      sheet = sheet_name,
      style = header_style,
      rows = 1,
      cols = seq_len(ncol(sheet_data)),
      gridExpand = TRUE
    )
    
    setColWidths(
      preparation_workbook,
      sheet = sheet_name,
      cols = seq_len(ncol(sheet_data)),
      widths = "auto"
    )
  }
  
  freezePane(
    preparation_workbook,
    sheet = sheet_name,
    firstRow = TRUE
  )
}

saveWorkbook(
  preparation_workbook,
  file = paste0(
    "results/tables/",
    "GSE274664_miRNA_preparation_summary.xlsx"
  ),
  overwrite = TRUE
)


# ------------------------------------------------------------
# 23. Print final confirmation
# ------------------------------------------------------------

cat("\n")
cat("============================================\n")
cat("GSE274664 microRNA preparation completed\n")
cat("============================================\n")

cat(
  "\nFinal microRNA matrix:",
  nrow(combined_mirna_counts),
  "human microRNAs and",
  ncol(combined_mirna_counts),
  "samples\n"
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
  "Zero-library samples:",
  length(zero_library_samples),
  "\n"
)

cat(
  "Minimum human-miRNA library size:",
  min(sample_statistics$Human_miRNA_library_size),
  "\n"
)

cat(
  "Maximum human-miRNA library size:",
  max(sample_statistics$Human_miRNA_library_size),
  "\n"
)

cat("\nFiles created:\n")
cat("results/R_objects/GSE274664_mature_miRNA_count_matrix.rds\n")
cat("results/R_objects/GSE274664_miRNA_metadata.rds\n")
cat("results/R_objects/GSE274664_miRNA_annotation.rds\n")
cat("data/processed/GSE274664_miRNA_metadata.csv\n")
cat("data/processed/GSE274664_miRNA_sample_statistics.csv\n")
cat("data/processed/GSE274664_miRNA_annotation.csv\n")
cat("results/tables/GSE274664_miRNA_preparation_summary.xlsx\n")