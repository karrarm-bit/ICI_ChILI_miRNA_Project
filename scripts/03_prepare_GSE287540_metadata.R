# ============================================================
# Script 03: Prepare GSE287540 Sample Metadata
# Project: ICI-Induced Liver Injury and miRNA Analysis
# ============================================================

# Clear the R environment
rm(list = ls())

# Load required packages
suppressPackageStartupMessages({
  library(GEOquery)
  library(openxlsx)
})

# Define the project directory
project_directory <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

# Define input and output directories
raw_data_directory <- file.path(
  project_directory,
  "data",
  "raw"
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

objects_directory <- file.path(
  project_directory,
  "results",
  "R_objects"
)

# Define the metadata input file
metadata_file <- file.path(
  raw_data_directory,
  "GSE287540_series_matrix.txt.gz"
)

# Check whether the metadata file exists
if (!file.exists(metadata_file)) {
  stop(
    "The GSE287540 metadata file was not found. ",
    "Run 02_download_GSE287540.R first."
  )
}

cat("\nReading GSE287540 metadata...\n")

# Read the GEO series matrix from the local file
gse_object <- GEOquery::getGEO(
  filename = metadata_file
)

# Extract the sample metadata
raw_metadata <- Biobase::pData(gse_object)

cat(
  "Number of samples detected:",
  nrow(raw_metadata),
  "\n"
)

# ============================================================
# Function to extract GEO characteristics
# ============================================================

extract_characteristic <- function(metadata, variable_name) {
  
  characteristic_columns <- grep(
    pattern = "^characteristics_ch1",
    x = colnames(metadata),
    value = TRUE
  )
  
  if (length(characteristic_columns) == 0) {
    stop("No GEO characteristic columns were detected.")
  }
  
  extracted_values <- apply(
    metadata[, characteristic_columns, drop = FALSE],
    MARGIN = 1,
    FUN = function(sample_values) {
      
      sample_values <- as.character(sample_values)
      
      matching_value <- sample_values[
        grepl(
          pattern = paste0("^", variable_name, ":"),
          x = sample_values,
          ignore.case = TRUE
        )
      ]
      
      if (length(matching_value) == 0) {
        return(NA_character_)
      }
      
      cleaned_value <- sub(
        pattern = "^[^:]+:[[:space:]]*",
        replacement = "",
        x = matching_value[1]
      )
      
      return(trimws(cleaned_value))
    }
  )
  
  return(unname(extracted_values))
}

# ============================================================
# Create the cleaned metadata table
# ============================================================

clean_metadata <- data.frame(
  Sample_ID = rownames(raw_metadata),
  
  Sample_title = as.character(
    raw_metadata$title
  ),
  
  Subject_ID = extract_characteristic(
    raw_metadata,
    "subject"
  ),
  
  Group = tolower(
    extract_characteristic(
      raw_metadata,
      "group"
    )
  ),
  
  Time_point = toupper(
    extract_characteristic(
      raw_metadata,
      "time"
    )
  ),
  
  Sex = tolower(
    extract_characteristic(
      raw_metadata,
      "sex"
    )
  ),
  
  Age = suppressWarnings(
    as.numeric(
      extract_characteristic(
        raw_metadata,
        "age"
      )
    )
  ),
  
  Treatment = extract_characteristic(
    raw_metadata,
    "treatment"
  ),
  
  Steroid_use = tolower(
    extract_characteristic(
      raw_metadata,
      "steroid"
    )
  ),
  
  Batch = extract_characteristic(
    raw_metadata,
    "batch"
  ),
  
  stringsAsFactors = FALSE
)

# ============================================================
# Add simplified descriptive variables
# ============================================================

clean_metadata$Group_label <- ifelse(
  clean_metadata$Group == "chili",
  "ChILI",
  ifelse(
    clean_metadata$Group == "control",
    "Control",
    NA
  )
)

clean_metadata$Time_label <- ifelse(
  clean_metadata$Time_point == "BC",
  "Before ICI treatment",
  ifelse(
    clean_metadata$Time_point == "AC",
    "After ICI treatment - control",
    ifelse(
      clean_metadata$Time_point == "TOL",
      "Time of liver injury",
      ifelse(
        clean_metadata$Time_point == "FU1",
        "First follow-up",
        ifelse(
          clean_metadata$Time_point == "FU2",
          "Second follow-up",
          NA
        )
      )
    )
  )
)

clean_metadata$Treatment_class <- ifelse(
  grepl(
    pattern = "Ipilimumab.*Nivolumab|Nivolumab.*Ipilimumab",
    x = clean_metadata$Treatment,
    ignore.case = TRUE
  ),
  "Dual Ipilimumab and Nivolumab",
  "Other ICI treatment"
)

# ============================================================
# Check the cleaned metadata
# ============================================================

required_columns <- c(
  "Sample_ID",
  "Subject_ID",
  "Group",
  "Time_point"
)

missing_required_values <- colSums(
  is.na(
    clean_metadata[, required_columns]
  )
)

cat("\nMissing values in required variables:\n")
print(missing_required_values)

# Check for duplicated sample identifiers
duplicated_samples <- sum(
  duplicated(clean_metadata$Sample_ID)
)

cat(
  "\nNumber of duplicated sample IDs:",
  duplicated_samples,
  "\n"
)

# Count unique subjects
number_of_subjects <- length(
  unique(clean_metadata$Subject_ID)
)

cat(
  "Number of unique subjects:",
  number_of_subjects,
  "\n"
)

# ============================================================
# Create summary tables
# ============================================================

group_summary <- as.data.frame(
  table(clean_metadata$Group_label)
)

colnames(group_summary) <- c(
  "Group",
  "Number_of_samples"
)

time_summary <- as.data.frame(
  table(clean_metadata$Time_point)
)

colnames(time_summary) <- c(
  "Time_point",
  "Number_of_samples"
)

group_time_summary <- as.data.frame(
  table(
    clean_metadata$Group_label,
    clean_metadata$Time_point
  )
)

colnames(group_time_summary) <- c(
  "Group",
  "Time_point",
  "Number_of_samples"
)

treatment_summary <- as.data.frame(
  table(clean_metadata$Treatment)
)

colnames(treatment_summary) <- c(
  "Treatment",
  "Number_of_samples"
)

sex_summary <- as.data.frame(
  table(clean_metadata$Sex)
)

colnames(sex_summary) <- c(
  "Sex",
  "Number_of_samples"
)

# Display the main summaries
cat("\nSamples by study group:\n")
print(group_summary)

cat("\nSamples by time point:\n")
print(time_summary)

cat("\nSamples by group and time point:\n")
print(group_time_summary)

cat("\nSamples by treatment:\n")
print(treatment_summary)

# ============================================================
# Save the cleaned metadata
# ============================================================

write.csv(
  clean_metadata,
  file = file.path(
    processed_data_directory,
    "GSE287540_metadata_clean.csv"
  ),
  row.names = FALSE
)

saveRDS(
  clean_metadata,
  file = file.path(
    objects_directory,
    "GSE287540_metadata_clean.rds"
  )
)

# Create an Excel workbook
metadata_workbook <- createWorkbook()

addWorksheet(
  metadata_workbook,
  "Clean metadata"
)

writeData(
  metadata_workbook,
  sheet = "Clean metadata",
  x = clean_metadata
)

addWorksheet(
  metadata_workbook,
  "Group summary"
)

writeData(
  metadata_workbook,
  sheet = "Group summary",
  x = group_summary
)

addWorksheet(
  metadata_workbook,
  "Time summary"
)

writeData(
  metadata_workbook,
  sheet = "Time summary",
  x = time_summary
)

addWorksheet(
  metadata_workbook,
  "Group and time"
)

writeData(
  metadata_workbook,
  sheet = "Group and time",
  x = group_time_summary
)

addWorksheet(
  metadata_workbook,
  "Treatment summary"
)

writeData(
  metadata_workbook,
  sheet = "Treatment summary",
  x = treatment_summary
)

addWorksheet(
  metadata_workbook,
  "Sex summary"
)

writeData(
  metadata_workbook,
  sheet = "Sex summary",
  x = sex_summary
)

# Apply automatic column widths
for (
  sheet_name in names(metadata_workbook)
) {
  setColWidths(
    metadata_workbook,
    sheet = sheet_name,
    cols = 1:20,
    widths = "auto"
  )
}

saveWorkbook(
  metadata_workbook,
  file = file.path(
    tables_directory,
    "GSE287540_metadata_summary.xlsx"
  ),
  overwrite = TRUE
)

# ============================================================
# Final validation
# ============================================================

expected_sample_number <- 128
expected_subject_number <- 61

if (nrow(clean_metadata) != expected_sample_number) {
  warning(
    "The detected sample number is different from the expected 128 samples."
  )
}

if (number_of_subjects != expected_subject_number) {
  warning(
    "The detected subject number is different from the expected 61 subjects."
  )
}

if (duplicated_samples > 0) {
  warning("Duplicated sample identifiers were detected.")
}

cat("\n============================================\n")
cat("GSE287540 metadata preparation completed\n")
cat("============================================\n")

cat("\nExpected results:\n")
cat("Total samples: 128\n")
cat("Unique subjects: 61\n")
cat("ChILI samples: 54\n")
cat("Control samples: 74\n")

cat("\nFiles created:\n")
cat("data/processed/GSE287540_metadata_clean.csv\n")
cat("results/R_objects/GSE287540_metadata_clean.rds\n")
cat("results/tables/GSE287540_metadata_summary.xlsx\n")