# ============================================================
# Script 02: Download GSE287540 Data
# Project: ICI-Induced Liver Injury and miRNA Analysis
# ============================================================

# Clear the R environment
rm(list = ls())

# Increase the download timeout
options(timeout = 1800)

# Define the project directory
project_directory <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

# Define the raw data directory
raw_data_directory <- file.path(
  project_directory,
  "data",
  "raw"
)

# Check whether the raw data folder exists
if (!dir.exists(raw_data_directory)) {
  stop(
    "The data/raw folder does not exist. ",
    "Run 01_create_project_folders.R first."
  )
}

# ============================================================
# File 1: Raw gene-count matrix
# ============================================================

counts_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  "GSE287nnn/GSE287540/suppl/",
  "GSE287540_ChILI_bulk_readcount.tsv.gz"
)

counts_file <- file.path(
  raw_data_directory,
  "GSE287540_ChILI_bulk_readcount.tsv.gz"
)

# ============================================================
# File 2: GEO series matrix containing sample metadata
# ============================================================

metadata_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  "GSE287nnn/GSE287540/matrix/",
  "GSE287540_series_matrix.txt.gz"
)

metadata_file <- file.path(
  raw_data_directory,
  "GSE287540_series_matrix.txt.gz"
)

# ============================================================
# Download function
# ============================================================

download_geo_file <- function(file_url, destination_file) {
  
  file_name <- basename(destination_file)
  
  if (file.exists(destination_file)) {
    
    current_size <- file.info(destination_file)$size
    
    if (!is.na(current_size) && current_size > 0) {
      
      cat(
        "File already exists:",
        file_name,
        "\n"
      )
      
      return(TRUE)
    }
  }
  
  cat(
    "\nDownloading:",
    file_name,
    "\n"
  )
  
  download.file(
    url = file_url,
    destfile = destination_file,
    mode = "wb",
    method = "libcurl",
    quiet = FALSE
  )
  
  download_success <- file.exists(destination_file) &&
    file.info(destination_file)$size > 0
  
  if (!download_success) {
    stop(
      "Download failed for: ",
      file_name
    )
  }
  
  cat(
    "Download completed:",
    file_name,
    "\n"
  )
  
  return(TRUE)
}

# ============================================================
# Download the files
# ============================================================

counts_downloaded <- download_geo_file(
  file_url = counts_url,
  destination_file = counts_file
)

metadata_downloaded <- download_geo_file(
  file_url = metadata_url,
  destination_file = metadata_file
)

# ============================================================
# Verify downloaded files
# ============================================================

download_manifest <- data.frame(
  Dataset = c(
    "GSE287540",
    "GSE287540"
  ),
  
  File_type = c(
    "Raw gene-count matrix",
    "GEO sample metadata"
  ),
  
  File_name = c(
    basename(counts_file),
    basename(metadata_file)
  ),
  
  File_size_MB = round(
    c(
      file.info(counts_file)$size,
      file.info(metadata_file)$size
    ) / 1024^2,
    3
  ),
  
  Download_status = c(
    counts_downloaded,
    metadata_downloaded
  ),
  
  stringsAsFactors = FALSE
)

# Display the download manifest
cat("\nDownloaded files:\n")
print(download_manifest)

# Save the download manifest
write.csv(
  download_manifest,
  file = file.path(
    project_directory,
    "logs",
    "GSE287540_download_manifest.csv"
  ),
  row.names = FALSE
)

# Final confirmation
if (
  all(download_manifest$Download_status) &&
  all(download_manifest$File_size_MB > 0)
) {
  
  cat("\n============================================\n")
  cat("GSE287540 files were downloaded successfully\n")
  cat("============================================\n")
  
  cat("\nFiles were saved in:\n")
  cat(raw_data_directory, "\n")
  
} else {
  
  stop("One or more GSE287540 files are missing.")
}