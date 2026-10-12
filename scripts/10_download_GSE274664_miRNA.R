# ============================================================
# Script 10: Download GSE274664 microRNA Data
# Project: ICI-Induced Liver Injury and miRNA Analysis
# ============================================================

# Clear the R environment
rm(list = ls())

# Increase the download timeout
options(
  stringsAsFactors = FALSE,
  timeout = 1800
)

# ============================================================
# Define project directories
# ============================================================

project_directory <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

raw_data_directory <- file.path(
  project_directory,
  "data",
  "raw"
)

logs_directory <- file.path(
  project_directory,
  "logs"
)

# Check required folders
if (!dir.exists(raw_data_directory)) {
  stop(
    "The data/raw folder does not exist. ",
    "Run 01_create_project_folders.R first."
  )
}

if (!dir.exists(logs_directory)) {
  
  dir.create(
    logs_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ============================================================
# Define GEO download URLs
# ============================================================

mirna_counts_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  "GSE274nnn/GSE274664/suppl/",
  "GSE274664_Processed_data_CD19_EVs_readcounts.txt.gz"
)

metadata_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  "GSE274nnn/GSE274664/matrix/",
  "GSE274664_series_matrix.txt.gz"
)

# ============================================================
# Define destination files
# ============================================================

mirna_counts_file <- file.path(
  raw_data_directory,
  "GSE274664_CD19_EV_miRNA_readcounts.txt.gz"
)

metadata_file <- file.path(
  raw_data_directory,
  "GSE274664_series_matrix.txt.gz"
)

# ============================================================
# Download function
# ============================================================

download_geo_file <- function(
    file_url,
    destination_file,
    file_description) {
  
  file_name <- basename(
    destination_file
  )
  
  cat("\n--------------------------------------------\n")
  cat("File:", file_description, "\n")
  cat("Destination:", file_name, "\n")
  
  if (file.exists(destination_file)) {
    
    existing_size <- file.info(
      destination_file
    )$size
    
    if (
      !is.na(existing_size) &&
      existing_size > 0
    ) {
      
      cat(
        "Status: file already exists and will not be downloaded again.\n"
      )
      
      return(TRUE)
    }
  }
  
  cat("Status: downloading...\n")
  
  download_result <- tryCatch(
    
    {
      download.file(
        url = file_url,
        destfile = destination_file,
        mode = "wb",
        method = "libcurl",
        quiet = FALSE
      )
      
      TRUE
    },
    
    error = function(error_message) {
      
      cat(
        "Download error:",
        error_message$message,
        "\n"
      )
      
      FALSE
    }
  )
  
  download_success <- download_result &&
    file.exists(destination_file) &&
    file.info(destination_file)$size > 0
  
  if (download_success) {
    
    downloaded_size_mb <- round(
      file.info(destination_file)$size /
        1024^2,
      4
    )
    
    cat("Status: download completed.\n")
    cat(
      "Downloaded size:",
      downloaded_size_mb,
      "MB\n"
    )
    
  } else {
    
    if (file.exists(destination_file)) {
      file.remove(destination_file)
    }
    
    stop(
      "The download failed for: ",
      file_description
    )
  }
  
  return(download_success)
}

# ============================================================
# Download the CD19-positive EV miRNA count matrix
# ============================================================

mirna_counts_downloaded <- download_geo_file(
  file_url = mirna_counts_url,
  destination_file = mirna_counts_file,
  file_description = paste(
    "CD19-positive extracellular-vesicle",
    "microRNA count matrix"
  )
)

# ============================================================
# Download the GEO sample metadata
# ============================================================

metadata_downloaded <- download_geo_file(
  file_url = metadata_url,
  destination_file = metadata_file,
  file_description = "GSE274664 sample metadata"
)

# ============================================================
# Validate the downloaded gzip files
# ============================================================

test_gzip_file <- function(file_path) {
  
  validation_result <- tryCatch(
    
    {
      connection_object <- gzfile(
        file_path,
        open = "rt"
      )
      
      first_line <- readLines(
        connection_object,
        n = 1,
        warn = FALSE
      )
      
      close(
        connection_object
      )
      
      length(first_line) == 1 &&
        nchar(first_line) > 0
    },
    
    error = function(error_message) {
      
      FALSE
    }
  )
  
  return(validation_result)
}

mirna_gzip_valid <- test_gzip_file(
  mirna_counts_file
)

metadata_gzip_valid <- test_gzip_file(
  metadata_file
)

# ============================================================
# Create download manifest
# ============================================================

download_manifest <- data.frame(
  Dataset = c(
    "GSE274664",
    "GSE274664"
  ),
  
  Data_type = c(
    "CD19-positive EV microRNA raw counts",
    "GEO sample metadata"
  ),
  
  File_name = c(
    basename(mirna_counts_file),
    basename(metadata_file)
  ),
  
  File_size_MB = round(
    c(
      file.info(mirna_counts_file)$size,
      file.info(metadata_file)$size
    ) / 1024^2,
    4
  ),
  
  Download_status = c(
    mirna_counts_downloaded,
    metadata_downloaded
  ),
  
  Valid_gzip_file = c(
    mirna_gzip_valid,
    metadata_gzip_valid
  ),
  
  stringsAsFactors = FALSE
)

cat("\n============================================\n")
cat("GSE274664 download manifest\n")
cat("============================================\n")

print(
  download_manifest
)

# Save the download manifest
write.csv(
  download_manifest,
  file = file.path(
    logs_directory,
    "GSE274664_miRNA_download_manifest.csv"
  ),
  row.names = FALSE
)

# ============================================================
# Final validation
# ============================================================

all_downloaded <- all(
  download_manifest$Download_status
)

all_valid <- all(
  download_manifest$Valid_gzip_file
)

all_nonempty <- all(
  download_manifest$File_size_MB > 0
)

if (
  !all_downloaded ||
  !all_valid ||
  !all_nonempty
) {
  
  stop(
    "One or more GSE274664 files failed validation."
  )
}

# ============================================================
# Final summary
# ============================================================

cat("\n============================================\n")
cat("GSE274664 miRNA data downloaded successfully\n")
cat("============================================\n")

cat("\nFiles saved:\n")

cat(
  "data/raw/",
  basename(mirna_counts_file),
  "\n",
  sep = ""
)

cat(
  "data/raw/",
  basename(metadata_file),
  "\n",
  sep = ""
)

cat("\nExpected biological comparison:\n")
cat("ICI responders: 8 samples\n")
cat("ICI non-responders: 5 samples\n")
cat("Biological material: CD19-positive extracellular vesicles\n")
cat("Data type: small-RNA sequencing read counts\n")