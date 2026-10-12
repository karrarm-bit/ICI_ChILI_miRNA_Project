# ============================================================
# Script 01: Create Project Folders
# Project: ICI-Induced Liver Injury and miRNA Analysis
# ============================================================

# Clear the R environment
rm(list = ls())

# Define the main project directory
project_directory <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

# Display the current project directory
cat("\nCurrent project directory:\n")
cat(project_directory, "\n\n")

# Define the required project folders
project_folders <- c(
  "data",
  "data/raw",
  "data/processed",
  "results",
  "results/tables",
  "results/figures",
  "results/network",
  "results/R_objects",
  "logs",
  "references",
  "manuscript"
)

# Create the folders
for (folder_name in project_folders) {
  
  folder_path <- file.path(
    project_directory,
    folder_name
  )
  
  if (!dir.exists(folder_path)) {
    
    dir.create(
      folder_path,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    cat("Created:", folder_name, "\n")
    
  } else {
    
    cat("Already exists:", folder_name, "\n")
  }
}

# Check whether all folders were created
folder_status <- data.frame(
  Folder = project_folders,
  Exists = dir.exists(
    file.path(
      project_directory,
      project_folders
    )
  )
)

# Display the folder status
cat("\nFolder creation status:\n")
print(folder_status)

# Save the folder status as a CSV file
write.csv(
  folder_status,
  file = file.path(
    project_directory,
    "logs",
    "project_folder_status.csv"
  ),
  row.names = FALSE
)

# Final confirmation message
if (all(folder_status$Exists)) {
  
  cat("\n============================================\n")
  cat("All project folders were created successfully\n")
  cat("============================================\n")
  
} else {
  
  warning("Some project folders could not be created.")
}