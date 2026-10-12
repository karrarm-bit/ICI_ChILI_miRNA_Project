# Run from the root folder ICI_ChILI_miRNA_Project.
# Each stage gets its own R process because original scripts clear the workspace.
# This regenerates outputs: run in a separate working copy of the release.
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "review"
if (!mode %in% c("review", "full")) stop("Usage: Rscript --vanilla RUN_ANALYSIS.R [review|full]")
if (!dir.exists("scripts") || !dir.exists("data/raw")) stop("Run from the repository root.")
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
run_stage <- function(script) {
  cat("\nRunning:", script, "\n")
  status <- system2(rscript, c("--vanilla", shQuote(file.path("scripts", script))))
  if (status != 0L) stop("Stage failed: ", script, "; later stages were not run.")
}
review <- c("23_GSE274664_matched_feature_review.R")
full <- c(
  "01_create_project_folders.R", "02_download_GSE287540.R",
  "03_prepare_GSE287540_metadata.R", "04_prepare_GSE287540_count_matrix.R",
  "05_GSE287540_quality_control.R", "06_review_QC_and_define_samples.R",
  "07_GSE287540_differential_expression.R", "08_GSE287540_paired_validation.R",
  "09_GSE287540_functional_enrichment.R", "10_download_GSE274664_miRNA.R",
  "11_prepare_combined_miRNA_counts.R", "13_GSE274664_miRNA_quality_control_final.R",
  "14_review_miRNA_QC_and_define_samples.R", "15_GSE274664_miRNA_differential_expression.R",
  "16_GSE274664_miRNA_sensitivity_analysis.R", "17B_extract_saved_miR509_targets.R",
  "17C_download_miR26b_targets.R", "17D_download_miR196b_targets.R",
  "17F_download_miR375_alias_targets.R", "18_integrate_miRNA_targets_with_mRNA_DEGs.R",
  "19_review_target_evidence.R", "20_test_miRNA_target_overlap_enrichment.R",
  "22_GSE274664_hairpin_only_review.R", "23_GSE274664_matched_feature_review.R")
for (script in if (mode == "full") full else review) run_stage(script)
cat("Completed:", mode, "mode. Review regenerated results before use.\n")
