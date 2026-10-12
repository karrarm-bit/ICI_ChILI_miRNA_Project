# Native R syntax and input preflight; no analysis and no network calls.
if (!file.exists("RUN_ANALYSIS.R")) stop("Run from the repository root.")
scripts <- c("RUN_ANALYSIS.R", "VERIFY_REPOSITORY.R", list.files("scripts", "\\.R$", full.names=TRUE))
for (script in scripts) {
  parse(file=script, keep.source=FALSE)
  cat("Syntax OK:", script, "\n")
}
required <- c("data/raw/GSE274664_CD19_EV_miRNA_readcounts.txt.gz",
 "data/raw/GSE287540_ChILI_bulk_readcount.tsv.gz",
 "results/R_objects/GSE274664_mature_miRNA_count_matrix.rds",
 "results/R_objects/GSE274664_final_miRNA_metadata.rds",
 "results/R_objects/GSE274664_miRNA_DE_results.rds")
if (any(!file.exists(required))) stop("Missing input files: ", paste(required[!file.exists(required)], collapse=", "))
y <- readRDS(required[3]); meta <- readRDS(required[4]); old <- readRDS(required[5])
stopifnot(identical(dim(y), c(275L, 13L)), nrow(meta)==12L, nrow(old)==67L,
 all(meta$Sample_ID %in% colnames(y)), !anyDuplicated(meta$Sample_ID),
 all(is.finite(y)), all(y >= 0), all(y == round(y)))
cat("Preflight passed. This does not constitute a complete analytical rerun.\n")
