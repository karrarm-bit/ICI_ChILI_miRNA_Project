# Install dependencies; run once from the repository root.
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 1800)
cran <- c("data.table", "ggplot2", "ggrepel", "pheatmap", "RColorBrewer",
          "openxlsx", "patchwork", "msigdbr", "statmod")
bioc <- c("GEOquery", "edgeR", "limma", "AnnotationDbi", "org.Hs.eg.db",
          "clusterProfiler", "enrichplot", "fgsea", "multiMiR")
missing <- function(pkgs) pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing(cran))) install.packages(missing(cran))
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
if (length(missing(bioc))) BiocManager::install(missing(bioc), ask=FALSE, update=FALSE)
absent <- missing(c(cran, bioc))
if (length(absent)) stop("Packages remain unavailable: ", paste(absent, collapse=", "))
cat("Dependencies available. Versions are not pinned; see docs/package_versions.csv.\n")
