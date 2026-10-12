# ICI-induced liver injury and miRNA-associated feature analysis

Review and teaching release prepared from the supplied graduation project.
This repository reanalyses two independent public cohorts: **GSE287540**
(PBMC mRNA in checkpoint inhibitor-induced liver injury) and **GSE274664**
(CD19-positive extracellular-vesicle miRNA-associated features in treatment
responders and non-responders). It supports inspection and rerunning; it is
not evidence of direct miRNA regulation across these separate cohorts.

## Start here

- Evaluators: [review findings](docs/REVIEW_REPORT.md), then
  [reviewed manuscript draft](reports/Graduation_Project_Reviewed_Draft.docx).
- Students: [script order and reproducibility](docs/RUNNING_AND_REVIEW.md).
- Upload and DOI instructions in Arabic: [GitHub and Zenodo guide](docs/UPLOAD_AND_DOI_AR.md).
- File integrity and provenance: [manifest](FILE_MANIFEST.csv) and [checksums](SHA256SUMS.txt).

## Principal and supplementary results

The stored primary mRNA model tests **15,977 genes**, identifying **203 DEGs**
(144 higher and 59 lower in ChILI). The paired comparison identifies 54 DEGs;
this is internal support, not independent cohort validation.

The later matched feature review is the principal miRNA-related result:
**65 hairpin-associated features**, **12 samples** (8 responders, 4
non-responders), and **one feature** meeting FDR < 0.05 and |log2FC| >= 1:
**hsa-miR-26b-5p**, log2FC -7.1406, FDR 0.04012. Positive logFC means higher
abundance in responders. Hairpin labels do not establish precursor abundance.

The earlier pooled mature + hairpin analysis tests 67 features and reports
four significant candidates. Its CSVs, workbooks and leave-one-out outputs
are **supplementary historical results**, not four confirmed principal findings.
The saved matched comparison reproduces that earlier model exactly.

The 17 miR-26b-5p target/DEG overlaps are **not enriched** (Fisher P = 0.71981;
BH FDR = 1). HGF is a hypothesis-generating overlap, not an established mechanism.

## Quick run

Extract/clone into a folder named `ICI_ChILI_miRNA_Project`, open the `.Rproj`,
and keep the working directory at the repository root. Work in a separate copy
to preserve the deposited result snapshot.

```r
source("scripts/00_install_packages.R")
source("VERIFY_REPOSITORY.R")
```

From a terminal in the same folder:

```bash
Rscript --vanilla RUN_ANALYSIS.R review
Rscript --vanilla RUN_ANALYSIS.R full
```

`review` reruns the matched feature comparison from supplied raw data and
checkpoints. `full` rebuilds the main analysis, supplementary pooled analysis,
target integration and both feature reviews. Enrichment may require online
services. Target scripts use the supplied cached queries where available.
The package installer requires internet access.

## What was verified

All retained R code was syntax checked with a tree-sitter R parser. Supplied
counts, metadata, stored results, FDR correction and target-overlap statistics
were checked independently in Python. **R/edgeR was not executed in the
preparation environment**; native R preflight and a full rerun remain to be
completed on the evaluator's or author's machine. See the review report for
the exact checks, fixes and limitations.

No installed Windows package library, RStudio cache, workspace dump, duplicate
PNG/PDF figure pairs or superseded manuscript versions are included.
The original uploaded project remains the provenance source.

## Attribution, release and reuse

The manuscript still contains student/supervisor placeholders. Complete the
actual author list and review the draft before public release. Do not infer
authorship from the person preparing the archive. A citation template is
provided as `CITATION.cff.example`; rename it to `CITATION.cff` only after
filling verified authors and repository details. No DOI has yet been minted.
Choose a code license as the rights holders; none is assigned by this preparation.
Public source data retain their original source attribution and terms.

Source datasets: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287540
and https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274664 .
