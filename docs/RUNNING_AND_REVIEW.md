# Execution and evaluator map

Run in a fresh copy of the repository root. Use `--vanilla` to avoid hidden
workspace state. Each pipeline stage runs in a separate process, so `rm(list=ls())`
inside original scripts cannot erase the pipeline controller.

| Stages | Purpose | Evaluator evidence |
|---|---|---|
| 00, 01 | Dependencies and folders | `docs/package_versions.csv`; session information |
| 02–04 | mRNA acquisition, metadata, count alignment | Raw counts and GEO matrix; preparation workbooks |
| 05–06 | mRNA QC and batch-aware primary sample selection | QC workbooks; primary sample CSV; PCA PDFs |
| 07 | Primary mRNA DE | Full 15,977-row CSV; DE workbook |
| 08 | Within-subject comparison | Paired workbook; primary/paired comparison CSV |
| 09 | GO, KEGG and Hallmark enrichment | Enrichment workbook and selected PDFs |
| 10, 11 | Original miRNA data; combined feature baseline | Raw file; 275 x 13 stored matrix; preparation workbook |
| 13–16 | Pooled-feature QC, final samples, DE, leave-one-out | QC decisions, 12-sample metadata, supplementary 67-feature CSV and workbooks |
| 17B, 17C, 17D, 17F | Cached target extraction | Saved query objects; source provenance in original scripts |
| 18–20 | Target overlap, evidence review, background-aware Fisher tests | Overlap CSV; evidence tables; enrichment CSV |
| 22 | Hairpin-only sensitivity review with its own settings | Writes a separate local review folder; does not replace baseline |
| 23 | Matched-settings combined/hairpin comparison | Saved 65-feature hairpin result; comparison and reproduction CSVs |

## Decisions affecting execution

- `11_prepare_combined_miRNA_counts.R` is the original script 11B, identical
  to the file called `11_2 prepare_GSE274664_miRNA.R`. The original mature-only
  script 11 is superseded and excluded. Its output filename contains `mature`
  for compatibility but the stored matrix is **combined mature + hairpin**.
- Script 13 replaces the earlier script 12 for final pre-exclusion QC.
- The generic original script 17 was not the completed target route; use
  the separate saved queries. Script 17A is a corrected optional fresh
  miR-509 query. Script 17E is optional identifier diagnosis and may stop if
  no exact miR-375-3p records exist. The historical alias query in 17F remains
  provisional. Neither 17A nor 17E is in the standard cached pipeline.
- Script 22 uses a sensitivity definition; script 23 matches the original
  filtering and statistical settings and is the principal comparison.
- The full run regenerates intermediate RDS/CSV files omitted from the
  deposited release to keep the upload below 100 files. The deposited review
  mode has all required inputs without these additional intermediates.
- Original document generators are omitted: they generate superseded claims
  or require local Word selections. The supplied reviewed draft is included
  unchanged and has not been declared publication ready.

## Reproducibility limits

The original archive has no `renv.lock`. Its Windows installed libraries are
not a portable environment. `docs/package_versions.csv` records package versions
from their original DESCRIPTION files; the installer retrieves currently
available compatible versions and is not a frozen environment. The archived
session file is from the matched review (edgeR 4.10.5 in the supplied record).
After a successful native R rerun, record `sessionInfo()` and make a reviewed
`renv.lock` if exact environment restoration is required. Do not copy Windows
package binaries to another OS.

Enrichment databases, annotations and fresh target services can change.
The archived outputs and cached target queries are the deposited evidence.
The full run should be compared with this snapshot before a new release.
The manuscript's Chapter One, references, time-point terminology and layout
still need supervisor review; author placeholders remain.
