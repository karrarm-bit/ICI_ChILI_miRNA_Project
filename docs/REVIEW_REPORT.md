# Code and deposited-evidence review

## Outcome

The release removes local RStudio state and installed Windows libraries,
retains the raw inputs, completed analysis code, selected review evidence,
cached target queries, reviewed draft and teaching instructions.

## Changes to code

1. **Script 17A S4 guard corrected**: original `isS4(target_query) || ...`
   rejects valid S4 query objects. Changed to `!isS4(target_query) || ...`.
   The default pipeline instead uses the deposited query through script 17B.
2. **Scripts 22 and 23 made portable**: hardcoded `D:/ICI_ChILI_miRNA_Project`
   replaced by the current repository root. Statistical settings unchanged.
3. **Canonical step 11 selected**: original 11B and `11_2 prepare...` are
   byte-identical; retained one as `11_prepare_combined_miRNA_counts.R`.
4. **Dependency installer completed** and a process-isolated runner added.
   These are preparation utilities, not original analysis provenance.
5. RStudio project now disables workspace/history restoration and saving.

## Checks passed during preparation

- Tree-sitter R syntax parsing of all 29 retained R files (not native R execution).
- Four raw gzip files decompressed completely; retained XLSX/DOCX containers
  passed ZIP integrity checks.
- Independently reconstructed combined counts from the original count file:
  **275 features x 13 samples**, exact values and row/column order agree with
  deposited RDS. Counts are finite, nonnegative integers.
- Final metadata: unique sample IDs, **8 responders + 4 non-responders**;
  NR5 (`Hg19_699_CD19`) excluded as recorded.
- Primary mRNA CSV: **15,977 tests; 203 DEGs; 144 up and 59 down** under
  FDR < 0.05 and |logFC| >= 1. BH FDR independently recomputed and agrees.
- Principal stored hairpin review: **65 tests; one primary significant feature;
  14 nominal candidates**. BH FDR independently recomputed and agrees.
- Archived matched reproduction table reports identical baseline features,
  zero maximum logFC difference and zero maximum FDR difference. This checks
  deposited evidence and is not a new edgeR rerun performed here.
- Independently recomputed one-sided Fisher probabilities and BH adjustment
  for the three deposited exact-name target tests; agree with the CSV.
- Original copied files have unchanged SHA256; modified files are explicitly
  marked in FILE_MANIFEST.csv. Every deposited file is below 25 MiB.

## Interpretation and remaining verification

Only miR-26b-5p is primary significant in the matched hairpin review.
miR-509-3p is not retained for hairpin testing; miR-196b-5p and miR-375-3p
have nominal, not FDR-significant, support. Pooled mature/hairpin comparability
and non-overlapping read assignment were not established. The datasets have
different subjects and contrasts, so overlap does not demonstrate an inverse
expression relationship or causal regulation. miR-375 alias targets remain
provisional. miR-26b's 17 overlaps are not enriched (P=0.71981, BH FDR=1).

**R and edgeR were not available in the preparation runtime.** No full native
analytical rerun is claimed. Run `VERIFY_REPOSITORY.R`, then `RUN_ANALYSIS.R
review`, and preferably `full`, in a fresh copy on a machine with R/Bioconductor.
Cached query extraction, package compatibility and online enrichment execution
therefore still require native R confirmation. Native verification utilities
themselves passed the R grammar check here.

No renv lockfile was supplied. Package versions are recorded from the original
installed files but are not a restoration guarantee. The manuscript remains
a supplied reviewed draft with author placeholders and supervisor review
requirements. Its layout and references were not newly validated here.

## Exclusions

Excluded installed libraries, workspace/history, RStudio state, redundant
step 11 copy, mature-only preparation, earlier QC script 12, generic target
script 17 without completed output provenance, superseded Word generators,
old drafts, duplicate figure formats and intermediate objects regenerable by
the ordered full pipeline. The uploaded original archive was not modified.
