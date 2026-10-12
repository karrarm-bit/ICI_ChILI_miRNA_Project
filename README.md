# ICI-Induced Liver Injury and miRNA-Associated Feature Analysis

**Graduation research project | College of Pharmacy, University of Karbala, Iraq**

This repository contains analysis code, public source data, and research outputs for an exploratory bioinformatics study of immune checkpoint inhibitor-related liver injury and treatment response. It provides a documented workflow for academic evaluation, reproducibility, and student training.

## Authors

- Karrar Mahmood Shaker Al-Khfaji — corresponding author
- Shahad Dakhil
- Tuktam Adnan
- Fatima Mohsin

## Study overview

The study reanalyses two independent Gene Expression Omnibus (GEO) datasets:

| Dataset | Molecular material | Analysis |
|---|---|---|
| [GSE287540](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287540) | Peripheral blood mononuclear cell mRNA | Primary comparison of ChILI TOL and control AC samples, with batch adjustment and a supplementary within-subject comparison |
| [GSE274664](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274664) | CD19-positive extracellular vesicle miRNA-associated features | Responder versus non-responder comparison, sample quality control, and matched evaluation of feature definitions |

The workflow includes count preparation, quality control, expression filtering, TMM normalization, edgeR differential expression, functional enrichment, sensitivity analyses, and integration with miRTarBase target records.

The datasets represent different participants, molecular materials, and clinical comparisons. Their integration is intended to generate biological hypotheses; it does not establish direct regulatory relationships between the cohorts.

## Main findings

| Analysis | Deposited result |
|---|---|
| Primary mRNA comparison | 15,977 genes tested; 203 differentially expressed genes, including 144 higher and 59 lower in ChILI |
| Within-subject mRNA comparison | 54 differentially expressed genes; provides internal supporting evidence |
| Principal miRNA-associated feature analysis | 65 hairpin-associated features tested in 12 samples: 8 responders and 4 non-responders |
| Principal significant feature | hsa-miR-26b-5p: log2 fold change = -7.1406; FDR = 0.04012 |
| miR-26b-5p target overlap | 17 targets overlapped with ChILI DEGs; the overlap was not enriched, with Fisher P = 0.71981 and BH-adjusted FDR = 1.00 |

Primary differential-expression selection uses **FDR < 0.05 and absolute log2 fold change >= 1**. In the response comparison, a positive log2 fold change indicates higher abundance in responders; the negative miR-26b-5p estimate indicates higher abundance in non-responders.

### Principal and supplementary feature definitions

The principal miRNA-related interpretation uses the matched **hairpin-associated** analysis. These count-file labels do not, by themselves, demonstrate precursor or mature miRNA abundance.

The earlier combined mature-plus-hairpin analysis tested 67 features and identified four significant candidates. Its outputs and leave-one-out analyses are retained as supplementary evidence. Combining the two feature categories requires verification of their comparability and read-assignment rules. The matched comparison reproduced the stored combined-feature results exactly.

HGF and other target overlaps are candidates for further investigation. Neither target overlap nor database targeting evidence establishes a mechanism in the study cohorts.

## Repository contents

| Location | Contents |
|---|---|
| `scripts/` | Ordered R scripts for data preparation and analysis |
| `data/raw/` | Original count files and GEO series matrices |
| `data/processed/` | Sample definitions, quality-control decisions, differential-expression tables, and feature-review outputs |
| `results/tables/` | Excel workbooks containing quality-control and analysis results |
| `results/figures/` | Selected mRNA figures and principal hairpin-analysis figures |
| `results/R_objects/` | Saved analysis objects and cached target queries |
| `docs/` | Execution instructions, technical review, and original package-version records |
| `ICI_ChILI_miRNA_Project.Rproj` | RStudio project file |
| `RUN_ANALYSIS.R` | Ordered analysis runner |
| `VERIFY_REPOSITORY.R` | Native R syntax and input checks |
| `CITATION.cff` | Author information and citation metadata |

For the script sequence and interpretation of individual outputs, see [Running and review instructions](docs/RUNNING_AND_REVIEW.md).

## Running the analysis

Download or clone the repository into a directory named **`ICI_ChILI_miRNA_Project`**. Open the RStudio project and keep the working directory at the repository root. Use a separate working copy to preserve the deposited outputs.

Install the required packages and perform the initial checks:

```r
source("scripts/00_install_packages.R")
source("VERIFY_REPOSITORY.R")
```

Run the matched feature comparison from a terminal in the repository root:

```bash
Rscript --vanilla RUN_ANALYSIS.R review
```

To run the complete workflow:

```bash
Rscript --vanilla RUN_ANALYSIS.R full
```

The full workflow regenerates intermediate files, including supplementary analyses. Package installation and some enrichment services require internet access. Target retrieval scripts use the deposited cached queries where available.

## Reproducibility and interpretation

Repository checks cover code syntax, count-matrix consistency, sample metadata, multiple-testing correction, and target-overlap calculations. A complete independent end-to-end execution in R has not been verified for this release.

Original installed package versions are recorded in [package_versions.csv](docs/package_versions.csv). The repository does not include a frozen `renv.lock` environment. Package compatibility, annotation updates, and online database changes should be considered when comparing a new run with the deposited results.

The study is exploratory. The paired mRNA analysis provides internal support, and the miRNA feature-definition review is a sensitivity assessment. Neither constitutes independent experimental validation. Findings require confirmation in appropriate biological studies.

## Citation and source attribution

Please cite the archived repository release using the author information in [CITATION.cff](CITATION.cff), and acknowledge the original GEO datasets and their associated publications. Once available, the Zenodo DOI should be included in citations to the archived version.

Public source data retain their original attribution and applicable terms. No software license has been specified for this release; contact the corresponding author regarding code reuse.
