# ==============================================================================
# Script 19: Review experimental evidence for overlapping miRNA target genes
#
# This script performs no online queries.
# Evidence flags describe reported assays, not clinical validation.
# Provisional miR-375 alias results remain separate.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Define inputs and output directories
# ------------------------------------------------------------------------------

integration_file <- paste0(
  "results/R_objects/",
  "GSE274664_GSE287540_target_integration.rds"
)

evidence_directory <- "results/R_objects/miRNA_targets"

evidence_files <- c(
  "hsa-miR-509-3p" = file.path(
    evidence_directory,
    "hsa_miR_509_3p_miRTarBase_evidence.rds"
  ),
  "hsa-miR-26b-5p" = file.path(
    evidence_directory,
    "hsa_miR_26b_5p_miRTarBase_evidence.rds"
  ),
  "hsa-miR-196b-5p" = file.path(
    evidence_directory,
    "hsa_miR_196b_5p_miRTarBase_evidence.rds"
  )
)

provisional_evidence_file <- file.path(
  evidence_directory,
  "hsa_miR_375_3p_alias_miRTarBase_evidence.rds"
)

csv_directory <- "data/processed/target_evidence_review"
rds_directory <- "results/R_objects"
table_directory <- "results/tables"

dir.create(csv_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  integration_file,
  unname(evidence_files),
  provisional_evidence_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    "Missing input files:\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 2. Load the integration results
# ------------------------------------------------------------------------------

cat("Loading exploratory target-integration results...\n")

integration_results <- readRDS(integration_file)

required_tables <- c(
  "Main_DEG_overlap",
  "Provisional_375_overlap"
)

if (!all(required_tables %in% names(integration_results))) {
  stop("Required integration tables were not found.")
}

main_overlap <- integration_results$Main_DEG_overlap
provisional_overlap <- integration_results$Provisional_375_overlap

# ------------------------------------------------------------------------------
# 3. Utility functions
# ------------------------------------------------------------------------------

collapse_values <- function(values) {
  
  values <- trimws(as.character(values))
  values <- values[!is.na(values) & nzchar(values)]
  
  if (length(values) == 0L) {
    return(NA_character_)
  }
  
  paste(sort(unique(values)), collapse = "; ")
}

extract_pmids <- function(values) {
  
  values <- as.character(values)
  values <- values[!is.na(values) & nzchar(trimws(values))]
  
  if (length(values) == 0L) {
    return(character(0))
  }
  
  # Extract numeric PMID tokens while retaining the original field elsewhere.
  tokens <- unlist(
    regmatches(values, gregexpr("[0-9]+", values))
  )
  
  sort(unique(tokens[nzchar(tokens)]))
}

load_evidence <- function(file, study_miRNA, mapping_status) {
  
  evidence <- readRDS(file)
  
  if (!is.data.frame(evidence) ||
      !"target_symbol" %in% names(evidence)) {
    stop("Unexpected evidence structure in: ", file)
  }
  
  optional_columns <- c(
    "experiment",
    "support_type",
    "pubmed_id",
    "mature_mirna_id",
    "database"
  )
  
  for (column in optional_columns) {
    if (!column %in% names(evidence)) {
      evidence[[column]] <- rep(NA_character_, nrow(evidence))
    }
  }
  
  evidence$Study_miRNA <- rep(study_miRNA, nrow(evidence))
  evidence$Mapping_status <- rep(mapping_status, nrow(evidence))
  
  evidence$Match_symbol <- toupper(
    trimws(as.character(evidence$target_symbol))
  )
  
  # Remove completely identical rows before counting evidence records.
  evidence <- unique(evidence)
  rownames(evidence) <- NULL
  
  evidence
}

review_overlap <- function(overlap, evidence) {
  
  pairs <- unique(overlap[, c("miRNA", "Match_symbol"), drop = FALSE])
  
  if (nrow(pairs) == 0L) {
    return(data.frame(
      miRNA = character(),
      Target_symbol = character(),
      Mapping_status = character(),
      Deduplicated_evidence_records = integer(),
      Unique_PMIDs = integer(),
      Reporter_assay_reported = logical(),
      Western_blot_reported = logical(),
      qPCR_reported = logical(),
      High_throughput_binding_reported = logical(),
      Evidence_group = character(),
      Experiments = character(),
      Support_types = character(),
      PubMed_IDs = character(),
      stringsAsFactors = FALSE
    ))
  }
  
  summaries <- lapply(
    seq_len(nrow(pairs)),
    function(index) {
      
      current_miRNA <- pairs$miRNA[index]
      current_gene <- pairs$Match_symbol[index]
      
      selected_rows <- (
        !is.na(evidence$Study_miRNA) &
          !is.na(evidence$Match_symbol) &
          evidence$Study_miRNA == current_miRNA &
          evidence$Match_symbol == current_gene
      )
      
      current_evidence <- evidence[
        selected_rows,
        ,
        drop = FALSE
      ]
      
      experiments <- as.character(current_evidence$experiment)
      experiments[is.na(experiments)] <- ""
      
      reporter_flag <- any(
        grepl("luciferase|reporter", experiments, ignore.case = TRUE)
      )
      
      western_flag <- any(
        grepl("western", experiments, ignore.case = TRUE)
      )
      
      qpcr_flag <- any(
        grepl(
          "qRT.?PCR|RT.?qPCR|qPCR|real.?time.?PCR",
          experiments,
          ignore.case = TRUE
        )
      )
      
      binding_flag <- any(
        grepl("CLIP|CLASH", experiments, ignore.case = TRUE)
      )
      
      has_assay_information <- any(nzchar(trimws(experiments)))
      
      evidence_group <- if (reporter_flag) {
        "Reporter_assay_reported"
      } else if (western_flag || qpcr_flag) {
        "Protein_or_RNA_measurement_without_reporter"
      } else if (binding_flag) {
        "High_throughput_binding_without_reporter"
      } else if (has_assay_information) {
        "Other_reported_assay"
      } else {
        "Assay_information_missing"
      }
      
      pmids <- extract_pmids(current_evidence$pubmed_id)
      
      data.frame(
        miRNA = current_miRNA,
        Target_symbol = current_gene,
        Mapping_status = collapse_values(
          current_evidence$Mapping_status
        ),
        Deduplicated_evidence_records = nrow(current_evidence),
        Unique_PMIDs = length(pmids),
        Reporter_assay_reported = reporter_flag,
        Western_blot_reported = western_flag,
        qPCR_reported = qpcr_flag,
        High_throughput_binding_reported = binding_flag,
        Evidence_group = evidence_group,
        Experiments = collapse_values(
          current_evidence$experiment
        ),
        Support_types = collapse_values(
          current_evidence$support_type
        ),
        PubMed_IDs = if (length(pmids) > 0L) {
          paste(pmids, collapse = "; ")
        } else {
          NA_character_
        },
        stringsAsFactors = FALSE
      )
    }
  )
  
  do.call(rbind, summaries)
}

attach_mrna_results <- function(summary_table, overlap) {
  
  mrna_columns <- c(
    "miRNA",
    "Match_symbol",
    "mRNA_Gene_symbol",
    "mRNA_Ensembl_ID",
    "mRNA_logFC",
    "mRNA_FDR",
    "ChILI_direction"
  )
  
  if (!all(mrna_columns %in% names(overlap))) {
    stop("Required mRNA overlap columns were not found.")
  }
  
  mrna_information <- unique(
    overlap[, mrna_columns, drop = FALSE]
  )
  
  merged <- merge(
    summary_table,
    mrna_information,
    by.x = c("miRNA", "Target_symbol"),
    by.y = c("miRNA", "Match_symbol"),
    all.x = TRUE,
    sort = FALSE
  )
  
  merged <- merged[
    order(
      !merged$Reporter_assay_reported,
      -merged$Unique_PMIDs,
      merged$mRNA_FDR,
      merged$Target_symbol
    ),
    ,
    drop = FALSE
  ]
  
  rownames(merged) <- NULL
  merged
}

get_overlapping_evidence <- function(evidence, overlap) {
  
  evidence_keys <- paste(
    evidence$Study_miRNA,
    evidence$Match_symbol,
    sep = "|"
  )
  
  overlap_keys <- paste(
    overlap$miRNA,
    overlap$Match_symbol,
    sep = "|"
  )
  
  evidence[
    evidence_keys %in% overlap_keys,
    ,
    drop = FALSE
  ]
}

# ------------------------------------------------------------------------------
# 4. Load the evidence tables
# ------------------------------------------------------------------------------

main_evidence_list <- lapply(
  names(evidence_files),
  function(miRNA_name) {
    load_evidence(
      evidence_files[[miRNA_name]],
      miRNA_name,
      "Exact_name_query"
    )
  }
)

# Use shared columns so differing optional fields do not prevent combination.
shared_columns <- Reduce(
  intersect,
  lapply(main_evidence_list, names)
)

main_evidence <- do.call(
  rbind,
  lapply(
    main_evidence_list,
    function(table) table[, shared_columns, drop = FALSE]
  )
)

provisional_evidence <- load_evidence(
  provisional_evidence_file,
  "hsa-miR-375-3p",
  "Alias_pending_verification"
)

# ------------------------------------------------------------------------------
# 5. Review evidence for the overlapping targets
# ------------------------------------------------------------------------------

main_summary <- review_overlap(
  main_overlap,
  main_evidence
)

main_summary <- attach_mrna_results(
  main_summary,
  main_overlap
)

provisional_summary <- review_overlap(
  provisional_overlap,
  provisional_evidence
)

provisional_summary <- attach_mrna_results(
  provisional_summary,
  provisional_overlap
)

reporter_supported <- main_summary[
  main_summary$Reporter_assay_reported,
  ,
  drop = FALSE
]

main_evidence_detail <- get_overlapping_evidence(
  main_evidence,
  main_overlap
)

provisional_evidence_detail <- get_overlapping_evidence(
  provisional_evidence,
  provisional_overlap
)

# ------------------------------------------------------------------------------
# 6. Save all review outputs
# ------------------------------------------------------------------------------

output_tables <- list(
  Main_evidence_summary = main_summary,
  Reporter_assay_targets = reporter_supported,
  Main_evidence_detail = main_evidence_detail,
  Provisional_375_summary = provisional_summary,
  Provisional_375_detail = provisional_evidence_detail
)

saveRDS(
  output_tables,
  file.path(
    rds_directory,
    "GSE274664_GSE287540_target_evidence_review.rds"
  )
)

for (table_name in names(output_tables)) {
  write.csv(
    output_tables[[table_name]],
    file.path(
      csv_directory,
      paste0(table_name, ".csv")
    ),
    row.names = FALSE,
    na = ""
  )
}

if (requireNamespace("openxlsx", quietly = TRUE)) {
  openxlsx::write.xlsx(
    output_tables,
    file = file.path(
      table_directory,
      "GSE274664_GSE287540_target_evidence_review.xlsx"
    ),
    overwrite = TRUE
  )
} else {
  cat("\nExcel export skipped; RDS and CSV outputs were saved.\n")
}

# ------------------------------------------------------------------------------
# 7. Completion report
# ------------------------------------------------------------------------------

cat("\n============================================\n")
cat("Target evidence review completed\n")
cat("============================================\n")

cat(
  "Main overlapping miRNA-target pairs:",
  nrow(main_summary),
  "\n"
)

cat(
  "Pairs with a reported reporter assay:",
  nrow(reporter_supported),
  "\n"
)

cat("\nEvidence groups:\n")
print(table(main_summary$Evidence_group))

cat("\nMain evidence summary:\n")
print(
  main_summary[, c(
    "miRNA",
    "Target_symbol",
    "Unique_PMIDs",
    "Reporter_assay_reported",
    "Evidence_group",
    "mRNA_logFC",
    "mRNA_FDR"
  )],
  row.names = FALSE
)

cat("\nProvisional miR-375 evidence summary:\n")
print(provisional_summary, row.names = FALSE)

cat("\nInterpretation notes:\n")
cat("Evidence records are not independent experimental replications.\n")
cat("Assay flags are based on database text and require paper-level review.\n")
cat("Reporter assays do not establish regulation in the present cohorts.\n")
cat("miR-375 alias mapping remains unverified.\n")
cat("No overlap-enrichment test was performed by this script.\n")
review_results <- readRDS(
  "results/R_objects/GSE274664_GSE287540_target_evidence_review.rds"
)

hgf_evidence <- review_results$Main_evidence_summary

print(
  hgf_evidence[
    hgf_evidence$Target_symbol == "HGF",
    c("miRNA", "Target_symbol", "Experiments",
      "Support_types", "PubMed_IDs")
  ],
  row.names = FALSE
)

cat("\nOutput directory:", csv_directory, "\n")