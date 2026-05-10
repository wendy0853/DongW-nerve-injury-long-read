#!/usr/bin/env Rscript

# =============================================================================
# Differential Transcript Usage Analysis with IsoformSwitchAnalyzeR + DEXSeq
# =============================================================================
#
# Purpose:
#   This script performs differential transcript usage (DTU) and isoform
#   switching analysis using IsoformSwitchAnalyzeR and DEXSeq.
#
# Inputs:
#   - Counts_filtered.tsv
#            Filtered for at least 10 counts for 2 biological samples from 
#            Counts Matrix from Generating_Count_Matrix.R
#   - all_transcripts_with_associated_genes_tpm.tsv 
#            TPM Matrix from Generating_Count_Matrix.R
#   - novel_isoform_translation.csv
#   - SQANTI3 *_corrected.gtf files
#
# Outputs:
#   - ALL_isoform_ratio_per_gene_condition.csv
#   - combined_longread_raw.gtf
#   - filtered_longread.gtf
#   - switch_list_C3C7C0.rds
#   - All_Isoform_Switching_Events_Annotated.csv
#   - All_Gene_Level_Switching_Events.csv
#
# Notes:
#   The GTF inputs should be SQANTI3 *_corrected.gtf files.
#
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(biomaRt)
  library(rtracklayer)
  library(IsoformSwitchAnalyzeR)
})

# -----------------------------
# User-defined files/directories
# -----------------------------

results_dir <- "/path/to/results_directory"                         # <-- MODIFY HERE
sqanti_gtf_dir <- "/path/to/SQANTI3_corrected_gtf_files"             # <-- MODIFY HERE

counts_path <- file.path(
  results_dir,
  "Counts_filtered.tsv"
)                                                                    # <-- MODIFY HERE IF NEEDED

tpm_path <- file.path(
  results_dir,
  "all_transcripts_with_associated_genes_tpm.tsv"
)                                                                    # <-- MODIFY HERE IF NEEDED

novel_lookup_path <- file.path(
  results_dir,
  "novel_isoform_translation.csv"
)                                                                    # <-- MODIFY HERE IF NEEDED

outdir <- "/path/to/isoform_switching_results"                       # <-- MODIFY HERE

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# User-adjustable parameters
# -----------------------------

IFcutoff <- 0.05                                                     # <-- MODIFY HERE IF NEEDED
dIF_cutoff <- 0.1                                                    # <-- MODIFY HERE IF NEEDED
padj_cutoff <- 0.05                                                  # <-- MODIFY HERE IF NEEDED

sample_names <- c(                                                   # <-- MODIFY HERE IF NEEDED
  "C0_Control_1", "C0_Control_2", "C0_Control_3",
  "C3_Injured_1", "C3_Injured_2", "C3_Injured_3",
  "C7_Injured_1", "C7_Injured_2", "C7_Injured_3"
)

condition_vector <- rep(c("C0", "C3", "C7"), each = 3)              # <-- MODIFY HERE IF NEEDED

# -----------------------------
# Input checks
# -----------------------------

required_files <- c(
  counts_path,
  tpm_path,
  novel_lookup_path
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required file(s):\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

if (!dir.exists(sqanti_gtf_dir)) {
  stop("Missing SQANTI3 corrected GTF directory: ", sqanti_gtf_dir)
}

# =============================================================================
# Step 1: Load filtered count and TPM matrices
# =============================================================================

counts_filtered <- readr::read_tsv(counts_path, show_col_types = FALSE)

filtered_ids <- counts_filtered %>%
  select(transcript_id)

tpm_filtered <- readr::read_tsv(tpm_path, show_col_types = FALSE) %>%
  filter(transcript_id %in% filtered_ids$transcript_id)

if ("gene_symbol" %in% colnames(tpm_filtered)) {
  tpm_filtered <- tpm_filtered %>%
    rename(gene_id = gene_symbol)
}

tpm_filtered <- tpm_filtered %>%
  select(
    transcript_id,
    gene_id,
    everything()
  )

# =============================================================================
# Step 2: Calculate isoform fractions per gene and sample
# =============================================================================

expr_long <- tpm_filtered %>%
  pivot_longer(
    cols = all_of(sample_names),
    names_to = "sample",
    values_to = "tpm"
  ) %>%
  group_by(sample, gene_id) %>%
  mutate(
    total_gene_tpm = sum(tpm, na.rm = TRUE),
    isoform_fraction = if_else(total_gene_tpm > 0, tpm / total_gene_tpm, 0)
  ) %>%
  ungroup() %>%
  mutate(
    condition = case_when(
      str_detect(sample, "^C0") ~ "C0",
      str_detect(sample, "^C3") ~ "C3",
      str_detect(sample, "^C7") ~ "C7",
      TRUE ~ "Other"
    )
  )

readr::write_csv(
  expr_long,
  file.path(outdir, "ALL_isoform_ratio_per_gene_condition.csv")
)

# =============================================================================
# Step 3: Combine SQANTI3 corrected GTF files
# =============================================================================
#
# IMPORTANT:
#   These should be SQANTI3-generated *_corrected.gtf files, not raw IsoQuant GTFs.

gtf_files <- list.files(
  sqanti_gtf_dir,
  pattern = "_corrected\\.gtf$",
  full.names = TRUE
)

if (length(gtf_files) == 0) {
  stop("No SQANTI3 *_corrected.gtf files found in: ", sqanti_gtf_dir)
}

message("Number of SQANTI3 corrected GTF files: ", length(gtf_files))

gtf_list <- lapply(gtf_files, rtracklayer::import)
gtf_combined <- do.call(c, gtf_list)

gtf_combined$gene_id <- stringr::str_remove(gtf_combined$gene_id, "\\.\\d+$")
gtf_combined$transcript_id <- stringr::str_remove(gtf_combined$transcript_id, "\\.\\d+$")

combined_gtf_path <- file.path(outdir, "combined_longread_raw.gtf")

rtracklayer::export(
  gtf_combined,
  combined_gtf_path
)

message("Wrote combined GTF: ", combined_gtf_path)

# =============================================================================
# Step 4: Prepare count and TPM matrices for IsoformSwitchAnalyzeR
# =============================================================================

tpm_matrix <- tpm_filtered %>%
  select(transcript_id, all_of(sample_names)) %>%
  rename(isoform_id = transcript_id)

count_matrix <- counts_filtered %>%
  select(transcript_id, all_of(sample_names)) %>%
  rename(isoform_id = transcript_id)

# IsoformSwitchAnalyzeR/DEXSeq generally handles cleaner IDs better.
# Remove Ensembl version suffixes for matching against cleaned GTF IDs.
tpm_matrix <- tpm_matrix %>%
  mutate(isoform_id = stringr::str_remove(isoform_id, "\\.\\d+$"))

count_matrix <- count_matrix %>%
  mutate(isoform_id = stringr::str_remove(isoform_id, "\\.\\d+$"))

# =============================================================================
# Step 5: Map unified novel IDs back to SQANTI3 transcript IDs for GTF matching
# =============================================================================

novel_to_id_map <- readr::read_csv(novel_lookup_path, show_col_types = FALSE) %>%
  distinct(unified_id, .keep_all = TRUE) %>%
  mutate(
    transcript_id = stringr::str_remove(transcript_id, "_sample\\d+$"),
    transcript_id = stringr::str_remove(transcript_id, "_\\d+$")
  )

tpm_matrix <- tpm_matrix %>%
  mutate(
    isoform_id = if_else(
      isoform_id %in% novel_to_id_map$unified_id,
      novel_to_id_map$transcript_id[match(isoform_id, novel_to_id_map$unified_id)],
      isoform_id
    )
  )

count_matrix <- count_matrix %>%
  mutate(
    isoform_id = if_else(
      isoform_id %in% novel_to_id_map$unified_id,
      novel_to_id_map$transcript_id[match(isoform_id, novel_to_id_map$unified_id)],
      isoform_id
    )
  )

tpm_matrix <- tpm_matrix[!duplicated(tpm_matrix$isoform_id), ]
count_matrix <- count_matrix[!duplicated(count_matrix$isoform_id), ]

# =============================================================================
# Step 6: Filter GTF to isoforms present in count matrix
# =============================================================================

gtf <- rtracklayer::import(combined_gtf_path)

gtf$transcript_id <- stringr::str_remove(gtf$transcript_id, "\\.\\d+$")
gtf$gene_id <- stringr::str_remove(gtf$gene_id, "\\.\\d+$")
mcols(gtf)$isoform_id <- mcols(gtf)$transcript_id

count_ids <- stringr::str_remove(count_matrix$isoform_id, "\\.\\d+$")

gtf_filtered <- gtf[gtf$transcript_id %in% count_ids]

filtered_gtf_path <- file.path(outdir, "filtered_longread.gtf")

rtracklayer::export(
  gtf_filtered,
  filtered_gtf_path
)

message("Wrote filtered GTF: ", filtered_gtf_path)

# =============================================================================
# Step 7: Create design matrix
# =============================================================================

design_matrix <- data.frame(
  sampleID = sample_names,
  condition = condition_vector
)

if (!all(sample_names %in% colnames(count_matrix))) {
  stop("Not all sample_names are present in count_matrix.")
}

counts_mat <- count_matrix %>%
  tibble::column_to_rownames("isoform_id") %>%
  as.matrix()

tpm_mat <- tpm_matrix %>%
  tibble::column_to_rownames("isoform_id") %>%
  as.matrix()

# =============================================================================
# Step 8: Import data into IsoformSwitchAnalyzeR
# =============================================================================

switch_list <- importRdata(
  isoformCountMatrix = counts_mat,
  designMatrix = design_matrix,
  isoformExonAnnoation = filtered_gtf_path,
  showProgress = TRUE
)

summary(switch_list)

# =============================================================================
# Step 9: Prefilter and run DEXSeq-based DTU testing
# =============================================================================

switch_list <- preFilter(
  switch_list,
  IFcutoff = IFcutoff,
  removeSingleIsoformGenes = FALSE
)

switch_list <- isoformSwitchTestDEXSeq(switch_list)

# =============================================================================
# Step 10: Annotate isoform features with gene symbols
# =============================================================================

message("Querying Ensembl BioMart for transcript annotations...")

ensembl <- biomaRt::useEnsembl(
  biomart = "ensembl",
  dataset = "mmusculus_gene_ensembl"
)

transcript_ids <- unique(switch_list$isoformFeatures$isoform_id)

transcript_lookup <- biomaRt::getBM(
  attributes = c(
    "ensembl_transcript_id",
    "ensembl_gene_id",
    "external_gene_name"
  ),
  filters = "ensembl_transcript_id",
  values = transcript_ids,
  mart = ensembl
) %>%
  mutate(
    ensembl_transcript_id = as.character(ensembl_transcript_id),
    ensembl_gene_id = as.character(ensembl_gene_id),
    external_gene_name = as.character(external_gene_name)
  )

switch_list$isoformFeatures <- switch_list$isoformFeatures %>%
  mutate(isoform_id = as.character(isoform_id)) %>%
  left_join(
    transcript_lookup,
    by = c("isoform_id" = "ensembl_transcript_id")
  ) %>%
  mutate(
    gene_id = coalesce(ensembl_gene_id, gene_id),
    gene_name = external_gene_name
  ) %>%
  select(-ensembl_gene_id, -external_gene_name)

# =============================================================================
# Step 11: Save switch object
# =============================================================================

switch_list_path <- file.path(outdir, "switch_list_C3C7C0.rds")

saveRDS(
  switch_list,
  switch_list_path
)

message("Saved IsoformSwitchAnalyzeR object: ", switch_list_path)

# =============================================================================
# Step 12: Export isoform-level DTU results
# =============================================================================

all_switches <- switch_list$isoformSwitchAnalysis

annotation_from_tpm <- tpm_filtered %>%
  select(
    transcript_id,
    transcript_symbol,
    associated_gene,
    gene_id
  ) %>%
  distinct() %>%
  mutate(isoform_id = transcript_id)

all_switches_annotated <- all_switches %>%
  left_join(annotation_from_tpm, by = "isoform_id") %>%
  mutate(
    transcript_symbol = if_else(
      isoform_id %in% novel_to_id_map$transcript_id,
      novel_to_id_map$unified_id[match(isoform_id, novel_to_id_map$transcript_id)],
      transcript_symbol
    ),
    gene_id = if_else(
      is.na(gene_id) | gene_id == "",
      stringr::str_remove(transcript_symbol, "-novel-\\d+.*$"),
      gene_id
    ),
    significant_DTU = !is.na(padj) & padj <= padj_cutoff & abs(dIF) >= dIF_cutoff
  ) %>%
  arrange(padj)

isoform_switch_path <- file.path(
  outdir,
  "All_Isoform_Switching_Events_Annotated.csv"
)

write.csv(
  all_switches_annotated,
  isoform_switch_path,
  row.names = FALSE
)

message("Wrote isoform-level switching results: ", isoform_switch_path)

# =============================================================================
# Step 13: Summarize gene-level switching events
# =============================================================================

all_gene_switches <- all_switches_annotated %>%
  mutate(
    padj = as.numeric(padj),
    dIF = as.numeric(dIF),
    gene_id = if_else(
      is.na(gene_id) | gene_id == "",
      stringr::str_remove(transcript_symbol, "-novel-\\d+.*$"),
      gene_id
    ),
    gene_name = gene_id
  ) %>%
  group_by(gene_ref, gene_name, condition_1, condition_2) %>%
  summarise(
    gene_switch_q_value = suppressWarnings(min(padj, na.rm = TRUE)),
    total_usage_change = 0.5 * sum(abs(dIF), na.rm = TRUE),
    max_abs_dIF = suppressWarnings(max(abs(dIF), na.rm = TRUE)),
    n_isoforms = n(),
    n_significant_isoforms = sum(significant_DTU, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(gene_switch_q_value) %>%
  mutate(Rank = row_number())

gene_switch_path <- file.path(
  outdir,
  "All_Gene_Level_Switching_Events.csv"
)

write.csv(
  all_gene_switches,
  gene_switch_path,
  row.names = FALSE
)

message("Wrote gene-level switching results: ", gene_switch_path)

# =============================================================================
# Step 14: Export significant DTU subset
# =============================================================================

significant_dtu <- all_switches_annotated %>%
  filter(significant_DTU)

write.csv(
  significant_dtu,
  file.path(outdir, "Significant_Isoform_Switching_Events.csv"),
  row.names = FALSE
)

message("DTU analysis complete. Results saved to: ", outdir)
