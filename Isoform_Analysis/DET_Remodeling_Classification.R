#!/usr/bin/env Rscript

# =============================================================================
# DET emergence and structural remodeling classification
# =============================================================================
#
# Purpose:
#   This script classifies differentially expressed transcripts (DETs) into:
#
#     1. Expressed at Baseline
#     2. Emerged: Activation
#     3. Emerged: Remodeling
#
#   For DETs classified as "Emerged: Remodeling", the script compares each DET
#   against the dominant baseline isoform of the same gene and assigns structural
#   remodeling categories:
#
#     - UTR Only
#     - ORF Modification
#     - ORF Gain
#     - ORF Loss
#     - Noncoding
#
# Inputs:
#   - Counts_filtered.tsv from DTE_Analysis_DESeq2.R
#   - C3_vs_C0_isoform_results.csv from DTE_Analysis_DESeq2.R
#   - C7_vs_C0_isoform_results.csv from DTE_Analysis_DESeq2.R
#   - Ensembl reference GTF
#   - renamed transcript GTF files from SQANTI3 to reflect novel transcripts
#   - SQANTI3 classification files with ORF prediction columns
#
# =============================================================================

suppressPackageStartupMessages({
  library(rtracklayer)
  library(tidyverse)
})

# -----------------------------
# User-defined files/directories
# -----------------------------

results_dir <- "/path/to/results_directory"                    # <-- MODIFY HERE

counts_filtered_path <- file.path(
  results_dir,
  "Counts_filtered.tsv"
)                                                              # <-- MODIFY HERE IF NEEDED

res_c3_path <- file.path(
  results_dir,
  "C3_vs_C0_isoform_results.csv"
)                                                              # <-- MODIFY HERE IF NEEDED

res_c7_path <- file.path(
  results_dir,
  "C7_vs_C0_isoform_results.csv"
)                                                              # <-- MODIFY HERE IF NEEDED

ref_gtf_path <- "/path/to/Mus_musculus.GRCm39.115.chr.gtf"      # <-- MODIFY HERE

renamed_gtf_dir <- "/path/to/renamed_or_unified_gtf_files"      # <-- MODIFY HERE

sqanti_class_dir <- "/path/to/SQANTI3_classification_files"     # <-- MODIFY HERE

outdir <- file.path(results_dir, "DET_Remodeling_Classification") # <-- MODIFY HERE IF NEEDED

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# User-adjustable parameters
# -----------------------------

min_count <- 10                                                 # <-- MODIFY HERE IF NEEDED
min_reps <- 2                                                   # <-- MODIFY HERE IF NEEDED
padj_cutoff <- 0.05                                             # <-- MODIFY HERE IF NEEDED
lfc_cutoff <- 1                                                 # <-- MODIFY HERE IF NEEDED

# Pattern for SQANTI3 classification files containing ORF prediction columns.
sqanti_class_pattern <- "classification.*\\.txt$"               # <-- MODIFY HERE IF NEEDED

# -----------------------------
# Input checks
# -----------------------------

required_files <- c(
  counts_filtered_path,
  res_c3_path,
  res_c7_path,
  ref_gtf_path
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required file(s):\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

if (!dir.exists(renamed_gtf_dir)) {
  stop("Missing renamed GTF directory: ", renamed_gtf_dir)
}

if (!dir.exists(sqanti_class_dir)) {
  stop("Missing SQANTI3 classification directory: ", sqanti_class_dir)
}

# =============================================================================
# Helper functions
# =============================================================================

get_baseline_status <- function(count_tbl, c0_cols, min_count = 10, min_reps = 2) {
  count_tbl %>%
    mutate(
      n_C0_ge_min = rowSums(across(all_of(c0_cols), ~ .x >= min_count)),
      iso_expressed_baseline = n_C0_ge_min >= min_reps
    ) %>%
    group_by(gene_symbol) %>%
    mutate(
      gene_has_any_baseline_isoform = any(iso_expressed_baseline, na.rm = TRUE),
      baseline_status = case_when(
        iso_expressed_baseline ~ "Expressed at Baseline",
        !iso_expressed_baseline & gene_has_any_baseline_isoform ~ "Emerged: Remodeling",
        !iso_expressed_baseline & !gene_has_any_baseline_isoform ~ "Emerged: Activation",
        TRUE ~ "Unknown"
      )
    ) %>%
    ungroup()
}

is_sig_det <- function(df) {
  df %>%
    filter(!is.na(padj)) %>%
    filter(padj <= padj_cutoff, abs(log2FoldChange) >= lfc_cutoff)
}

make_interval_key <- function(df, region) {
  df %>%
    filter(region_type == region) %>%
    arrange(seqnames, start, end, strand) %>%
    transmute(key = paste(seqnames, start, end, strand, sep = ":")) %>%
    pull(key) %>%
    paste(collapse = ";")
}

split_exon_by_orf <- function(exon_row) {
  exon_start <- exon_row$start
  exon_end <- exon_row$end
  cds_start <- exon_row$CDS_genomic_start
  cds_end <- exon_row$CDS_genomic_end

  if (is.na(cds_start) || is.na(cds_end)) {
    return(tibble(
      seqnames = exon_row$seqnames,
      start = exon_start,
      end = exon_end,
      strand = exon_row$strand,
      transcript_id = exon_row$transcript_id,
      region_type = "UTR"
    ))
  }

  cds_min <- min(cds_start, cds_end)
  cds_max <- max(cds_start, cds_end)

  pieces <- list()

  if (exon_start < cds_min) {
    pieces[[length(pieces) + 1]] <- tibble(
      seqnames = exon_row$seqnames,
      start = exon_start,
      end = min(exon_end, cds_min - 1),
      strand = exon_row$strand,
      transcript_id = exon_row$transcript_id,
      region_type = "UTR"
    )
  }

  if (exon_end >= cds_min && exon_start <= cds_max) {
    pieces[[length(pieces) + 1]] <- tibble(
      seqnames = exon_row$seqnames,
      start = max(exon_start, cds_min),
      end = min(exon_end, cds_max),
      strand = exon_row$strand,
      transcript_id = exon_row$transcript_id,
      region_type = "CDS"
    )
  }

  if (exon_end > cds_max) {
    pieces[[length(pieces) + 1]] <- tibble(
      seqnames = exon_row$seqnames,
      start = max(exon_start, cds_max + 1),
      end = exon_end,
      strand = exon_row$strand,
      transcript_id = exon_row$transcript_id,
      region_type = "UTR"
    )
  }

  bind_rows(pieces) %>%
    filter(start <= end)
}

compare_tx_to_ref <- function(tx_id, ref_tx, tx_structure_all, tx_coding_annot) {
  tx_structure <- tx_structure_all %>% filter(tx_id == !!tx_id)
  ref_structure <- tx_structure_all %>% filter(tx_id == !!ref_tx)

  tx_coding <- tx_coding_annot %>%
    filter(tx_id == !!tx_id) %>%
    pull(coding_class) %>%
    first()

  ref_coding <- tx_coding_annot %>%
    filter(tx_id == !!ref_tx) %>%
    pull(coding_class) %>%
    first()

  tx_cds_key <- make_interval_key(tx_structure, "CDS")
  ref_cds_key <- make_interval_key(ref_structure, "CDS")

  tx_utr_key <- make_interval_key(tx_structure, "UTR")
  ref_utr_key <- make_interval_key(ref_structure, "UTR")

  structure_change_class <- case_when(
    is.na(tx_coding) | is.na(ref_coding) ~ "Unknown",
    tx_coding != ref_coding ~ "Coding-Noncoding Switch",
    tx_coding == "coding" & ref_coding == "coding" & tx_cds_key != ref_cds_key ~ "CDS change",
    tx_cds_key == ref_cds_key & tx_utr_key != ref_utr_key ~ "UTR only",
    tx_coding == "noncoding" & ref_coding == "noncoding" ~ "Noncoding",
    TRUE ~ "No detected structural difference"
  )

  tibble(
    tx_id = tx_id,
    ref_tx = ref_tx,
    tx_coding_class = tx_coding,
    ref_coding_class = ref_coding,
    tx_cds_key = tx_cds_key,
    ref_cds_key = ref_cds_key,
    tx_utr_key = tx_utr_key,
    ref_utr_key = ref_utr_key,
    structure_change_class = structure_change_class
  )
}

# =============================================================================
# Step 1: Load DTE results and identify significant DETs
# =============================================================================

counts_filtered <- readr::read_tsv(counts_filtered_path, show_col_types = FALSE)

res_c3 <- readr::read_csv(res_c3_path, show_col_types = FALSE) %>%
  mutate(contrast = "C3_vs_C0", contrast_pretty = "C3 vs C0")

res_c7 <- readr::read_csv(res_c7_path, show_col_types = FALSE) %>%
  mutate(contrast = "C7_vs_C0", contrast_pretty = "C7 vs C0")

det_tbl <- bind_rows(res_c3, res_c7) %>%
  is_sig_det() %>%
  distinct(contrast, contrast_pretty, transcript_id, .keep_all = TRUE)

readr::write_csv(
  det_tbl,
  file.path(outdir, "significant_DETs_combined.csv")
)

# =============================================================================
# Step 2: Classify DET emergence relative to baseline
# =============================================================================

c0_cols <- names(counts_filtered)[str_detect(names(counts_filtered), "^C0")]

if (length(c0_cols) == 0) {
  stop("No C0 columns detected in counts table.")
}

baseline_tbl <- get_baseline_status(
  count_tbl = counts_filtered,
  c0_cols = c0_cols,
  min_count = min_count,
  min_reps = min_reps
) %>%
  select(
    transcript_id,
    transcript_symbol,
    gene_symbol,
    n_C0_ge_min,
    iso_expressed_baseline,
    gene_has_any_baseline_isoform,
    baseline_status
  ) %>%
  distinct(transcript_id, .keep_all = TRUE)

det_annot <- det_tbl %>%
  left_join(
    baseline_tbl,
    by = c("transcript_id", "transcript_symbol", "gene_symbol")
  ) %>%
  mutate(
    baseline_status = factor(
      baseline_status,
      levels = c(
        "Emerged: Remodeling",
        "Emerged: Activation",
        "Expressed at Baseline"
      )
    )
  )

readr::write_csv(
  det_annot,
  file.path(outdir, "DET_emergence_classification.csv")
)

emergence_summary <- det_annot %>%
  count(contrast_pretty, baseline_status) %>%
  group_by(contrast_pretty) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

readr::write_csv(
  emergence_summary,
  file.path(outdir, "DET_emergence_summary.csv")
)

det_emerge_remodel <- det_annot %>%
  filter(baseline_status == "Emerged: Remodeling") %>%
  select(transcript_id, transcript_symbol, gene_symbol, contrast, contrast_pretty)

det_emerge_activation <- det_annot %>%
  filter(baseline_status == "Emerged: Activation") %>%
  select(transcript_id, transcript_symbol, gene_symbol, contrast, contrast_pretty)

readr::write_csv(
  det_emerge_remodel,
  file.path(outdir, "DET_emerged_remodeling.csv")
)

readr::write_csv(
  det_emerge_activation,
  file.path(outdir, "DET_emerged_activation.csv")
)

# =============================================================================
# Step 3: Identify dominant baseline isoform per gene
# =============================================================================

count_long <- counts_filtered %>%
  select(transcript_id, transcript_symbol, gene_symbol, starts_with("C0_")) %>%
  pivot_longer(
    cols = starts_with("C0_"),
    names_to = "sample",
    values_to = "count"
  )

c0_dominant_isoform <- count_long %>%
  group_by(gene_symbol, transcript_symbol) %>%
  summarise(
    mean_count_C0 = mean(count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    baseline_tbl %>%
      select(transcript_symbol, gene_symbol, iso_expressed_baseline),
    by = c("transcript_symbol", "gene_symbol")
  ) %>%
  filter(iso_expressed_baseline) %>%
  group_by(gene_symbol) %>%
  slice_max(mean_count_C0, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(ref_tx = transcript_symbol)

# =============================================================================
# Step 4: Determine transcripts requiring structural comparison
# =============================================================================

isoform_structure_input <- det_emerge_remodel %>%
  mutate(tx_id = transcript_symbol) %>%
  left_join(c0_dominant_isoform, by = "gene_symbol") %>%
  filter(!is.na(ref_tx), tx_id != ref_tx)

tx_keep <- unique(c(
  isoform_structure_input$tx_id,
  isoform_structure_input$ref_tx
)) %>%
  na.omit()

known_tx_keep <- tx_keep[!str_detect(tx_keep, "novel")]
novel_tx_keep <- tx_keep[str_detect(tx_keep, "novel")]

# =============================================================================
# Step 5: Get known transcript structure from Ensembl reference GTF
# =============================================================================

message("Importing reference GTF: ", ref_gtf_path)

ref_gtf <- rtracklayer::import(ref_gtf_path) %>%
  as.data.frame()

known_structure <- ref_gtf %>%
  mutate(
    transcript_name = if_else(
      is.na(transcript_name),
      transcript_id,
      transcript_name
    )
  ) %>%
  filter(
    transcript_name %in% known_tx_keep,
    type %in% c("CDS", "five_prime_utr", "three_prime_utr")
  ) %>%
  transmute(
    tx_id = transcript_name,
    seqnames,
    start,
    end,
    strand,
    region_type = if_else(type == "CDS", "CDS", "UTR"),
    source = "Ensembl_reference"
  ) %>%
  distinct()

known_coding_by_cds <- ref_gtf %>%
  mutate(
    transcript_name = if_else(
      is.na(transcript_name),
      transcript_id,
      transcript_name
    )
  ) %>%
  filter(transcript_name %in% known_tx_keep) %>%
  group_by(transcript_name) %>%
  summarise(
    has_CDS = any(type == "CDS"),
    transcript_biotype = first(na.omit(transcript_biotype)),
    .groups = "drop"
  ) %>%
  mutate(
    tx_id = transcript_name,
    coding_class = if_else(has_CDS, "coding", "noncoding"),
    coding_source = "Ensembl_CDS_presence"
  ) %>%
  select(tx_id, coding_class, coding_source, transcript_biotype)

# =============================================================================
# Step 6: Get novel transcript structure and ORF predictions
# =============================================================================

class_files <- list.files(
  sqanti_class_dir,
  pattern = sqanti_class_pattern,
  full.names = TRUE
)

if (length(class_files) == 0) {
  stop("No SQANTI3 classification files found in: ", sqanti_class_dir)
}

orf_sqanti <- map_dfr(class_files, ~ readr::read_tsv(.x, show_col_types = FALSE)) %>%
  transmute(
    transcript_id,
    strand,
    coding,
    CDS_genomic_start_raw = suppressWarnings(as.numeric(CDS_genomic_start)),
    CDS_genomic_end_raw = suppressWarnings(as.numeric(CDS_genomic_end)),
    cds_len = abs(CDS_genomic_end_raw - CDS_genomic_start_raw)
  ) %>%
  group_by(transcript_id) %>%
  arrange(desc(replace_na(cds_len, -Inf)), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    transcript_id,
    strand,
    coding,
    CDS_genomic_start = CDS_genomic_start_raw,
    CDS_genomic_end = CDS_genomic_end_raw
  ) %>%
  distinct()

gtf_files <- list.files(
  renamed_gtf_dir,
  pattern = "\\.gtf$",
  full.names = TRUE
)

if (length(gtf_files) == 0) {
  stop("No GTF files found in: ", renamed_gtf_dir)
}

gtf_all <- map_dfr(gtf_files, function(f) {
  rtracklayer::import(f) %>%
    as.data.frame() %>%
    filter(transcript_id %in% tx_keep) %>%
    mutate(source_gtf = basename(f))
})

cds_from_gtf <- gtf_all %>%
  filter(type == "CDS") %>%
  group_by(transcript_id) %>%
  summarise(
    CDS_genomic_start = min(start),
    CDS_genomic_end = max(end),
    .groups = "drop"
  )

orf_combined <- full_join(
  orf_sqanti,
  cds_from_gtf,
  by = "transcript_id"
) %>%
  transmute(
    transcript_id,
    coding,
    CDS_genomic_start = coalesce(CDS_genomic_start.x, CDS_genomic_start.y),
    CDS_genomic_end = coalesce(CDS_genomic_end.x, CDS_genomic_end.y)
  )

gtf_exons <- gtf_all %>%
  filter(type == "exon") %>%
  left_join(orf_combined, by = "transcript_id")

novel_selected <- gtf_exons %>%
  filter(transcript_id %in% novel_tx_keep) %>%
  arrange(transcript_id, start) %>%
  distinct(
    seqnames,
    start,
    end,
    strand,
    transcript_id,
    gene_id,
    coding,
    CDS_genomic_start,
    CDS_genomic_end,
    .keep_all = TRUE
  )

novel_structure <- novel_selected %>%
  rowwise() %>%
  do(split_exon_by_orf(.)) %>%
  ungroup() %>%
  mutate(
    tx_id = transcript_id,
    region_type = case_when(
      region_type == "CDS" ~ "CDS",
      region_type %in% c("5UTR", "3UTR", "UTR") ~ "UTR",
      TRUE ~ NA_character_
    ),
    source = "SQANTI3_ORF_prediction"
  ) %>%
  select(tx_id, seqnames, start, end, strand, region_type, source) %>%
  filter(!is.na(region_type)) %>%
  distinct()

novel_coding <- orf_combined %>%
  filter(transcript_id %in% novel_tx_keep) %>%
  mutate(
    tx_id = transcript_id,
    has_predicted_orf = !is.na(CDS_genomic_start) & !is.na(CDS_genomic_end),
    coding_class = if_else(coding == "coding", "coding", "noncoding"),
    coding_source = "SQANTI3_ORF_prediction",
    transcript_biotype = if_else(
      has_predicted_orf,
      "predicted_protein_coding",
      "predicted_noncoding"
    )
  ) %>%
  select(tx_id, coding_class, coding_source, transcript_biotype)

# =============================================================================
# Step 7: Compare DETs to dominant baseline isoforms
# =============================================================================

tx_structure_all <- bind_rows(
  known_structure,
  novel_structure
) %>%
  filter(!is.na(region_type)) %>%
  distinct()

tx_coding_annot <- bind_rows(
  known_coding_by_cds,
  novel_coding
) %>%
  distinct(tx_id, .keep_all = TRUE)

compare_df <- isoform_structure_input %>%
  distinct(tx_id, ref_tx, gene_symbol, contrast, contrast_pretty) %>%
  mutate(
    compare = map2(
      tx_id,
      ref_tx,
      ~ compare_tx_to_ref(.x, .y, tx_structure_all, tx_coding_annot)
    )
  ) %>%
  select(-tx_id, -ref_tx) %>%
  unnest(compare)

isoform_structure_gr <- isoform_structure_input %>%
  left_join(
    compare_df,
    by = c("tx_id", "ref_tx", "gene_symbol", "contrast", "contrast_pretty")
  ) %>%
  mutate(
    coding_transition = case_when(
      ref_coding_class == "coding" & tx_coding_class == "coding" ~ "coding_to_coding",
      ref_coding_class == "coding" & tx_coding_class == "noncoding" ~ "coding_to_noncoding",
      ref_coding_class == "noncoding" & tx_coding_class == "coding" ~ "noncoding_to_coding",
      ref_coding_class == "noncoding" & tx_coding_class == "noncoding" ~ "noncoding_to_noncoding",
      TRUE ~ "unknown"
    ),
    comparison_type = case_when(
      !str_detect(tx_id, "novel") & !str_detect(ref_tx, "novel") ~ "known_vs_known",
      str_detect(tx_id, "novel") & !str_detect(ref_tx, "novel") ~ "novel_vs_known",
      !str_detect(tx_id, "novel") & str_detect(ref_tx, "novel") ~ "known_vs_novel",
      str_detect(tx_id, "novel") & str_detect(ref_tx, "novel") ~ "novel_vs_novel",
      TRUE ~ "unknown"
    ),
    final_class = case_when(
      structure_change_class == "UTR only" ~ "UTR Only",
      structure_change_class == "CDS change" ~ "ORF Modification",
      structure_change_class == "Coding-Noncoding Switch" & coding_transition == "noncoding_to_coding" ~ "ORF Gain",
      structure_change_class == "Coding-Noncoding Switch" & coding_transition == "coding_to_noncoding" ~ "ORF Loss",
      structure_change_class == "Coding-Noncoding Switch" & coding_transition == "noncoding_to_noncoding" ~ "Noncoding",
      structure_change_class == "Noncoding" ~ "Noncoding",
      TRUE ~ "unknown"
    )
  )

readr::write_csv(
  isoform_structure_gr,
  file.path(outdir, "DET_structural_remodeling_classification.csv")
)

structure_summary <- isoform_structure_gr %>%
  count(contrast_pretty, final_class) %>%
  group_by(contrast_pretty) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

readr::write_csv(
  structure_summary,
  file.path(outdir, "DET_structural_remodeling_summary.csv")
)

message("DET remodeling classification complete. Results saved to: ", outdir)
