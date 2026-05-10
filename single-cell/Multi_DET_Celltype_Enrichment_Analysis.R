#!/usr/bin/env Rscript

# =============================================================================
# Cell-type enrichment of multi-DET genes in pseudobulk DEGs
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# -----------------------------
# User-defined files/directories
# -----------------------------

single_cell_results_dir <- "/path/to/single_cell_results"        # <-- MODIFY HERE
isoform_results_dir <- "/path/to/isoform_analysis_results"       # <-- MODIFY HERE
outdir <- single_cell_results_dir                                # <-- MODIFY HERE IF NEEDED

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# User-defined input files
# -----------------------------

counts_filtered_path <- file.path(
  isoform_results_dir,
  "Counts_filtered.tsv"
)                                                                 # <-- MODIFY HERE IF NEEDED

c3_isoform_results_path <- file.path(
  isoform_results_dir,
  "C3_vs_C0_isoform_results.csv"
)                                                                 # <-- MODIFY HERE IF NEEDED

c7_isoform_results_path <- file.path(
  isoform_results_dir,
  "C7_vs_C0_isoform_results.csv"
)                                                                 # <-- MODIFY HERE IF NEEDED

# -----------------------------
# User-adjustable thresholds
# -----------------------------

padj_cutoff <- 0.05                                               # <-- MODIFY HERE IF NEEDED
lfc_cutoff <- 1                                                   # <-- MODIFY HERE IF NEEDED

# -----------------------------
# Helper functions
# -----------------------------

get_multi_det_genes <- function(det_path, padj_cutoff = 0.05, lfc_cutoff = 1) {
  readr::read_csv(det_path, show_col_types = FALSE) %>%
    filter(!is.na(padj)) %>%
    filter(padj <= padj_cutoff, abs(log2FoldChange) >= lfc_cutoff) %>%
    filter(!is.na(gene_symbol), gene_symbol != "") %>%
    group_by(gene_symbol) %>%
    summarise(n_DET_isoforms = n_distinct(transcript_id), .groups = "drop") %>%
    filter(n_DET_isoforms >= 2)
}

read_pseudobulk_deg <- function(path, padj_cutoff = 0.05, lfc_cutoff = 1) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(
      gene_symbol = case_when(
        "gene_symbol" %in% colnames(.) ~ gene_symbol,
        "gene" %in% colnames(.) ~ gene,
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(gene_symbol), gene_symbol != "") %>%
    filter(!is.na(padj)) %>%
    filter(padj <= padj_cutoff, abs(log2FoldChange) >= lfc_cutoff) %>%
    pull(gene_symbol) %>%
    unique()
}

run_fisher_enrichment <- function(query_genes, background_genes, target_genes) {
  query_genes <- intersect(unique(query_genes), background_genes)
  target_genes <- intersect(unique(target_genes), background_genes)

  a <- length(intersect(query_genes, target_genes))
  b <- length(setdiff(query_genes, target_genes))
  c <- length(setdiff(target_genes, query_genes))
  d <- length(setdiff(background_genes, union(query_genes, target_genes)))

  mat <- matrix(c(a, b, c, d), nrow = 2)

  test <- fisher.test(mat, alternative = "greater")

  tibble(
    overlap = a,
    query_genes = length(query_genes),
    target_genes = length(target_genes),
    background_genes = length(background_genes),
    odds_ratio = unname(test$estimate),
    p_value = test$p.value
  )
}

# -----------------------------
# Load background universe
# -----------------------------

counts_filtered <- readr::read_tsv(counts_filtered_path, show_col_types = FALSE)

background_genes <- counts_filtered %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  pull(gene_symbol) %>%
  unique()

# -----------------------------
# Get multi-DET genes
# -----------------------------

multi_det_c3 <- get_multi_det_genes(c3_isoform_results_path, padj_cutoff, lfc_cutoff)
multi_det_c7 <- get_multi_det_genes(c7_isoform_results_path, padj_cutoff, lfc_cutoff)

readr::write_csv(
  multi_det_c3,
  file.path(outdir, "C3_multi_DET_genes.csv")
)

readr::write_csv(
  multi_det_c7,
  file.path(outdir, "C7_multi_DET_genes.csv")
)

multi_det_list <- list(
  C3_vs_C0 = multi_det_c3$gene_symbol,
  C7_vs_C0 = multi_det_c7$gene_symbol
)

# -----------------------------
# Locate pseudobulk DESeq2 result files
# -----------------------------

pseudobulk_files <- list.files(
  single_cell_results_dir,
  pattern = "^DE_.*_C[37]_vs_C0\\.csv$",
  full.names = TRUE
)

if (length(pseudobulk_files) == 0) {
  stop("No pseudobulk DESeq2 files found in: ", single_cell_results_dir)
}

# -----------------------------
# Run enrichment
# -----------------------------

enrichment_results <- purrr::map_dfr(pseudobulk_files, function(file) {
  file_name <- basename(file)

  contrast <- case_when(
    str_detect(file_name, "C3_vs_C0") ~ "C3_vs_C0",
    str_detect(file_name, "C7_vs_C0") ~ "C7_vs_C0",
    TRUE ~ NA_character_
  )

  cell_type <- file_name %>%
    str_remove("^DE_") %>%
    str_remove("_C3_vs_C0\\.csv$") %>%
    str_remove("_C7_vs_C0\\.csv$") %>%
    str_replace_all("_", " ")

  deg_genes <- read_pseudobulk_deg(file, padj_cutoff, lfc_cutoff)

  run_fisher_enrichment(
    query_genes = multi_det_list[[contrast]],
    background_genes = background_genes,
    target_genes = deg_genes
  ) %>%
    mutate(
      cell_type = cell_type,
      contrast = contrast,
      source_file = file_name
    )
}) %>%
  group_by(contrast) %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH"),
    neg_log10_FDR = -log10(FDR)
  ) %>%
  ungroup() %>%
  select(
    contrast,
    cell_type,
    overlap,
    query_genes,
    target_genes,
    background_genes,
    odds_ratio,
    p_value,
    FDR,
    neg_log10_FDR,
    source_file
  )

readr::write_csv(
  enrichment_results,
  file.path(outdir, "multi_DET_cell_type_enrichment_results.csv")
)

message("Cell-type enrichment analysis complete.")
