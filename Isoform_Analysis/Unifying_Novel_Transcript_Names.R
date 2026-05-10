#!/usr/bin/env Rscript

# =============================================================================
# Unify SQANTI3 novel transcript names across samples
# =============================================================================
#
# Purpose:
#   SQANTI3 assigns novel transcript IDs independently for each sample.
#   This script identifies structurally identical novel transcripts across
#   sample-specific GTF files and assigns stable unified names:
#
#     <gene_symbol>-novel-<index>
#
# Output:
#   novel_isoform_translation.csv
#
# IMPORTANT:
#   The input GTFs should be the SQANTI3-generated *_corrected.gtf files
#   produced by sqanti3_qc.py, NOT the original IsoQuant GTFs.
#
# =============================================================================

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(stringr)
  library(biomaRt)
  library(readr)
  library(tidyr)
})

# -----------------------------
# User-defined files/directories
# -----------------------------

gtf_list_file <- "/path/to/sqanti3_corrected_gtf_list.txt"   # <-- MODIFY HERE
outdir <- "/path/to/output_directory"                        # <-- MODIFY HERE
output_file <- "novel_isoform_translation.csv"               # <-- MODIFY HERE IF NEEDED

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Read input GTF list
# -----------------------------
# gtf_list_file should contain one SQANTI3 *_corrected.gtf file per line.
#
# Example:
# /path/to/C0_Sciatic_1_corrected.gtf
# /path/to/C0_Sciatic_2_corrected.gtf
# /path/to/C0_Sciatic_3_corrected.gtf

gtf_paths <- readr::read_lines(gtf_list_file)
gtf_paths <- gtf_paths[gtf_paths != ""]

if (length(gtf_paths) == 0) {
  stop("No GTF paths found in gtf_list_file")
}

missing_files <- gtf_paths[!file.exists(gtf_paths)]

if (length(missing_files) > 0) {
  stop(
    "The following GTF files do not exist:\n",
    paste(missing_files, collapse = "\n")
  )
}

message("Number of SQANTI3 corrected GTF files: ", length(gtf_paths))

# -----------------------------
# Function to read novel transcripts
# -----------------------------

read_novel_transcripts <- function(gtf_path, sample_index) {
  message("Importing: ", gtf_path)

  gr <- rtracklayer::import(gtf_path, format = "gtf")

  exons <- gr[
    gr$type == "exon" &
      stringr::str_starts(gr$transcript_id, "transcript")
  ]

  if (length(exons) == 0) {
    warning("No SQANTI3 novel transcript IDs found in: ", gtf_path)
    return(GRangesList())
  }

  tx_grl <- split(exons, exons$transcript_id)

  # Append sample index so sample-specific SQANTI3 transcript IDs remain unique
  names(tx_grl) <- paste0(names(tx_grl), "_sample", sample_index)

  GRangesList(lapply(tx_grl, function(x) {
    mcols(x)$gene_id <- x$gene_id[1]
    x
  }))
}

# -----------------------------
# Import and combine transcripts
# -----------------------------

all_tx_grl <- purrr::imap(gtf_paths, read_novel_transcripts) |>
  do.call(what = c)

message("Total novel transcript instances collected: ", length(all_tx_grl))

if (length(all_tx_grl) == 0) {
  stop("No novel transcript models were detected.")
}

# -----------------------------
# Build transcript fingerprints
# -----------------------------

tx_meta <- tibble(
  transcript_id = names(all_tx_grl),
  gr = as.list(all_tx_grl)
) %>%
  mutate(
    gene_id = map_chr(gr, ~ as.character(mcols(.x)$gene_id[1])),
    gene_id_core = sub("\\.\\d+$", "", gene_id),
    fingerprint = map_chr(gr, ~ paste0(
      as.character(seqnames(.x)),
      ":",
      start(.x),
      "-",
      end(.x),
      "(",
      as.character(strand(.x)),
      ")",
      collapse = ";"
    ))
  ) %>%
  select(transcript_id, gene_id, gene_id_core, fingerprint)

# -----------------------------
# Map Ensembl gene IDs to gene symbols
# -----------------------------

message("Querying Ensembl BioMart for gene symbols...")

ensembl <- biomaRt::useMart(
  biomart = "ensembl",
  dataset = "mmusculus_gene_ensembl"
)

gene_map <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = unique(tx_meta$gene_id_core),
  mart = ensembl
) %>%
  rename(
    gene_id_core = ensembl_gene_id,
    gene_symbol = external_gene_name
  ) %>%
  right_join(
    tibble(gene_id_core = unique(tx_meta$gene_id_core)),
    by = "gene_id_core"
  ) %>%
  mutate(
    gene_symbol = ifelse(
      is.na(gene_symbol) | gene_symbol == "",
      gene_id_core,
      gene_symbol
    )
  )

tx_meta <- tx_meta %>%
  left_join(gene_map, by = "gene_id_core") %>%
  mutate(
    gene_symbol = coalesce(gene_symbol, gene_id_core)
  )

# -----------------------------
# Assign unified novel transcript IDs
# -----------------------------

lookup <- tx_meta %>%
  distinct(gene_id, gene_symbol, fingerprint) %>%
  group_by(gene_id, gene_symbol) %>%
  arrange(fingerprint, .by_group = TRUE) %>%
  mutate(
    novel_idx = row_number(),
    unified_id = paste0(gene_symbol, "-novel-", novel_idx)
  ) %>%
  ungroup() %>%
  inner_join(
    tx_meta %>% select(transcript_id, gene_id, fingerprint),
    by = c("gene_id", "fingerprint")
  ) %>%
  select(unified_id, transcript_id) %>%
  arrange(unified_id, transcript_id)

# -----------------------------
# Write output
# -----------------------------

out_path <- file.path(outdir, output_file)
readr::write_csv(lookup, out_path)

message("Wrote novel isoform translation table: ", out_path)
