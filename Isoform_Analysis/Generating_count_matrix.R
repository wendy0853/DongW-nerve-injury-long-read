
# ---------------------------------------------------------------
# SQANTI3 Long-Read RNA-seq Isoform Processing and Annotation
#
# Author: Wendy Dong
# Affiliation: Washington University School of Medicine
#
# Description:
#   This script filters SQANTI3 isoforms, rescues isoforms supported by
#   automatic inclusion lists, unifies sample-specific novel transcript IDs,
#   merges per-sample abundance estimates into a count matrix, and annotates
#   transcripts/genes using SQANTI3 classification output and BioMart.
#
# Important transcript-ID handling:
#   - Transcript IDs are preserved with their Ensembl version suffixes
#     throughout the count table whenever they are present.
#   - Version suffixes are stripped ONLY in temporary columns used for
#     BioMart lookup, because BioMart expects unversioned Ensembl IDs.
# ---------------------------------------------------------------


# --- Libraries --------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(biomaRt)
})


# --- User parameters --------------------------------------------
# Edit these paths before running.

project_dir <- "/Users/wendydong/Documents/WUSM PhD/Long Read RNA Seq/Analysis_Long_Short_updated"
data_dir    <- "/Users/wendydong/Documents/WUSM PhD/Long Read RNA Seq/Long_read_results_updated"
results_dir <- file.path(project_dir, "C3_C7_vs_C0_results_FINAL")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


# --- Sample metadata --------------------------------------------

sample_tbl <- tibble::tribble(
  ~sample_name,     ~sample_dir,     ~sqanti_prefix,
  "C0_Control_1",   "C0_Sciatic_1",  "C0_Sciatic_1",
  "C0_Control_2",   "C0_Sciatic_2",  "C0_Sciatic_2",
  "C0_Control_3",   "C0_Sciatic_3",  "C0_Sciatic_3",
  "C3_Injured_1",   "C3_Injured_1",  "C3_Injured_Sciatic_1",
  "C3_Injured_2",   "C3_Injured_2",  "C3_Injured_Sciatic_2",
  "C3_Injured_3",   "C3_Injured_3",  "C3_Injured_Sciatic_3",
  "C7_Injured_1",   "C7_Injured_1",  "C7_Injured_Sciatic_1",
  "C7_Injured_2",   "C7_Injured_2",  "C7_Injured_Sciatic_2",
  "C7_Injured_3",   "C7_Injured_3",  "C7_Injured_Sciatic_3"
) %>%
  mutate(
    sqanti_dir       = file.path(data_dir, sample_dir, "SQANTI3_Results"),
    abundance_file   = file.path(sqanti_dir, "abundance.tsv"),
    class_file       = file.path(
      sqanti_dir,
      paste0(sqanti_prefix, "_filtered_RulesFilter_result_classification.txt")
    ),
    rescue_file      = file.path(
      sqanti_dir,
      paste0(sqanti_prefix, "_rescue_automatic_inclusion_list.tsv")
    ),
    rescued_class_file = file.path(
      sqanti_dir,
      paste0(sqanti_prefix, "_filtered_RulesFilter_result_classification_with_rescue_flag.txt")
    ),
    unified_abundance_file = file.path(
      results_dir,
      paste0(sample_name, "_abundance_unified.tsv")
    )
  )


# --- Helper functions -------------------------------------------

check_files_exist <- function(paths, label) {
  missing_paths <- paths[!file.exists(paths)]
  
  if (length(missing_paths) > 0) {
    stop(
      "Missing ", label, " file(s):\n",
      paste(missing_paths, collapse = "\n"),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


strip_ensembl_version <- function(x) {
  # Use only for lookup keys. Do not use this to overwrite original IDs.
  stringr::str_remove(x, "\\.\\d+$")
}


read_rescue_ids <- function(path) {
  readr::read_tsv(path, col_names = FALSE, show_col_types = FALSE)[[1]] %>%
    unique()
}


add_rescue_flag <- function(class_path, rescue_path, out_path) {
  cls <- readr::read_tsv(class_path, show_col_types = FALSE)
  rescued_ids <- read_rescue_ids(rescue_path)
  
  cls_rescued <- cls %>%
    mutate(
      filter_rescued_result = if_else(
        filter_result == "Isoform" | isoform %in% rescued_ids,
        "Isoform",
        "Artifact"
      )
    )
  
  readr::write_tsv(cls_rescued, out_path)
  message("Saved rescue-flagged classification: ", out_path)
  
  invisible(out_path)
}


unify_sample_abundance <- function(abundance_file,
                                   rescued_class_file,
                                   lookup_tbl,
                                   sample_index,
                                   out_file) {
  valid_isoforms <- readr::read_tsv(rescued_class_file, show_col_types = FALSE) %>%
    filter(filter_rescued_result == "Isoform") %>%
    pull(isoform) %>%
    unique()
  
  abundance_unified <- readr::read_tsv(abundance_file, show_col_types = FALSE) %>%
    mutate(
      base_id = target_id,
      lookup_id = paste0(base_id, "_", sample_index)
    ) %>%
    filter(base_id %in% valid_isoforms) %>%
    left_join(lookup_tbl, by = c("lookup_id" = "transcript_id")) %>%
    mutate(
      # Preserve original SQANTI3/IsoQuant transcript IDs unless a unified
      # novel-isoform ID is available from the lookup table.
      transcript_id = coalesce(unified_id, base_id)
    ) %>%
    select(transcript_id, everything(), -unified_id, -lookup_id, -base_id)
  
  readr::write_tsv(abundance_unified, out_file)
  message("Wrote unified abundance file: ", out_file, " (", nrow(abundance_unified), " isoforms)")
  
  invisible(out_file)
}


# --- Input checks -----------------------------------------------

check_files_exist(sample_tbl$abundance_file, "abundance")
check_files_exist(sample_tbl$class_file, "SQANTI3 classification")
check_files_exist(sample_tbl$rescue_file, "rescue inclusion-list")

lookup_file <- file.path(results_dir, "novel_isoform_translation_FINAL.csv")
check_files_exist(lookup_file, "novel isoform translation lookup")


# --- Step 1: Add rescue flags to SQANTI3 classification files ----

pwalk(
  sample_tbl %>% select(class_file, rescue_file, rescued_class_file),
  add_rescue_flag
)


# --- Step 2: Filter and unify each sample abundance file ---------

lookup_tbl <- readr::read_csv(lookup_file, show_col_types = FALSE) %>%
  select(unified_id, transcript_id)

pwalk(
  sample_tbl %>%
    mutate(sample_index = row_number()) %>%
    select(
      abundance_file,
      rescued_class_file,
      sample_index,
      unified_abundance_file
    ),
  function(abundance_file, rescued_class_file, sample_index, unified_abundance_file) {
    unify_sample_abundance(
      abundance_file      = abundance_file,
      rescued_class_file  = rescued_class_file,
      lookup_tbl          = lookup_tbl,
      sample_index        = sample_index,
      out_file            = unified_abundance_file
    )
  }
)


# --- Step 3: Merge unified abundance tables ----------------------
# Change `est_counts` to `tpm` below if you want a TPM matrix instead.

sample_dfs <- map2(
  sample_tbl$unified_abundance_file,
  sample_tbl$sample_name,
  function(path, sample_name) {
    readr::read_tsv(path, show_col_types = FALSE) %>%
      filter(!str_starts(transcript_id, "transcript")) %>%
      select(transcript_id, est_counts) %>%
      rename(!!sample_name := est_counts)
  }
)

combined_counts <- reduce(sample_dfs, full_join, by = "transcript_id") %>%
  mutate(across(-transcript_id, as.numeric)) %>%
  replace_na(as.list(setNames(rep(0, nrow(sample_tbl)), sample_tbl$sample_name)))


# --- Step 4: Extract transcript-to-gene mappings -----------------

class_df <- map_dfr(sample_tbl$rescued_class_file, readr::read_tsv, show_col_types = FALSE) %>%
  filter(filter_rescued_result == "Isoform") %>%
  transmute(
    transcript_id = isoform,
    associated_gene = associated_gene
  ) %>%
  distinct(transcript_id, associated_gene)

annotated <- combined_counts %>%
  left_join(class_df, by = "transcript_id") %>%
  mutate(
    associated_gene = strip_ensembl_version(associated_gene),
    associated_gene = if_else(
      str_starts(transcript_id, "ENSMUST"),
      associated_gene,
      str_remove(transcript_id, "-novel.*$")
    )
  )


# --- Step 5: Add gene symbols from BioMart -----------------------

genes_to_lookup <- annotated %>%
  filter(str_starts(associated_gene, "ENSMUSG")) %>%
  pull(associated_gene) %>%
  unique()

ensembl <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")

gene_map <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters    = "ensembl_gene_id",
  values     = genes_to_lookup,
  mart       = ensembl
) %>%
  as_tibble() %>%
  rename(
    associated_gene = ensembl_gene_id,
    gene_symbol     = external_gene_name
  )

final_tbl <- annotated %>%
  left_join(gene_map, by = "associated_gene") %>%
  mutate(gene_symbol = coalesce(gene_symbol, associated_gene))


# --- Step 6: Add transcript symbols from BioMart -----------------
final_tbl <- final_tbl %>%
  mutate(
    transcript_core_for_biomart = if_else(
      str_starts(transcript_id, "ENSMUST"),
      strip_ensembl_version(transcript_id),
      NA_character_
    )
  )

tx_map <- biomaRt::getBM(
  attributes = c("ensembl_transcript_id", "external_transcript_name"),
  filters    = "ensembl_transcript_id",
  values     = unique(na.omit(final_tbl$transcript_core_for_biomart)),
  mart       = ensembl
) %>%
  as_tibble() %>%
  rename(
    transcript_core_for_biomart = ensembl_transcript_id,
    transcript_symbol           = external_transcript_name
  )

final_tbl <- final_tbl %>%
  left_join(tx_map, by = "transcript_core_for_biomart") %>%
  mutate(
    transcript_symbol = coalesce(transcript_symbol, transcript_id),
    gene_symbol = if_else(
      is.na(gene_symbol) | gene_symbol == "",
      associated_gene,
      gene_symbol
    )
  ) %>%
  select(
    transcript_id,
    transcript_symbol,
    associated_gene,
    gene_symbol,
    everything(),
    -transcript_core_for_biomart
  )


# --- Step 7: Fill remaining missing gene symbols/IDs -------------

final_tbl <- final_tbl %>%
  mutate(
    transcript_symbol = if_else(
      is.na(transcript_symbol) | transcript_symbol == "",
      transcript_id,
      transcript_symbol
    ),
    gene_symbol = if_else(
      is.na(gene_symbol) | gene_symbol == "",
      str_remove(transcript_symbol, "-[0-9]+$"),
      gene_symbol
    ),
    associated_gene = na_if(trimws(associated_gene), "")
  ) %>%
  group_by(gene_symbol) %>%
  mutate(
    associated_gene = if_else(
      is.na(associated_gene),
      first(associated_gene[!is.na(associated_gene)]),
      associated_gene
    )
  ) %>%
  ungroup()


# Reverse BioMart lookup for rows where associated_gene is a symbol

symbols_to_lookup <- final_tbl %>%
  filter(is.na(associated_gene) | !str_starts(associated_gene, "ENSMUSG")) %>%
  filter(!str_starts(coalesce(associated_gene, ""), "novel_gene_")) %>%
  pull(associated_gene) %>%
  na.omit() %>%
  unique()

symbol_map <- biomaRt::getBM(
  attributes = c("external_gene_name", "ensembl_gene_id"),
  filters    = "external_gene_name",
  values     = symbols_to_lookup,
  mart       = ensembl
) %>%
  as_tibble() %>%
  rename(
    associated_gene_symbol = external_gene_name,
    ensmusg_id             = ensembl_gene_id
  ) %>%
  distinct(associated_gene_symbol, .keep_all = TRUE)

final_tbl <- final_tbl %>%
  left_join(symbol_map, by = c("associated_gene" = "associated_gene_symbol")) %>%
  mutate(
    associated_gene = if_else(
      !str_starts(coalesce(associated_gene, ""), "ENSMUSG") & !is.na(ensmusg_id),
      ensmusg_id,
      associated_gene
    )
  ) %>%
  select(-ensmusg_id)


# --- Step 8: Save final annotated count table --------------------

out_file <- file.path(results_dir, "all_transcripts_with_associated_genes_count.tsv")
readr::write_tsv(final_tbl, out_file)

