# ---------------------------------------------------------------
# SQANTI3 Long-Read RNA-seq Isoform Processing and Annotation
# Author: Wendy Dong
# Affiliation: Washington University School of Medicine
# Description: Filters SQANTI3 isoforms, unifies transcript IDs,
#              merges abundance counts, and annotates genes/transcripts.
# ---------------------------------------------------------------

# --- Setup -----------------------------------------------------
setwd("/Users/wendydong/Documents/WUSM PhD/Long Read RNA Seq/Analysis_Long_Short_updated/")

data_dir    <- "/Users/wendydong/Documents/WUSM PhD/Long Read RNA Seq/Long_read_results_updated/"
results_dir <- file.path(getwd(), "C3_C7_vs_C0_results_FINAL")

# Load libraries
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(purrr)
library(tibble)
library(stringr)
library(biomaRt)
library(readr)
library(tidyr)

# --- File Definitions ------------------------------------------
sample_names <- c(
  "C0_Control_1", "C0_Control_2", "C0_Control_3",
  "C3_Injured_1", "C3_Injured_2", "C3_Injured_3",
  "C7_Injured_1", "C7_Injured_2", "C7_Injured_3"
)

ab_files <- file.path(data_dir, c(
  "C0_Sciatic_1/SQANTI3_Results/abundance.tsv",
  "C0_Sciatic_2/SQANTI3_Results/abundance.tsv",
  "C0_Sciatic_3/SQANTI3_Results/abundance.tsv",
  "C3_Injured_1/SQANTI3_Results/abundance.tsv",
  "C3_Injured_2/SQANTI3_Results/abundance.tsv",
  "C3_Injured_3/SQANTI3_Results/abundance.tsv",
  "C7_Injured_1/SQANTI3_Results/abundance.tsv",
  "C7_Injured_2/SQANTI3_Results/abundance.tsv",
  "C7_Injured_3/SQANTI3_Results/abundance.tsv"
))

class_files <- file.path(data_dir, c(
  "C0_Sciatic_1/SQANTI3_Results/C0_Sciatic_1_filtered_RulesFilter_result_classification.txt",
  "C0_Sciatic_2/SQANTI3_Results/C0_Sciatic_2_filtered_RulesFilter_result_classification.txt",
  "C0_Sciatic_3/SQANTI3_Results/C0_Sciatic_3_filtered_RulesFilter_result_classification.txt",
  "C3_Injured_1/SQANTI3_Results/C3_Injured_Sciatic_1_filtered_RulesFilter_result_classification.txt",
  "C3_Injured_2/SQANTI3_Results/C3_Injured_Sciatic_2_filtered_RulesFilter_result_classification.txt",
  "C3_Injured_3/SQANTI3_Results/C3_Injured_Sciatic_3_filtered_RulesFilter_result_classification.txt",
  "C7_Injured_1/SQANTI3_Results/C7_Injured_Sciatic_1_filtered_RulesFilter_result_classification.txt",
  "C7_Injured_2/SQANTI3_Results/C7_Injured_Sciatic_2_filtered_RulesFilter_result_classification.txt",
  "C7_Injured_3/SQANTI3_Results/C7_Injured_Sciatic_3_filtered_RulesFilter_result_classification.txt"
))

# --- File Definitions ------------------------------------------
sample_names <- c(
  "C0_Control_1", "C0_Control_2", "C0_Control_3",
  "C3_Injured_1", "C3_Injured_2", "C3_Injured_3",
  "C7_Injured_1", "C7_Injured_2", "C7_Injured_3"
)

ab_files <- file.path(data_dir, c(
  "C0_Sciatic_1/SQANTI3_Results/abundance.tsv",
  "C0_Sciatic_2/SQANTI3_Results/abundance.tsv",
  "C0_Sciatic_3/SQANTI3_Results/abundance.tsv",
  "C3_Injured_1/SQANTI3_Results/abundance.tsv",
  "C3_Injured_2/SQANTI3_Results/abundance.tsv",
  "C3_Injured_3/SQANTI3_Results/abundance.tsv",
  "C7_Injured_1/SQANTI3_Results/abundance.tsv",
  "C7_Injured_2/SQANTI3_Results/abundance.tsv",
  "C7_Injured_3/SQANTI3_Results/abundance.tsv"
))

class_files <- file.path(data_dir, c(
  "C0_Sciatic_1/SQANTI3_Results/C0_Sciatic_1_filtered_RulesFilter_result_classification.txt",
  "C0_Sciatic_2/SQANTI3_Results/C0_Sciatic_2_filtered_RulesFilter_result_classification.txt",
  "C0_Sciatic_3/SQANTI3_Results/C0_Sciatic_3_filtered_RulesFilter_result_classification.txt",
  "C3_Injured_1/SQANTI3_Results/C3_Injured_Sciatic_1_filtered_RulesFilter_result_classification.txt",
  "C3_Injured_2/SQANTI3_Results/C3_Injured_Sciatic_2_filtered_RulesFilter_result_classification.txt",
  "C3_Injured_3/SQANTI3_Results/C3_Injured_Sciatic_3_filtered_RulesFilter_result_classification.txt",
  "C7_Injured_1/SQANTI3_Results/C7_Injured_Sciatic_1_filtered_RulesFilter_result_classification.txt",
  "C7_Injured_2/SQANTI3_Results/C7_Injured_Sciatic_2_filtered_RulesFilter_result_classification.txt",
  "C7_Injured_3/SQANTI3_Results/C7_Injured_Sciatic_3_filtered_RulesFilter_result_classification.txt"
))


rescue_inclusion_list <- file.path(data_dir, c(
  "C0_Sciatic_1/SQANTI3_Results/C0_Sciatic_1_rescue_automatic_inclusion_list.tsv",
  "C0_Sciatic_2/SQANTI3_Results/C0_Sciatic_2_rescue_automatic_inclusion_list.tsv",
  "C0_Sciatic_3/SQANTI3_Results/C0_Sciatic_3_rescue_automatic_inclusion_list.tsv",
  "C3_Injured_1/SQANTI3_Results/C3_Injured_Sciatic_1_rescue_automatic_inclusion_list.tsv",
  "C3_Injured_2/SQANTI3_Results/C3_Injured_Sciatic_2_rescue_automatic_inclusion_list.tsv",
  "C3_Injured_3/SQANTI3_Results/C3_Injured_Sciatic_3_rescue_automatic_inclusion_list.tsv",
  "C7_Injured_1/SQANTI3_Results/C7_Injured_Sciatic_1_rescue_automatic_inclusion_list.tsv",
  "C7_Injured_2/SQANTI3_Results/C7_Injured_Sciatic_2_rescue_automatic_inclusion_list.tsv",
  "C7_Injured_3/SQANTI3_Results/C7_Injured_Sciatic_3_rescue_automatic_inclusion_list.tsv"
))

lookup <- read_csv(file.path(results_dir, "novel_isoform_translation.csv"), show_col_types = FALSE) %>%
  dplyr::select(unified_id, transcript_id)

out_files <- paste0(sample_names, "_abundance_unified.tsv")

#stopifnot(length(class_files) == length(rescue_inclusion_list))

walk2(class_files, rescue_inclusion_list, function(class_path, rescue_path) {
  cls <- read_tsv(class_path, show_col_types = FALSE)
  rescued_ids <- read_tsv(rescue_path,col_names = FALSE,show_col_types = FALSE)[[1]]
  cls <- cls %>%
    mutate(
      filter_rescued_result = if_else(filter_result == "Isoform" | isoform %in% rescued_ids,"Isoform","Artifact" )
    )
  # Write new file
  out_path <- sub("\\.txt$", "_with_rescue_flag.txt", class_path)
  write_tsv(cls, out_path)
  message("Saved: ", out_path)
})

# --- Step 1: Filter and unify each sample’s abundance file -------
class_files <- file.path(data_dir, c(
  "C0_Sciatic_1/SQANTI3_Results/C0_Sciatic_1_filtered_RulesFilter_result_classification_with_rescue_flag.txt",
  "C0_Sciatic_2/SQANTI3_Results/C0_Sciatic_2_filtered_RulesFilter_result_classification_with_rescue_flag.txt",
  "C0_Sciatic_3/SQANTI3_Results/C0_Sciatic_3_filtered_RulesFilter_result_classification_with_rescue_flag.txt",
  "C3_Injured_1/SQANTI3_Results/C3_Injured_Sciatic_1_filtered_RulesFilter_result_classification_with_rescue_flag.txt",
  "C3_Injured_2/SQANTI3_Results/C3_Injured_Sciatic_2_filtered_RulesFilter_result_classification_with_rescue_flag.txt",
  "C3_Injured_3/SQANTI3_Results/C3_Injured_Sciatic_3_filtered_RulesFilter_result_classification_with_rescue_flag.txt",
  "C7_Injured_1/SQANTI3_Results/C7_Injured_Sciatic_1_filtered_RulesFilter_result_classification_with_rescue_flag.txt",
  "C7_Injured_2/SQANTI3_Results/C7_Injured_Sciatic_2_filtered_RulesFilter_result_classification_with_rescue_flag.txt",
  "C7_Injured_3/SQANTI3_Results/C7_Injured_Sciatic_3_filtered_RulesFilter_result_classification_with_rescue_flag.txt"
))

for (i in seq_along(ab_files)) {
  valid_isoforms <- read_tsv(class_files[i], show_col_types = FALSE) %>%
    filter(filter_rescued_result == "Isoform") %>%
    mutate(transcript_id = str_remove(isoform, "\\.\\d+$")) %>%
    pull(transcript_id) %>%
    unique()
  
  df2 <- read_tsv(ab_files[i], show_col_types = FALSE) %>%
    mutate(base_id = str_remove(target_id, "\\.\\d+$")) %>%
    filter(base_id %in% valid_isoforms) %>%
    mutate(lookup_id = paste0(base_id, "_", i)) %>%
    left_join(lookup, by = c("lookup_id" = "transcript_id")) %>%
    mutate(transcript_id = coalesce(unified_id, base_id)) %>%
    dplyr::select(transcript_id, everything(), -unified_id, -lookup_id, -base_id)
  
  write_tsv(df2, file.path(results_dir, out_files[i]))
  message("✔ Wrote unified file for ", sample_names[i], " (", nrow(df2), " isoforms)")
}

# --- Step 2: Merge all unified abundance tables ------------------
sample_dfs <- map2(file.path(results_dir, out_files), sample_names, function(f, s) {
  read_tsv(f, show_col_types = FALSE) %>%
    filter(!str_starts(transcript_id, "transcript")) %>%
    dplyr::select(transcript_id, est_counts) %>%     # Change to "tpm" for TPM matrix, "est_counts" for Count matrix
    rename(!!s := est_counts)     # Change to "tpm" for TPM matrix, "est_counts" for Count matrix
})

combined <- reduce(sample_dfs, full_join, by = "transcript_id") %>%
  mutate(across(-transcript_id, as.numeric)) %>%
  replace_na(as.list(setNames(rep(0, length(sample_names)), sample_names)))

# --- Step 3: Extract transcript–gene mappings --------------------
class_df <- map_dfr(class_files, ~ read_tsv(.x, show_col_types = FALSE)) %>%
  filter(filter_result == "Isoform") %>%
  mutate(transcript_id = str_remove(isoform, "\\.\\d+$")) %>%
  dplyr::select(transcript_id, associated_gene) %>%
  distinct(transcript_id, associated_gene)

annotated <- combined %>%
  mutate(tmp = str_remove(transcript_id, "\\.\\d+$")) %>%
  left_join(class_df %>% rename(tmp = transcript_id), by = "tmp") %>%
  dplyr::select(-tmp) %>%
  mutate(
    associated_gene = str_remove(associated_gene, "\\.\\d+$"),
    associated_gene = if_else(
      str_starts(transcript_id, "ENSMUST"),
      associated_gene,
      str_remove(transcript_id, "-novel.*$")
    )
  )

# --- Step 4: Add gene symbols from biomaRt -----------------------
genes_to_lookup <- unique(annotated$associated_gene[str_starts(annotated$associated_gene, "ENSMUSG")])
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

gene_map <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = genes_to_lookup,
  mart = ensembl
) %>%
  as_tibble() %>%
  rename(associated_gene = ensembl_gene_id, gene_symbol = external_gene_name)

final_tbl <- annotated %>%
  left_join(gene_map, by = "associated_gene") %>%
  mutate(gene_symbol = coalesce(gene_symbol, associated_gene))

# --- Step 5: Add transcript symbols --------------------------------
final_tbl <- final_tbl %>%
  mutate(transcript_core = if_else(
    str_starts(transcript_id, "ENSMUST"),
    str_remove(transcript_id, "\\.\\d+$"),
    NA_character_
  ))

bm_tx_map <- getBM(
  attributes = c("ensembl_transcript_id", "external_transcript_name"),
  filters = "ensembl_transcript_id",
  values = unique(na.omit(final_tbl$transcript_core)),
  mart = ensembl
) %>%
  as_tibble() %>%
  rename(transcript_core = ensembl_transcript_id, transcript_symbol = external_transcript_name) %>%
  mutate(transcript_core = str_remove(transcript_core, "\\.\\d+$"))

final_tbl <- final_tbl %>%
  left_join(bm_tx_map, by = "transcript_core") %>%
  mutate(
    transcript_symbol = coalesce(transcript_symbol, transcript_id),
    gene_symbol = if_else(is.na(gene_symbol) | gene_symbol == "", associated_gene, gene_symbol)
  ) %>%
  dplyr::select(transcript_id, transcript_symbol, associated_gene, gene_symbol, everything(), -transcript_core)

final_tbl <- final_tbl %>%
  mutate(
    transcript_symbol = if_else(
      is.na(transcript_symbol) | transcript_symbol == "",
      transcript_id,
      transcript_symbol
    ),
    gene_symbol = if_else(
      is.na(gene_symbol) | gene_symbol == "",
      sub("-[0-9]+$", "", transcript_symbol),
      gene_symbol
    )
  )

final_tbl <- final_tbl %>%
  mutate(
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

# --- Save final annotated table -----------------------------------
write_tsv(final_tbl, file.path(results_dir, "all_transcripts_with_associated_genes_count.tsv"))

