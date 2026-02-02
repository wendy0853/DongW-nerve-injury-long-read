# Set working directories
setwd("/Users/wendydong/Documents/WUSM PhD/Long Read RNA Seq/Analysis_Long_Short_updated/")
data_dir <- "/Users/wendydong/Documents/WUSM PhD/Long Read RNA Seq/Long_read_results_updated/Filtered_SQANTI3_Results/"
results_dir <- "/Users/wendydong/Documents/WUSM PhD/Long Read RNA Seq/Long_read_results_updated/Filtered_SQANTI3_Results/"
meta_dir <- "/Users/wendydong/Documents/WUSM PhD/Long Read RNA Seq/Analysis_Long_Short_updated/C3_C7_vs_C0_results_FINAL/"

# install.packages(c("rtracklayer","GenomicRanges","dplyr","purrr","tibble","stringr"))
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(purrr)
library(tibble)
library(stringr)
library(biomaRt)
library(readr)
library(tidyr)

# 1. List of your full GTF paths
file_names <- c(
  "C0_Sciatic_1_filtered_RulesFilter_result_classification_with_rescue_flag.txt",
  "C0_Sciatic_2_filtered_RulesFilter_result_classification_with_rescue_flag.txt",
  "C0_Sciatic_3_filtered_RulesFilter_result_classification_with_rescue_flag.txt",
  "C3_Injured_Sciatic_1_filtered_RulesFilter_result_classification_with_rescue_flag.txt",
  "C3_Injured_Sciatic_2_filtered_RulesFilter_result_classification_with_rescue_flag.txt",
  "C3_Injured_Sciatic_3_filtered_RulesFilter_result_classification_with_rescue_flag.txt",
  "C7_Injured_Sciatic_1_filtered_RulesFilter_result_classification_with_rescue_flag.txt",
  "C7_Injured_Sciatic_2_filtered_RulesFilter_result_classification_with_rescue_flag.txt",
  "C7_Injured_Sciatic_3_filtered_RulesFilter_result_classification_with_rescue_flag.txt"
)
# Full paths to classification.txt files
file_paths <- file.path(data_dir, file_names)

ens_map <- read_tsv(file.path(meta_dir, "all_transcripts_with_associated_genes_count.tsv")) %>%
  dplyr::select(transcript_id, transcript_symbol) %>% distinct()

novel_map <- read_csv(file.path(meta_dir, "novel_isoform_translation.csv")) %>%
  dplyr::select(transcript_id, unified_id) %>% distinct()

add_transcript_id <- function(classfile, sample_idx, results_dir) {
  df <- read_tsv(classfile, show_col_types = FALSE)
  
  df <- df %>%
    mutate(
      # Only strip version for Ensembl transcripts
      isoform_nover = if_else(
        str_starts(isoform, "ENSMUST"),
        str_remove(isoform, "\\.\\d+$"),
        isoform
      ),
      
      transcript_id_full = if_else(
        str_starts(isoform_nover, "transcript"),
        paste0(isoform_nover, "_", sample_idx),
        isoform_nover
      )
    ) %>%
    left_join(novel_map, by = c("transcript_id_full" = "transcript_id")) %>%
    left_join(ens_map,   by = c("isoform_nover" = "transcript_id")) %>%
    mutate(
      transcript_id = coalesce(unified_id, transcript_symbol, isoform_nover)
    ) %>%
    dplyr::select(
      transcript_id,
      everything(),
      -unified_id, -transcript_symbol, -transcript_id_full, -isoform_nover
    )
  
  fname_out <- file.path(
    results_dir,
    paste0(str_replace(basename(classfile), "\\.txt$", ""), "_with_transcriptID.txt")
  )
  write_tsv(df, fname_out)
  message("Saved: ", fname_out)
}

# Loop over files
walk2(file_paths, seq_along(file_paths), ~add_transcript_id(.x, .y, results_dir))
