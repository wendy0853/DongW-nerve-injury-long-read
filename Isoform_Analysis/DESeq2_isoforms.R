# Set working directories
setwd("/Users/wendydong/Documents/WUSM PhD/Long Read RNA Seq/Analysis_Long_Short_updated/")
data_dir <- "/Users/wendydong/Documents/WUSM PhD/Long Read RNA Seq/Long_read_results_updated/Filtered_SQANTI3_Results/"
results_dir <- "/Users/wendydong/Documents/WUSM PhD/Long Read RNA Seq/Analysis_Long_Short_updated/sample_metrics/"

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

# 1. List of your full Classification Files
file_names <- c(
  "C0_Sciatic_1_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt",
  "C0_Sciatic_2_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt",
  "C0_Sciatic_3_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt",
  "C3_Injured_Sciatic_1_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt",
  "C3_Injured_Sciatic_2_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt",
  "C3_Injured_Sciatic_3_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt",
  "C7_Injured_Sciatic_1_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt",
  "C7_Injured_Sciatic_2_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt",
  "C7_Injured_Sciatic_3_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt"
)

# These are under a common data_dir, prepend it:
file_paths <- file.path(data_dir, file_names)

# --- TPM Filtering (≥1 TPM in ≥2 samples) -----------------------
tpm <- read_tsv("/Users/wendydong/Documents/WUSM PhD/Long Read RNA Seq/Analysis_Long_Short_updated/C3_C7_vs_C0_results_FINAL/all_transcripts_with_associated_genes_tpm.tsv", show_col_types = FALSE)
tpm <- tpm %>%
  mutate(transcript_symbol = if_else(is.na(transcript_symbol) | transcript_symbol == "",transcript_id,transcript_symbol ))
tpm_filtered <- tpm %>%
  mutate(samples_ge1 = rowSums(across(5:ncol(.), ~ .x >= 1))) %>%
  filter(samples_ge1 >= 2)
filtered_ids <- tpm_filtered$transcript_symbol

# ---- Read + filter classification consistently (Isoform + TPM-filtered IDs) ----
read_class_filtered <- function(file, filtered_ids) {
  df <- read_tsv(file, show_col_types = FALSE)
   
  # If transcript_symbol column doesn't exist, fail early with a clear message
  if (!("transcript_id" %in% colnames(df))) {
    stop("Missing column 'transcript_id' in: ", file,
         "\nAvailable columns:\n", paste(colnames(df), collapse = ", "))
  }
  
  df %>%
    filter(filter_rescued_result == "Isoform") %>%
    filter(transcript_id %in% filtered_ids)
}


# Helper to extract sample name
get_sample_name <- function(path) {
  basename(path) %>%
    str_replace("_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt", "")
}

# Summarize QC Metrics for one file
summarize_result <- function(file) {
  
  # NEW: read + filter once
  true_isoforms <- read_class_filtered(file, filtered_ids)
  
  # If you still want totals BEFORE filtering, read raw separately (optional)
  df_raw <- read_tsv(file, show_col_types = FALSE)
  
  # Totals (raw)
  n_total <- nrow(df_raw)
  n_iso_raw <- sum(df_raw$filter_result == "Isoform", na.rm = TRUE)
  n_art_raw <- sum(df_raw$filter_result != "Isoform", na.rm = TRUE)
  
  # Totals AFTER your combined filters (Isoform + TPM list)
  n_iso <- nrow(true_isoforms)
  
  # --- Known vs Novel genes (based on filtered true_isoforms) ---
  n_total_genes <- true_isoforms %>%
    distinct(associated_gene) %>%
    nrow()
  
  n_known <- true_isoforms %>%
    filter(!is.na(associated_gene), associated_gene != "") %>%
    filter(str_detect(associated_gene, "^ENSMUSG")) %>%
    distinct(associated_gene) %>%
    nrow()
  
  n_novel <- true_isoforms %>%
    filter(!is.na(associated_gene), associated_gene != "") %>%
    filter(!str_detect(associated_gene, "^ENSMUSG")) %>%
    distinct(associated_gene) %>%
    nrow()
  
  denom_genes <- n_known + n_novel
  
  gene_summary <- tibble(
    Sample_ID        = get_sample_name(file),
    Detected_Genes   = n_total_genes,
    Known_Genes      = n_known,
    Known_Genes_Pct  = ifelse(denom_genes > 0, round(100 * n_known / denom_genes, 2), NA_real_),
    Novel_Genes      = n_novel,
    Novel_Genes_Pct  = ifelse(denom_genes > 0, round(100 * n_novel / denom_genes, 2), NA_real_)
  )
  
  # --- Isoform vs Artifact (raw, so you can still report filtering impact) ---
  isoform_summary <- tibble(
    Total_Isoforms_raw     = n_total,
    True_Isoform_Count_raw = n_iso_raw,
    Pct_Isoform_raw        = ifelse(n_total > 0, round(100 * n_iso_raw / n_total, 2), NA_real_),
    Artifact_Count_raw     = n_art_raw,
    Pct_Artifact_raw       = ifelse(n_total > 0, round(100 * n_art_raw / n_total, 2), NA_real_),
    
    # NEW: post-filter count (Isoform + TPM)
    True_Isoform_Count_postTPM = n_iso
  )
  
  # --- Structural categories (based on filtered true_isoforms) ---
  expected_cats <- c("full-splice_match", "incomplete-splice_match",
                     "novel_in_catalog", "novel_not_in_catalog","genic",
                     "antisense", "fusion", "intergenic","genic_intron")
  
  cat_counts <- true_isoforms %>%
    group_by(structural_category) %>%
    tally(name = "Count") %>%
    ungroup() %>%
    complete(structural_category = expected_cats, fill = list(Count = 0)) %>%
    pivot_wider(names_from = structural_category, values_from = Count) %>%
    rename(
      FSM          = `full-splice_match`,
      ISM          = `incomplete-splice_match`,
      NIC          = `novel_in_catalog`,
      NNC          = `novel_not_in_catalog`,
      Genic        = genic,
      Antisense    = antisense,
      Fusion       = fusion,
      Intergenic   = intergenic,
      Genic_Intron = genic_intron
    )
  
  cat_summary <- cat_counts %>%
    mutate(
      FSM_Pct        = ifelse(n_iso > 0, round(100 * FSM / n_iso, 2), NA_real_),
      ISM_Pct        = ifelse(n_iso > 0, round(100 * ISM / n_iso, 2), NA_real_),
      NIC_Pct        = ifelse(n_iso > 0, round(100 * NIC / n_iso, 2), NA_real_),
      NNC_Pct        = ifelse(n_iso > 0, round(100 * NNC / n_iso, 2), NA_real_),
      Genic_Pct      = ifelse(n_iso > 0, round(100 * Genic / n_iso, 2), NA_real_),
      Antisense_Pct  = ifelse(n_iso > 0, round(100 * Antisense / n_iso, 2), NA_real_),
      Fusion_Pct     = ifelse(n_iso > 0, round(100 * Fusion / n_iso, 2), NA_real_),
      Intergenic_Pct = ifelse(n_iso > 0, round(100 * Intergenic / n_iso, 2), NA_real_),
      Genic_Intron_Pct = ifelse(n_iso > 0, round(100 * Genic_Intron / n_iso, 2), NA_real_)
    ) %>%
    dplyr::select(FSM, FSM_Pct, ISM, ISM_Pct, NIC, NIC_Pct, NNC, NNC_Pct,
           Genic, Genic_Pct, Antisense, Antisense_Pct,
           Fusion, Fusion_Pct, Intergenic, Intergenic_Pct,
           Genic_Intron, Genic_Intron_Pct)
  
  # --- Coding vs non-coding (filtered true_isoforms) ---
  n_coding <- sum(true_isoforms$coding == "coding", na.rm = TRUE)
  n_non_coding <- sum(true_isoforms$coding == "non_coding", na.rm = TRUE)
  
  coding_summary <- tibble(
    Protein_Coding_Counts     = n_coding,
    Protein_Coding_Pct        = ifelse(n_iso > 0, round(100 * n_coding / n_iso, 2), NA_real_),
    Protein_Non_Coding_Counts = n_non_coding,
    Protein_Non_Coding_Pct    = ifelse(n_iso > 0, round(100 * n_non_coding / n_iso, 2), NA_real_)
  )
  
  # --- QC summary (filtered true_isoforms) ---
  qc_summary <- tibble(
    Canonical_Junctions_Pct = ifelse(n_iso > 0,
                                     round(100 * sum(true_isoforms$all_canonical == "canonical", na.rm = TRUE) / n_iso, 2), NA_real_
    ),
    Noncanonical_Junction_Count = sum(true_isoforms$all_canonical != "canonical", na.rm = TRUE)
  )
  
  min_cov_summary <- tibble(
    Mean_MinCov   = ifelse(n_iso > 0, mean(as.numeric(true_isoforms$min_cov), na.rm = TRUE), NA_real_),
    Median_MinCov = ifelse(n_iso > 0, median(as.numeric(true_isoforms$min_cov), na.rm = TRUE), NA_real_)
  )
  
  bind_cols(gene_summary, isoform_summary, cat_summary, coding_summary, qc_summary, min_cov_summary)
}


# Run on all samples
summary <- map_dfr(file_paths, summarize_result)
View(summary)

write_csv(summary, file.path(results_dir, "Filtered_SQANTI3_QC_summary.csv"))



# --- Sample Metric Plot Summary ---
# Summarize Sample Metrics for plotting in one file
summarize_result <- function(file) {
  true_isoforms <- read_class_filtered(file, filtered_ids)
  n_iso <- nrow(true_isoforms)
  
  
  # --- Known vs Novel genes (unique ENSMUSG vs novel) ---
  n_total_genes <- true_isoforms %>%
    distinct(associated_gene) %>%
    nrow()
  
  n_known <- true_isoforms %>%
    filter(!is.na(associated_gene), associated_gene != "") %>%
    filter(str_detect(associated_gene, "^ENSMUSG")) %>%
    distinct(associated_gene) %>%
    nrow()
  
  n_novel <- true_isoforms %>%
    filter(!is.na(associated_gene), associated_gene != "") %>%
    filter(!str_detect(associated_gene, "^ENSMUSG")) %>%
    distinct(associated_gene) %>%
    nrow()
  
  gene_summary <- tibble(
    Sample_ID          = get_sample_name(file),
    Unique_Genes   = n_total_genes,
    Unique_Isoforms = n_iso,
    Annotated_Genes      = n_known,
    Annotated_Genes_Pct  = ifelse(n_total_genes > 0, round(100 * n_known / n_total_genes, 2), NA_real_),
    Novel_Genes      = n_novel,
    Novel_Genes_Pct  = ifelse(n_total_genes > 0, round(100 * n_novel / n_total_genes, 2), NA_real_)
  )
  
  # --- Protein Coding vs. Non-coding (counts + %) ---
  n_coding <- sum(true_isoforms$coding == "coding", na.rm = TRUE)
  n_non_coding <- sum(true_isoforms$coding == "non_coding", na.rm = TRUE)
  
  # Add percentage columns right after each count
  coding_summary <- tibble(
    Protein_Coding_Transcripts = n_coding,
    Protein_Coding_Transcripts_Pct = ifelse(n_iso > 0, round(100 * n_coding / n_iso, 2), NA_real_),
    Protein_Non_Coding_Transcripts = n_non_coding,
    Protein_Non_Coding_Transcripts_Pct = ifelse(n_iso > 0, round(100 * n_non_coding / n_iso, 2), NA_real_)
  )
  
  
  # --- Isoform Annotated vs. Novel (counts + %) ---
  n_t_known <- true_isoforms %>%
    filter(!is.na(associated_transcript), associated_transcript != "") %>%
    filter(str_detect(associated_transcript, "^ENSMUST")) %>%
    distinct(associated_transcript) %>%
    nrow()
  
  n_t_novel <- true_isoforms %>%
    filter(!is.na(associated_transcript), associated_transcript != "") %>%
    filter(str_detect(associated_transcript, "novel")) %>%
    nrow()
  
  isoform_summary <- tibble(
    Annotated_Transcripts      = n_t_known,
    Annotated_Transcripts_Pct  = ifelse(n_iso > 0, round(100 * n_t_known / n_iso, 2), NA_real_),
    Novel_Transcripts      = n_t_novel,
    Novel_Transcripts_Pct  = ifelse(n_iso > 0, round(100 * n_t_novel / n_iso, 2), NA_real_)
  )
  
  # --- Structural categories (counts + %) ---
  expected_cats <- c("full-splice_match", "incomplete-splice_match",
                     "novel_in_catalog", "novel_not_in_catalog","genic",
                     "antisense", "fusion", "intergenic","genic_intron")
  
  # Tally, complete with zeros for missing categories
  cat_counts <- true_isoforms %>%
    group_by(structural_category) %>%
    tally(name = "Count") %>%
    ungroup() %>%
    complete(structural_category = expected_cats, fill = list(Count = 0)) %>%
    pivot_wider(names_from = structural_category,
                values_from = Count) %>%
    # Rename to simpler labels
    rename(
      FSM          = `full-splice_match`,
      ISM          = `incomplete-splice_match`,
      NIC          = `novel_in_catalog`,
      NNC          = `novel_not_in_catalog`,
      Genic        = genic,
      Antisense    = antisense,
      Fusion       = fusion,
      Intergenic   = intergenic,
      Genic_Intron = genic_intron
    )
  
  # Add percentage columns right after each count
  cat_summary <- cat_counts %>%
    mutate(
      FSM_Pct        = ifelse(n_iso > 0, round(100 * FSM / n_iso, 2), NA_real_),
      ISM_Pct        = ifelse(n_iso > 0, round(100 * ISM / n_iso, 2), NA_real_),
      NIC_Pct        = ifelse(n_iso > 0, round(100 * NIC / n_iso, 2), NA_real_),
      NNC_Pct        = ifelse(n_iso > 0, round(100 * NNC / n_iso, 2), NA_real_),
      Genic_Pct  = ifelse(n_iso > 0, round(100 * Genic / n_iso, 2), NA_real_),
      Antisense_Pct  = ifelse(n_iso > 0, round(100 * Antisense / n_iso, 2), NA_real_),
      Fusion_Pct  = ifelse(n_iso > 0, round(100 * Fusion / n_iso, 2), NA_real_),
      Intergenic_Pct = ifelse(n_iso > 0, round(100 * Intergenic / n_iso, 2), NA_real_),
      Genic_Intron_Pct = ifelse(n_iso > 0, round(100 * Genic_Intron / n_iso, 2), NA_real_)
    ) %>%
    dplyr::select(FSM, FSM_Pct,
                  ISM, ISM_Pct,
                  NIC, NIC_Pct,
                  NNC, NNC_Pct,
                  Genic, Genic_Pct,
                  Antisense, Antisense_Pct,
                  Fusion, Fusion_Pct,
                  Intergenic, Intergenic_Pct,
                  Genic_Intron, Genic_Intron_Pct)
  

  
  # --- Final combined row ---
  bind_cols(gene_summary, coding_summary, isoform_summary, cat_summary)
}

# Run on all samples
summary <- map_dfr(file_paths, summarize_result)
View(summary)

write_csv(summary, file.path(results_dir, "Filtered_SQANTI3_Plot_summary.csv"))






# --- Sample Metric Plot Summary ---
# Summarize Sample Metrics for plotting in one file
summarize_result <- function(file) {
  true_isoforms <- read_class_filtered(file, filtered_ids)
  n_iso <- nrow(true_isoforms)
  
  expected_cats <- c("full-splice_match", "incomplete-splice_match",
                     "novel_in_catalog", "novel_not_in_catalog",
                     "genic", "antisense", "fusion", "intergenic", "genic_intron")
  
  # Count first (robust)
  counts_long <- true_isoforms %>%
    mutate(
      structural_category = str_trim(as.character(structural_category)),
      coding_group = ifelse(coding == "coding", "coding", "non_coding")
    ) %>%
    count(structural_category, coding_group, name = "Count")
  
  wide <- counts_long %>%
    complete(structural_category = expected_cats,
             coding_group = c("coding", "non_coding"),
             fill = list(Count = 0)) %>%
    pivot_wider(names_from = c(structural_category, coding_group),
                values_from = Count,
                values_fill = 0)
  
  # Define the exact wide columns you want (counts)
  count_cols <- c(
    "full-splice_match_coding","full-splice_match_non_coding",
    "incomplete-splice_match_coding","incomplete-splice_match_non_coding",
    "novel_in_catalog_coding","novel_in_catalog_non_coding",
    "novel_not_in_catalog_coding","novel_not_in_catalog_non_coding",
    "genic_coding","genic_non_coding",
    "antisense_coding","antisense_non_coding",
    "fusion_coding","fusion_non_coding",
    "intergenic_coding","intergenic_non_coding",
    "genic_intron_coding","genic_intron_non_coding"
  )
  
  # Ensure all exist (some samples won't have some categories)
  for (nm in count_cols) if (!nm %in% names(wide)) wide[[nm]] <- 0L
  
  # Percent columns for THESE ONLY
  wide <- wide %>%
    mutate(across(
      all_of(count_cols),
      ~ ifelse(n_iso > 0, round(100 * .x / n_iso, 2), NA_real_),
      .names = "{.col}_pct"
    ))
  
  ordered_cols <- c(rbind(count_cols, paste0(count_cols, "_pct")))
  
  bind_cols(
    tibble(Sample_ID = get_sample_name(file), n_iso = n_iso),
    wide %>% dplyr::select(all_of(ordered_cols))
  )
}


# Run on all samples
summary <- purrr::map_dfr(file_paths, summarize_result)
View(summary)
write_csv(summary, file.path(results_dir, "Filtered_SQANTI3_Category_summary.csv"))


