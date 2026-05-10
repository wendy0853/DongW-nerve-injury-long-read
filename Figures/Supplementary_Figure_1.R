#!/usr/bin/env Rscript

################################################################################
# Supplementary Figure 1
#
# Purpose:
#   Generate long-read inferred gene-level DESeq2 results and concordance plots
#   comparing short-read gene-level DGE with long-read inferred gene-level DGE.
#
# Panels:
#   S1a: srRNA-seq vs lrRNA-seq gene-level concordance, C3 vs C0
#   S1b: srRNA-seq vs lrRNA-seq gene-level concordance, C7 vs C0
#   S1c: srRNA-seq vs lrRNA-seq gene-level concordance, C7 vs C3
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(ggrepel)
  library(grid)
})

# ==============================================================================
# User-defined directories
# ==============================================================================

sqanti_results_dir <- "/path/to/long_read_sqanti_results"       # <-- MODIFY HERE
short_read_dir <- "/path/to/short_read_results"                 # <-- MODIFY HERE
isoform_results_dir <- "/path/to/isoform_analysis_results"      # <-- MODIFY HERE
figure_dir <- "/path/to/supplementary_figure_output"            # <-- MODIFY HERE

dir.create(isoform_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Input files
# ==============================================================================

classification_files <- c(                                      # <-- MODIFY HERE IF NEEDED
  "Filtered_SQANTI3_Results/C0_Sciatic_1_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt",
  "Filtered_SQANTI3_Results/C0_Sciatic_2_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt",
  "Filtered_SQANTI3_Results/C0_Sciatic_3_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt",
  "Filtered_SQANTI3_Results/C3_Injured_Sciatic_1_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt",
  "Filtered_SQANTI3_Results/C3_Injured_Sciatic_2_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt",
  "Filtered_SQANTI3_Results/C3_Injured_Sciatic_3_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt",
  "Filtered_SQANTI3_Results/C7_Injured_Sciatic_1_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt",
  "Filtered_SQANTI3_Results/C7_Injured_Sciatic_2_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt",
  "Filtered_SQANTI3_Results/C7_Injured_Sciatic_3_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt"
)

classification_paths <- file.path(sqanti_results_dir, classification_files)

missing_class_files <- classification_paths[!file.exists(classification_paths)]

if (length(missing_class_files) > 0) {
  stop(
    "Missing classification file(s):\n",
    paste(missing_class_files, collapse = "\n"),
    call. = FALSE
  )
}

# ==============================================================================
# Parameters
# ==============================================================================

comparison_config <- list(
  C3_vs_C0 = list(
    short_read_file = "C3_Injured_vs_C0_results.csv",
    long_read_file = "C3_vs_C0_inferred_gene_results.csv",
    title = "C3 vs. C0",
    color = "#E69F00",
    output = "FigS1_C3_Scatterplot.png"
  ),
  C7_vs_C0 = list(
    short_read_file = "C7_Injured_vs_C0_results.csv",
    long_read_file = "C7_vs_C0_inferred_gene_results.csv",
    title = "C7 vs. C0",
    color = "#66BD63",
    output = "FigS1_C7_Scatterplot.png"
  ),
  C7_vs_C3 = list(
    short_read_file = "C7_Injured_vs_C3_Injured_results.csv",
    long_read_file = "C7_vs_C3_inferred_gene_results.csv",
    title = "C7 vs. C3",
    color = "#0072B2",
    output = "FigS1_C7C3_Scatterplot.png"
  )
)

# ==============================================================================
# Helper functions
# ==============================================================================

get_sample_name <- function(path) {
  basename(path) %>%
    str_remove("_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt")
}

get_condition <- function(path) {
  case_when(
    str_detect(path, "C0") ~ "C0_Control",
    str_detect(path, "C3") ~ "C3_Injured",
    str_detect(path, "C7") ~ "C7_Injured",
    TRUE ~ NA_character_
  )
}

extract_results <- function(dds, num, den) {
  res <- results(dds, contrast = c("group", num, den)) %>%
    as.data.frame()

  norm_counts <- counts(dds, normalized = TRUE)

  res %>%
    mutate(
      gene_symbol = rownames(res),
      mean_num = rowMeans(norm_counts[, colData(dds)$group == num, drop = FALSE]),
      mean_den = rowMeans(norm_counts[, colData(dds)$group == den, drop = FALSE])
    ) %>%
    arrange(padj)
}

plot_concordance <- function(short_read_file,
                             long_read_file,
                             title,
                             point_color,
                             output_name) {
  sr_path <- file.path(short_read_dir, short_read_file)
  lr_path <- file.path(isoform_results_dir, long_read_file)

  if (!file.exists(sr_path)) stop("Missing short-read file: ", sr_path)
  if (!file.exists(lr_path)) stop("Missing long-read inferred gene file: ", lr_path)

  sr <- readr::read_csv(sr_path, show_col_types = FALSE)
  lr <- readr::read_csv(lr_path, show_col_types = FALSE)

  lr_gene <- lr %>%
    filter(!is.na(gene_symbol), !is.na(log2FoldChange), !is.na(padj)) %>%
    group_by(gene_symbol) %>%
    summarise(
      lr_gene_log2FC = log2FoldChange[which.min(ifelse(is.na(padj), Inf, padj))][1],
      .groups = "drop"
    )

  sr_gene <- sr %>%
    transmute(
      gene_symbol,
      sr_gene_log2FC = log2FoldChange,
      sr_padj = padj
    ) %>%
    filter(!is.na(gene_symbol), !is.na(sr_gene_log2FC), !is.na(sr_padj)) %>%
    distinct(gene_symbol, .keep_all = TRUE)

  plot_df <- inner_join(sr_gene, lr_gene, by = "gene_symbol")

  cor_res <- cor.test(
    plot_df$sr_gene_log2FC,
    plot_df$lr_gene_log2FC,
    method = "spearman",
    exact = FALSE
  )

  r_lab <- sprintf("%.3f", unname(cor_res$estimate))
  p_val <- cor_res$p.value
  p_lab <- if (p_val == 0) {
    " P < 1e-300"
  } else {
    sprintf(" P = %.3e", p_val)
  }

  p <- ggplot(plot_df, aes(x = sr_gene_log2FC, y = lr_gene_log2FC)) +
    geom_point(alpha = 0.25, size = 0.1, color = point_color) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.2, color = "red") +
    geom_hline(yintercept = 0, linewidth = 0.2, color = "grey60") +
    geom_vline(xintercept = 0, linewidth = 0.2, color = "grey60") +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      hjust = 0,
      vjust = 1,
      label = paste0(" Spearman r = ", r_lab, "\n", p_lab),
      size = 2,
      color = "black"
    ) +
    scale_x_continuous(limits = range(plot_df$sr_gene_log2FC, na.rm = TRUE)) +
    scale_y_continuous(limits = range(plot_df$lr_gene_log2FC, na.rm = TRUE)) +
    labs(
      title = title,
      x = expression("srRNA-seq Gene Expression Log"[2] * "FC"),
      y = expression("lrRNA-seq Gene Expression Log"[2] * "FC")
    ) +
    theme_classic(base_size = 7) +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(0.12, "cm"),
      axis.text.x = element_text(size = 5),
      axis.text.y = element_text(size = 5),
      axis.title.x = element_text(size = 6, margin = margin(t = 2)),
      axis.title.y = element_text(size = 6, margin = margin(r = 2)),
      plot.title = element_text(size = 7, face = "bold", hjust = 0.5, margin = margin(b = 4))
    )

  ggsave(
    filename = file.path(figure_dir, output_name),
    plot = p,
    width = 70,
    height = 70,
    units = "mm",
    dpi = 600,
    bg = "white"
  )

  readr::write_csv(
    plot_df,
    file.path(figure_dir, str_replace(output_name, "\\.png$", "_source_data.csv"))
  )
}

# ==============================================================================
# Generate long-read inferred gene-level DESeq2 results
# ==============================================================================

counts <- purrr::map_dfr(classification_paths, function(path) {
  readr::read_tsv(path, show_col_types = FALSE) %>%
    filter(filter_rescued_result == "Isoform") %>%
    select(associated_gene, gene_symbol, gene_exp) %>%
    mutate(
      Sample = get_sample_name(path),
      Condition = get_condition(path)
    )
})

counts_gene_sample <- counts %>%
  mutate(
    gene = case_when(
      !is.na(gene_symbol) & str_trim(gene_symbol) != "" ~ gene_symbol,
      !is.na(associated_gene) & str_trim(associated_gene) != "" ~ associated_gene,
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(gene)) %>%
  group_by(gene, Sample) %>%
  summarise(gene_exp = sum(gene_exp, na.rm = TRUE), .groups = "drop")

count_matrix <- counts_gene_sample %>%
  pivot_wider(
    names_from = Sample,
    values_from = gene_exp,
    values_fill = 0
  )

sample_cols <- colnames(count_matrix)[2:ncol(count_matrix)]

count_data <- count_matrix[, sample_cols] %>%
  mutate(across(everything(), round)) %>%
  as.data.frame()

rownames(count_data) <- count_matrix$gene

sample_info <- data.frame(
  sample = sample_cols,
  condition = sample_cols
) %>%
  mutate(
    time_point = sub("^(C[0-9]+).*", "\\1", condition),
    treatment = case_when(
      grepl("Injured", condition) ~ "Injured",
      grepl("Sciatic", condition) ~ "Control",
      TRUE ~ "Baseline"
    ),
    group = factor(
      paste0(time_point, "_", treatment),
      levels = c("C0_Control", "C3_Injured", "C7_Injured")
    )
  )

rownames(sample_info) <- sample_info$sample

dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = sample_info,
  design = ~ group
)

dds <- DESeq(dds)

result_list <- list(
  C3_vs_C0 = extract_results(dds, "C3_Injured", "C0_Control"),
  C7_vs_C0 = extract_results(dds, "C7_Injured", "C0_Control"),
  C7_vs_C3 = extract_results(dds, "C7_Injured", "C3_Injured")
)

for (nm in names(result_list)) {
  write_csv(
    result_list[[nm]],
    file.path(isoform_results_dir, paste0(nm, "_inferred_gene_results.csv"))
  )
}

# ==============================================================================
# Figure S1a-c: Concordance scatterplots
# ==============================================================================

purrr::iwalk(comparison_config, function(cfg, nm) {
  plot_concordance(
    short_read_file = cfg$short_read_file,
    long_read_file = cfg$long_read_file,
    title = cfg$title,
    point_color = cfg$color,
    output_name = cfg$output
  )
})

message("Supplementary Figure 1 complete.")
