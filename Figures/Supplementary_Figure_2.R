#!/usr/bin/env Rscript

################################################################################
# Supplementary Figure 2
#
# Purpose:
#   Generate supplementary quality-control and differential transcript expression
#   summary plots.
#
# Panels:
#   S2a: PCA plot from variance-stabilized transcript counts
#   S2b: Heatmap of top differentially expressed transcripts
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(pheatmap)
  library(grid)
})

# ==============================================================================
# User-defined directories
# ==============================================================================

isoform_results_dir <- "/path/to/isoform_analysis_results"      # <-- MODIFY HERE
figure_dir <- "/path/to/supplementary_figure_output"            # <-- MODIFY HERE

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Input files
# ==============================================================================

dds_path <- file.path(isoform_results_dir, "dds_full_object.rds")          # <-- MODIFY HERE IF NEEDED
vsd_path <- file.path(isoform_results_dir, "vsd_full_object.rds")          # <-- MODIFY HERE IF NEEDED
res_c3_path <- file.path(isoform_results_dir, "C3_vs_C0_isoform_results.csv")
res_c7_path <- file.path(isoform_results_dir, "C7_vs_C0_isoform_results.csv")
res_c7_c3_path <- file.path(isoform_results_dir, "C7_vs_C3_isoform_results.csv")

required_files <- c(dds_path, vsd_path, res_c3_path, res_c7_path, res_c7_c3_path)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required file(s):\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

# ==============================================================================
# Parameters
# ==============================================================================

padj_cutoff <- 0.05                                             # <-- MODIFY HERE IF NEEDED
lfc_cutoff <- 1                                                 # <-- MODIFY HERE IF NEEDED

condition_colors <- c(
  "C0" = "#0072B2",
  "C3" = "#E69F00",
  "C7" = "#66BD63"
)

# ==============================================================================
# Helper functions
# ==============================================================================

theme_pub <- function(base_size = 7) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.3),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(0.15, "cm"),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 5.5, vjust = 0.8),
      axis.text.y = element_text(size = 5),
      axis.title.x = element_text(size = 6, margin = margin(t = 3)),
      axis.title.y = element_text(size = 6, margin = margin(r = 3)),
      plot.title = element_text(size = 7, face = "bold", hjust = 0.5, margin = margin(b = 4)),
      legend.key.size = unit(3, "mm"),
      legend.text = element_text(size = 5),
      legend.title = element_blank(),
      plot.margin = margin(t = 8, r = 4, b = 4, l = 4, unit = "mm")
    )
}

save_panel <- function(plot, filename, width_mm, height_mm) {
  ggsave(
    filename = file.path(figure_dir, filename),
    plot = plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 600,
    bg = "white"
  )
}

pick_top_tx <- function(res, n = 25) {
  res %>%
    filter(!is.na(padj), padj <= padj_cutoff, abs(log2FoldChange) >= lfc_cutoff) %>%
    mutate(score = abs(log2FoldChange) * -log10(padj)) %>%
    arrange(desc(score)) %>%
    slice_head(n = n) %>%
    pull(transcript_id)
}

# ==============================================================================
# Load data
# ==============================================================================

dds <- readRDS(dds_path)
vsd <- readRDS(vsd_path)

res_list <- list(
  C3_vs_C0 = readr::read_csv(res_c3_path, show_col_types = FALSE),
  C7_vs_C0 = readr::read_csv(res_c7_path, show_col_types = FALSE),
  C7_vs_C3 = readr::read_csv(res_c7_c3_path, show_col_types = FALSE)
)

# ==============================================================================
# Figure S2a: PCA plot
# ==============================================================================

pca_data <- plotPCA(vsd, intgroup = c("time_point", "treatment"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(
  pca_data,
  aes(PC1, PC2, color = time_point, shape = treatment)
) +
  geom_point(size = 3) +
  scale_color_manual(values = condition_colors) +
  labs(
    x = paste0("PC1: ", percent_var[1], "%"),
    y = paste0("PC2: ", percent_var[2], "%"),
    title = "PCA: Nerve Crush Injury Isoforms"
  ) +
  theme_pub()

save_panel(
  pca_plot,
  "FigS2A_PCA_Plot.png",
  width_mm = 70,
  height_mm = 70
)

# ==============================================================================
# Figure S2b: Heatmap of top DETs
# ==============================================================================

top_tx <- unique(c(
  pick_top_tx(res_list$C3_vs_C0),
  pick_top_tx(res_list$C7_vs_C0),
  pick_top_tx(res_list$C7_vs_C3)
))

tx_map <- bind_rows(
  res_list$C3_vs_C0 %>% select(transcript_id, transcript_symbol),
  res_list$C7_vs_C0 %>% select(transcript_id, transcript_symbol),
  res_list$C7_vs_C3 %>% select(transcript_id, transcript_symbol)
) %>%
  distinct(transcript_id, .keep_all = TRUE) %>%
  filter(transcript_id %in% top_tx) %>%
  mutate(
    transcript_symbol = if_else(
      is.na(transcript_symbol) | transcript_symbol == "",
      transcript_id,
      transcript_symbol
    )
  )

heat_mat <- assay(vsd)[tx_map$transcript_id, , drop = FALSE]
heat_mat <- heat_mat - rowMeans(heat_mat)
rownames(heat_mat) <- make.unique(tx_map$transcript_symbol)

ann_col <- data.frame(
  Time = dds$time_point,
  Treatment = dds$treatment,
  row.names = colnames(heat_mat)
)

heat_plot <- pheatmap(
  heat_mat,
  annotation_col = ann_col,
  main = "Top Differentially Expressed Isoforms",
  clustering_method = "ward.D2",
  angle_col = 315,
  silent = TRUE
)

save_panel(
  heat_plot$gtable,
  "FigS2B_Heatmap.png",
  width_mm = 150,
  height_mm = 300
)

message("Supplementary Figure 2 complete.")
