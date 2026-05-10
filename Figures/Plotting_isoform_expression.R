#!/usr/bin/env Rscript

# =============================================================================
# Isoform-level expression visualization across injury conditions
# =============================================================================
#
# Purpose:
#   Plot DESeq2 normalized isoform expression for a gene of interest across C0, C3,
#   and C7 injury conditions.
#
# Inputs:
#   - dds_full_object.rds from Isoform_Analysis/DTE_Analysis_DESeq2.R
#   - Counts_filtered.tsv from Isoform_Analysis/DTE_Analysis_DESeq2.R
#
# Output:
#   - Isoform expression plot for the selected gene (log10 axis)
#
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(RColorBrewer)
  library(data.table)
  library(readr)
  library(grid)
})

data.table::setDTthreads()

# -----------------------------
# User-defined files/directories
# -----------------------------

results_dir <- "/path/to/isoform_analysis_results"             # <-- MODIFY HERE
figure_dir <- "/path/to/figure_output_directory"               # <-- MODIFY HERE

dds_path <- file.path(results_dir, "dds_full_object.rds")       # <-- MODIFY HERE IF NEEDED
counts_path <- file.path(results_dir, "Counts_filtered.tsv")    # <-- MODIFY HERE IF NEEDED

gene_of_interest <- "Nrxn3"                                    # <-- MODIFY HERE

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Input checks
# -----------------------------

required_files <- c(dds_path, counts_path)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required file(s):\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

# -----------------------------
# Load processed data
# -----------------------------

dds <- readRDS(dds_path)

count_matrix <- readr::read_tsv(
  counts_path,
  show_col_types = FALSE
)

expression_matrix_table <- data.table::data.table(count_matrix)

# -----------------------------
# Extract normalized counts
# -----------------------------

norm_counts <- counts(dds, normalized = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("transcript_id")

sample_info <- as.data.frame(colData(dds)) %>%
  rownames_to_column("sample_id")

annotated_counts <- norm_counts %>%
  left_join(
    expression_matrix_table %>%
      dplyr::select(transcript_id, transcript_symbol, gene_symbol),
    by = "transcript_id"
  )

long_counts <- annotated_counts %>%
  pivot_longer(
    cols = all_of(sample_info$sample_id),
    names_to = "sample_id",
    values_to = "normalized_count"
  ) %>%
  left_join(sample_info, by = "sample_id")

# -----------------------------
# Filter to gene of interest
# -----------------------------

plot_data <- long_counts %>%
  filter(gene_symbol == gene_of_interest)

if (nrow(plot_data) == 0) {
  stop("No isoforms found for gene: ", gene_of_interest)
}

# Keep isoforms expressed in at least one condition
isoforms_to_keep <- plot_data %>%
  group_by(transcript_symbol, group) %>%
  summarise(
    mean_expr = mean(normalized_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(transcript_symbol) %>%
  summarise(
    n_groups_expressed = sum(mean_expr > 0),
    .groups = "drop"
  ) %>%
  filter(n_groups_expressed >= 1) %>%
  pull(transcript_symbol)

plot_data <- plot_data %>%
  filter(transcript_symbol %in% isoforms_to_keep) %>%
  mutate(
    norm_count_adj = normalized_count + 1
  )

# -----------------------------
# Summarize expression
# -----------------------------

sum_tbl <- plot_data %>%
  group_by(group, transcript_symbol) %>%
  summarise(
    mean_norm = mean(norm_count_adj, na.rm = TRUE),
    se_norm = sd(norm_count_adj, na.rm = TRUE) / sqrt(sum(!is.na(norm_count_adj))),
    .groups = "drop"
  )

# -----------------------------
# Helper: extended color palette
# -----------------------------

make_extended_palette <- function(n) {
  n <- as.integer(n)
  base <- RColorBrewer::brewer.pal(8, "Set2")

  if (n <= 8) {
    return(base[1:n])
  }

  extra_needed <- n - 8

  c(
    base,
    grDevices::hcl.colors(extra_needed, palette = "Harmonic")
  )
}

pal <- make_extended_palette(length(unique(sum_tbl$transcript_symbol)))

# -----------------------------
# Plot isoform expression
# -----------------------------

p_isoform_log10 <- ggplot(
  sum_tbl,
  aes(
    x = group,
    y = mean_norm,
    color = transcript_symbol,
    group = transcript_symbol
  )
) +
  geom_errorbar(
    aes(
      ymin = pmax(mean_norm - se_norm, 1e-6),
      ymax = mean_norm + se_norm
    ),
    width = 0.15,
    linewidth = 0.4,
    position = position_dodge(0.1)
  ) +
  geom_line(
    linewidth = 0.6,
    position = position_dodge(0.1)
  ) +
  geom_point(
    size = 1.6,
    shape = 16,
    stroke = 0.25,
    position = position_dodge(0.1)
  ) +
  scale_x_discrete(
    labels = c(
      "C0_Control" = "C0",
      "C3_Injured" = "C3",
      "C7_Injured" = "C7"
    ),
    expand = c(0, 0.2)
  ) +
  scale_y_log10(
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_color_manual(values = pal) +
  labs(
    title = gene_of_interest,
    x = NULL,
    y = "Normalized Count (Mean ± SEM)"
  ) +
  theme_minimal(base_size = 7) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.ticks.length = unit(0.12, "cm"),
    axis.text.x = element_text(size = 5, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 5),
    axis.title.y = element_text(size = 6, margin = margin(r = 2)),
    plot.title = element_text(
      size = 7,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 4)
    ),
    legend.position = "right",
    legend.background = element_rect(
      color = "black",
      fill = "white",
      linewidth = 0.25
    ),
    legend.key.size = unit(3, "mm"),
    legend.text = element_text(size = 5),
    legend.title = element_blank(),
    plot.margin = margin(t = 6, r = 3, b = 4, l = 1, unit = "mm")
  )

# -----------------------------
# Save figure
# -----------------------------

out_file <- file.path(
  figure_dir,
  paste0(gene_of_interest, "_IsoformExpression.png")
)

ggsave(
  out_file,
  p_isoform_log10,
  width = 55,
  height = 55,
  units = "mm",
  dpi = 600
)

message("Saved isoform expression plot: ", out_file)
