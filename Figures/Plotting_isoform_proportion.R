#!/usr/bin/env Rscript

# =============================================================================
# Isoform-level proportion visualization across injury conditions
# =============================================================================
#
# Purpose:
#   Plot relative isoform usage/proportion for a gene of interest across C0, C3,
#   and C7 injury conditions.
#
# Input:
#   - ALL_isoform_ratio_per_gene_condition.csv from
#     Isoform_Analysis/DTU_Analysis_IsoformSwitchAnalyzeR_DEXSeq.R
#
# Output:
#   - Stacked bar plot of mean isoform fraction by condition
#
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(RColorBrewer)
  library(data.table)
  library(readr)
  library(scales)
  library(grid)
})

data.table::setDTthreads()

# -----------------------------
# User-defined files/directories
# -----------------------------

isoform_ratio_path <- "/path/to/ALL_isoform_ratio_per_gene_condition.csv" # <-- MODIFY HERE
figure_dir <- "/path/to/figure_output_directory"                         # <-- MODIFY HERE

gene_of_interest <- "Vcan"                                               # <-- MODIFY HERE

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Input checks
# -----------------------------

if (!file.exists(isoform_ratio_path)) {
  stop("Isoform ratio file not found: ", isoform_ratio_path)
}

# -----------------------------
# Load processed data
# -----------------------------

expr_long <- readr::read_csv(
  isoform_ratio_path,
  show_col_types = FALSE
)

required_cols <- c(
  "sample",
  "gene_id",
  "transcript_symbol",
  "condition",
  "tpm",
  "isoform_fraction"
)

missing_cols <- setdiff(required_cols, colnames(expr_long))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns in isoform ratio file: ",
    paste(missing_cols, collapse = ", ")
  )
}

# -----------------------------
# Filter to gene of interest
# -----------------------------

plot_data <- expr_long %>%
  filter(gene_id == gene_of_interest) %>%
  group_by(sample, gene_id) %>%
  mutate(
    total_tpm = sum(tpm, na.rm = TRUE)
  ) %>%
  filter(total_tpm > 1) %>%
  mutate(
    isoform_fraction = replace_na(isoform_fraction, 0)
  ) %>%
  ungroup() %>%
  arrange(transcript_symbol)

if (nrow(plot_data) == 0) {
  stop("No isoform proportion data found for gene: ", gene_of_interest)
}

# -----------------------------
# Color palette
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

isoform_levels <- unique(plot_data$transcript_symbol)

color_palette <- make_extended_palette(length(isoform_levels))
names(color_palette) <- isoform_levels

# -----------------------------
# Plot isoform proportions
# -----------------------------

p <- ggplot(
  plot_data,
  aes(
    x = condition,
    y = isoform_fraction,
    fill = transcript_symbol
  )
) +
  stat_summary(
    fun = mean,
    geom = "bar",
    position = "stack",
    color = "black",
    linewidth = 0.15
  ) +
  scale_fill_manual(values = color_palette) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 10L),
    expand = c(0, 0)
  ) +
  scale_x_discrete(
    expand = expansion(mult = c(0.3, 0.3))
  ) +
  labs(
    title = gene_of_interest,
    x = NULL,
    y = "Proportion of Isoforms",
    fill = "Transcript"
  ) +
  theme_minimal(base_size = 7) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.3),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.ticks.length = unit(0.15, "cm"),
    axis.text.x = element_text(size = 5, color = "black"),
    axis.text.y = element_text(size = 5, color = "black"),
    axis.title.y = element_text(size = 6, margin = margin(r = 2)),
    plot.title = element_text(
      size = 7,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 4)
    ),
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 5),
    legend.key.size = unit(2.5, "mm"),
    legend.spacing.x = unit(1.5, "mm"),
    legend.margin = margin(t = 0, b = 0, unit = "mm"),
    panel.grid = element_blank(),
    plot.margin = margin(t = 6, r = 6, b = 4, l = 4, unit = "mm")
  )

# -----------------------------
# Save figure
# -----------------------------

out_file <- file.path(
  figure_dir,
  paste0(gene_of_interest, "_IsoformProportion.png")
)

ggsave(
  out_file,
  p,
  width = 57,
  height = 55,
  units = "mm",
  dpi = 600
)

message("Saved isoform proportion plot: ", out_file)
