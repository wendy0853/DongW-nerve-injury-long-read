#!/usr/bin/env Rscript

################################################################################
# Supplementary Figure 3
#
#
# Panels:
#   S3a: Trem2 isoform expression
#   S3b: Trem2 isoform track plot
#   S3c-d: DTE volcano plots
#   S3e: Isoform expression plots for selected injury genes
#
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
})

# ==============================================================================
# User-defined directories
# ==============================================================================

figure_scripts_dir <- "/path/to/Figures"                         # <-- MODIFY HERE
figure_dir <- "/path/to/supplementary_figure_output"             # <-- MODIFY HERE

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Inputs required by shared plotting scripts
# ==============================================================================

# The following scripts are expected to be edited/run with gene-specific inputs:
#
#   Plotting_isoform_expression.R
#   Plotting_isoform_trackplot.R
#   Figure_3.R
#
# This file records which shared scripts were used for each Supplementary Figure 3 panel.

panel_sources <- tibble::tribble(
  ~panel, ~description, ~script, ~gene_or_input,
  "S3A", "Trem2 isoform expression", "Plotting_isoform_expression.R", "Trem2",
  "S3B", "Trem2 isoform structure", "Plotting_isoform_trackplot.R", "Trem2 isoforms",
  "S3C-D", "DTE volcano plots", "Figure_3.R", "C3/C7 DTE volcano panels",
  "S3E", "Isoform expression of injury-associated genes", "Plotting_isoform_expression.R", "Pmp22, Spp1, Atf3, Lgals3"
)

