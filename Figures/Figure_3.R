#!/usr/bin/env Rscript

################################################################################
# Figure 3: Isoform-level analysis reveals differential transcript expression 
#  independent of gene-level changes after peripheral nerve injury. 
#
#  Panels:
#   A & D: Short-read differential gene expression volcano plots for C3/C7 vs. C0
#   B & E: Long-read differential transcript expression volcano plots for C3/C7 vs. C0
#   C & F: GO enrichment plots for DTE-only genes
#
#  Panels G - L are generated using:
#     Plotting_isoform_expression.R
#     Plotting_isoform_proportion.R
#     Plotting_isoform_trackplot.R
#
################################################################################

# ==============================================================================
# Setup
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(EnhancedVolcano)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(DOSE)
  library(scales)
  library(grid)
})

# ==============================================================================
# User-defined directories
# ==============================================================================

short_read_dir <- "/path/to/short_read_results"                 # <-- MODIFY HERE
long_read_dir  <- "/path/to/isoform_analysis_results"           # <-- MODIFY HERE
figure_dir     <- "/path/to/Figure_3_output"                    # <-- MODIFY HERE

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# User-adjustable parameters
# ==============================================================================

padj_cutoff <- 0.05                                             # <-- MODIFY HERE IF NEEDED
lfc_cutoff  <- 1                                                # <-- MODIFY HERE IF NEEDED

# Choose comparison: "C3" or "C7"
comparison <- "C7"                                              # <-- MODIFY HERE

# ==============================================================================
# Comparison-specific settings
# ==============================================================================

comparison_config <- list(
  C3 = list(
    short_read_file = "C3_Injured_vs_C0_results.csv",           # <-- MODIFY HERE IF NEEDED
    long_read_file  = "C3_vs_C0_isoform_results.csv",           # <-- MODIFY HERE IF NEEDED
    contrast_label  = "C3 Injured vs. C0 Control",
    color           = "#E69F00",
    dge_labels      = c("Cdkn2c", "Pdlim7", "Cyth2"),          # <-- MODIFY HERE IF NEEDED
    dte_labels      = c(                                        # <-- MODIFY HERE IF NEEDED
      "Cdkn2c-201", "Cdkn2c-202",
      "Pdlim7-202", "Pdlim7-203", "Pdlim7-204", "Pdlim7-209",
      "Cyth2-202", "Cyth2-204"
    ),
    dte_xlim        = c(-30, 30)                                # <-- MODIFY HERE IF NEEDED
  ),

  C7 = list(
    short_read_file = "C7_Injured_vs_C0_results.csv",           # <-- MODIFY HERE IF NEEDED
    long_read_file  = "C7_vs_C0_isoform_results.csv",           # <-- MODIFY HERE IF NEEDED
    contrast_label  = "C7 Injured vs. C0 Control",
    color           = "#66BD63",
    dge_labels      = c("Itgb5", "Sema4c", "Schip1"),          # <-- MODIFY HERE IF NEEDED
    dte_labels      = c(                                        # <-- MODIFY HERE IF NEEDED
      "Itgb5-202", "Itgb5-201", "Itgb5-207",
      "Schip1-203", "Schip1-205",
      "Sema4c-202", "Sema4c-208"
    ),
    dte_xlim        = c(-30, 30)                                # <-- MODIFY HERE IF NEEDED
  )
)

cfg <- comparison_config[[comparison]]

if (is.null(cfg)) {
  stop("comparison must be one of: ", paste(names(comparison_config), collapse = ", "))
}

# ==============================================================================
# Helper functions
# ==============================================================================

check_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }

  invisible(path)
}

theme_pub <- function(base_size = 7) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.3),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(0.15, "cm"),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 5, hjust = 0.5, vjust = 0.5),
      axis.text.y = element_text(size = 5),
      axis.title.x = element_text(size = 6, margin = margin(t = 3)),
      axis.title.y = element_text(size = 6, margin = margin(r = 3)),
      plot.title = element_text(
        size = 7,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 4)
      ),
      legend.key.size = unit(3, "mm"),
      legend.text = element_text(size = 5),
      legend.title = element_blank(),
      plot.margin = margin(t = 8, r = 4, b = 4, l = 4, unit = "mm")
    )
}

theme_volcano <- function(base_size = 7) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.text = element_text(size = 6, color = "black"),
      axis.title = element_text(size = 6, face = "bold"),
      plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 6, hjust = 0.5),
      plot.caption = element_text(size = 5, color = "gray30"),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.3),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(0.15, "cm"),
      legend.position = "none"
    )
}

save_panel <- function(plot, filename, width_mm = 60, height_mm = 60) {
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

filter_significant <- function(res, padj_cut = padj_cutoff, lfc_cut = lfc_cutoff) {
  res %>%
    filter(!is.na(padj)) %>%
    filter(padj < padj_cut, abs(log2FoldChange) > lfc_cut)
}

to_entrez <- function(symbols) {
  symbols <- unique(symbols[!is.na(symbols) & symbols != ""])

  if (length(symbols) == 0) {
    return(character())
  }

  suppressMessages(
    clusterProfiler::bitr(
      symbols,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Mm.eg.db
    )
  ) %>%
    pull(ENTREZID) %>%
    unique()
}

parse_gene_ratio <- function(x) {
  purrr::map_dbl(x, function(ratio) {
    parts <- stringr::str_split(ratio, "/", simplify = TRUE)
    as.numeric(parts[1]) / as.numeric(parts[2])
  })
}

# ==============================================================================
# Load input data
# ==============================================================================

gene_file <- file.path(short_read_dir, cfg$short_read_file)
transcript_file <- file.path(long_read_dir, cfg$long_read_file)

check_file_exists(gene_file)
check_file_exists(transcript_file)

gene_res <- readr::read_csv(gene_file, show_col_types = FALSE)
transcript_res <- readr::read_csv(transcript_file, show_col_types = FALSE)

gene_sig <- filter_significant(gene_res)
transcript_sig <- filter_significant(transcript_res)

# ==============================================================================
# Define DTE-only genes
# ==============================================================================

dge_ids <- unique(gene_sig$gene_symbol)
dte_ids <- unique(transcript_sig$gene_symbol)

dte_only_ids <- setdiff(dte_ids, dge_ids)

dte_only_df <- transcript_sig %>%
  filter(gene_symbol %in% dte_only_ids)

readr::write_csv(
  dte_only_df,
  file.path(figure_dir, paste0("Fig3_", comparison, "_DTE_only_genes_source_data.csv"))
)

# ==============================================================================
# Volcano plot function
# ==============================================================================

make_volcano <- function(
    res,
    label_col,
    selected_labels,
    title,
    subtitle,
    highlight_color,
    xlim = NULL,
    ylim = c(0, 100),
    point_size = 0.8,
    max_overlaps = 20
) {

  plot_df <- res %>%
    mutate(
      padj = replace_na(padj, 1),
      label = if_else(
        .data[[label_col]] %in% selected_labels,
        .data[[label_col]],
        NA_character_
      )
    )

  if (is.null(xlim)) {
    lfc_max <- max(abs(plot_df$log2FoldChange), na.rm = TRUE)
    xlim <- c(
      -ceiling(lfc_max + 0.5),
      ceiling(lfc_max + 0.5)
    )
  }

  EnhancedVolcano(
    plot_df,
    lab = plot_df$label,
    selectLab = selected_labels,
    x = "log2FoldChange",
    y = "padj",
    pCutoff = padj_cutoff,
    FCcutoff = lfc_cutoff,
    xlim = xlim,
    ylim = ylim,
    col = c("gray50", "gray50", "gray50", highlight_color),
    boxedLabels = TRUE,
    max.overlaps = max_overlaps,
    labCol = "black",
    labSize = 1.5,
    drawConnectors = TRUE,
    colConnectors = "black",
    widthConnectors = 0.2,
    pointSize = point_size,
    cutoffLineWidth = 0.2,
    xlab = NULL,
    ylab = NULL,
    title = NULL,
    subtitle = NULL,
    caption = NULL,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    legendPosition = "none"
  ) +
    theme_volcano() +
    labs(
      x = bquote(~Log[2]~italic(FC)),
      y = bquote(~-Log[10]~italic(adj.~P)),
      title = title,
      subtitle = bquote(italic(.(subtitle)))
    )
}

# ==============================================================================
# Figure 3 volcano panels
# ==============================================================================

p_dge <- make_volcano(
  res = gene_res,
  label_col = "gene_symbol",
  selected_labels = cfg$dge_labels,
  title = "Differential Gene Expression",
  subtitle = cfg$contrast_label,
  highlight_color = cfg$color,
  point_size = 0.8,
  max_overlaps = 20
)

save_panel(
  p_dge,
  paste0("Fig3_", comparison, "_DGE.png"),
  width_mm = 60,
  height_mm = 60
)

p_dte <- make_volcano(
  res = transcript_res,
  label_col = "transcript_symbol",
  selected_labels = cfg$dte_labels,
  title = "Differential Transcript Expression",
  subtitle = cfg$contrast_label,
  highlight_color = cfg$color,
  xlim = cfg$dte_xlim,
  point_size = 0.7,
  max_overlaps = Inf
)

save_panel(
  p_dte,
  paste0("Fig3_", comparison, "_DTE.png"),
  width_mm = 67,
  height_mm = 60
)

# ==============================================================================
# GO enrichment for DTE-only genes
# ==============================================================================

run_dte_only_go <- function(
    dte_only_table,
    universe_table,
    top_n = 6,
    simplify_cutoff = 0.8
) {

  dte_genes <- dte_only_table %>%
    distinct(gene_symbol) %>%
    pull(gene_symbol)

  universe_genes <- universe_table %>%
    pull(gene_symbol) %>%
    unique()

  gene_entrez <- to_entrez(dte_genes)
  universe_entrez <- to_entrez(universe_genes)

  if (length(gene_entrez) == 0 || length(universe_entrez) == 0) {
    warning("No Entrez IDs available for GO enrichment.")
    return(tibble())
  }

  ego <- enrichGO(
    gene = gene_entrez,
    universe = universe_entrez,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = padj_cutoff,
    qvalueCutoff = padj_cutoff,
    readable = TRUE
  )

  if (is.null(ego) || nrow(ego@result) == 0) {
    warning("No significant GO terms found.")
    return(tibble())
  }

  ego <- simplify(
    ego,
    cutoff = simplify_cutoff,
    by = "p.adjust",
    select_fun = min
  )

  ego@result %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    as_tibble()
}

plot_go_dot_from_df <- function(
    go_df,
    title,
    dot_color,
    fdr_max = padj_cutoff,
    top_n = 6
) {

  if (nrow(go_df) == 0) {
    warning("GO result table is empty; skipping GO plot.")
    return(NULL)
  }

  df <- go_df %>%
    filter(p.adjust <= fdr_max) %>%
    mutate(
      gene_ratio_num = parse_gene_ratio(GeneRatio),
      neglog10_fdr = -log10(p.adjust + 1e-300)
    ) %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    arrange(gene_ratio_num) %>%
    mutate(
      Description = factor(Description, levels = Description)
    )

  if (nrow(df) == 0) {
    warning("No GO terms passed plotting threshold; skipping GO plot.")
    return(NULL)
  }

  ggplot(df, aes(x = gene_ratio_num, y = Description)) +
    geom_point(aes(size = Count, color = neglog10_fdr)) +
    scale_y_discrete(labels = label_wrap(20)) +
    scale_size(
      range = c(0.8, 3),
      breaks = pretty_breaks(n = 4)
    ) +
    scale_color_gradient(
      low = "grey30",
      high = dot_color,
      breaks = pretty_breaks(n = 4)
    ) +
    labs(
      title = title,
      x = "Gene ratio",
      y = NULL,
      color = expression(-log[10]("adj. P")),
      size = "\nCount"
    ) +
    theme_pub() +
    theme(
      panel.grid.major = element_line(linewidth = 0.15, color = "grey90"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 5),
      axis.text.y = element_text(size = 5),
      axis.title = element_text(size = 5),
      legend.title = element_text(size = 4),
      legend.text = element_text(size = 4),
      legend.key.height = unit(0.35, "cm"),
      legend.key.width = unit(0.2, "cm"),
      legend.spacing.y = unit(0.1, "cm")
    ) +
    guides(
      color = guide_colorbar(
        order = 1,
        barheight = unit(1.5, "cm"),
        barwidth = unit(0.25, "cm")
      ),
      size = guide_legend(
        order = 2,
        keyheight = unit(0.35, "cm"),
        keywidth = unit(0.35, "cm")
      )
    )
}

go_df <- run_dte_only_go(
  dte_only_table = dte_only_df,
  universe_table = transcript_res,
  top_n = 6
)

readr::write_csv(
  go_df,
  file.path(figure_dir, paste0("Fig3_", comparison, "_DTE_only_GO_source_data.csv"))
)

p_go <- plot_go_dot_from_df(
  go_df,
  title = paste0("GO: ", comparison, " DTE-Only Genes"),
  dot_color = cfg$color,
  top_n = 6
)

if (!is.null(p_go)) {
  save_panel(
    p_go,
    paste0("Fig3_", comparison, "_DTE_GO.png"),
    width_mm = 62,
    height_mm = 70
  )
}

message("Figure 3 plotting complete. Outputs saved to: ", figure_dir)
