#!/usr/bin/env Rscript

################################################################################
# Supplementary Figure 3
#
# Purpose:
#   Generate supplementary DGE/DTE overlap outputs, volcano plots and GO
#   enrichment of DET-only genes.
#
# Panels:
#   S3A: DGE/DTE Venn diagrams
#   S3B-C: DGE and DTE volcano plots for C3/C7 vs. C0
#   S3D: GO enrichment of DET-only genes (background: multi-isoform genes)
#
# Inputs:
#   - all_transcripts_with_associated_genes_count.tsv
#       from Isoform_Analysis/Generating_Isoform_Count_Matrix.R
#   - *_isoform_results.csv
#       from Isoform_Analysis/DTE_Analysis_DESeq2.R
#   - *_inferred_gene_results.csv
#       from Gene_Analysis/DGE_Analysis_Long_Read.R
#
################################################################################

suppressPackageStartupMessages({
  library(EnhancedVolcano)
  library(VennDiagram)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(tidyverse)
  library(scales)
  library(grid)
})

# ==============================================================================
# User-defined files/directories
# ==============================================================================

results_dir <- "/path/to/isoform_analysis_results"               # <-- MODIFY HERE

# Output directory for DGE-only, shared and DTE-only gene tables
overlap_dir <- "/path/to/gene_transcript_overlap_results"        # <-- MODIFY HERE

figure_dir <- "/path/to/supplementary_figure_output"             # <-- MODIFY HERE

dir.create(overlap_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

count_path <- file.path(
  results_dir,
  "all_transcripts_with_associated_genes_count.tsv"
)                                                                # <-- MODIFY HERE IF NEEDED

# ==============================================================================
# Parameters
# ==============================================================================

padj_cutoff <- 0.05                                              # <-- MODIFY HERE IF NEEDED
lfc_cutoff <- 1                                                  # <-- MODIFY HERE IF NEEDED

go_ontology <- "BP"                                              # <-- MODIFY HERE IF NEEDED
go_top_n <- 6                                                    # <-- MODIFY HERE IF NEEDED
go_simplify_cutoff <- 0.6                                        # <-- MODIFY HERE IF NEEDED

# Comparisons to run: c("C3", "C7") for both, or a single one such as "C7"
comparisons <- c("C3", "C7")                                     # <-- MODIFY HERE IF NEEDED

# ==============================================================================
# Comparison-specific settings
# ==============================================================================

comparison_config <- list(
  C3 = list(
    gene_file      = "C3_vs_C0_inferred_gene_results.csv",       # <-- MODIFY HERE IF NEEDED
    isoform_file   = "C3_vs_C0_isoform_results.csv",             # <-- MODIFY HERE IF NEEDED
    contrast_label = "C3 Injured vs. C0 Control",
    color          = "#E69F00",
    gene_labels    = c(                                          # <-- MODIFY HERE IF NEEDED
      "Arg1", "Rrm2", "Cd9", "Cntf", "Top2a", "Chil3"
    ),
    transcript_labels = c(                                       # <-- MODIFY HERE IF NEEDED
      "Arg1-201", "Rrm2-201", "Cd9-201", "Cntf-201", "Top2a-201", "Chil3-201"
    )
  ),

  C7 = list(
    gene_file      = "C7_vs_C0_inferred_gene_results.csv",       # <-- MODIFY HERE IF NEEDED
    isoform_file   = "C7_vs_C0_isoform_results.csv",             # <-- MODIFY HERE IF NEEDED
    contrast_label = "C7 Injured vs. C0 Control",
    color          = "#66BD63",
    gene_labels    = c(                                          # <-- MODIFY HERE IF NEEDED
      "H2ac11", "Art3", "Slc25a1", "Cntf", "Top2a", "Trem2"
    ),
    transcript_labels = c(                                       # <-- MODIFY HERE IF NEEDED
      "H2ac11-201", "Art3-201", "Slc25a1-201", "Cntf-201", "Top2a-201",
      "Trem2-201", "Trem2-202", "Trem2-203"
    )
  )
)

# ==============================================================================
# Input checks
# ==============================================================================

required_files <- c(
  count_path,
  unlist(purrr::map(comparison_config, ~ c(
    file.path(results_dir, .x$gene_file),
    file.path(results_dir, .x$isoform_file)
  )))
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required file(s):\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

# ==============================================================================
# Helper functions
# ==============================================================================

theme_volcano <- function(base_size = 7) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.text = element_text(size = 7, color = "black"),
      axis.title = element_text(size = 7, face = "bold"),
      plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 7, hjust = 0.5),
      plot.caption = element_text(size = 5, color = "gray30"),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.3),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(0.15, "cm"),
      legend.position = "none"
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

filter_significant <- function(res) {
  res %>%
    filter(!is.na(padj), padj <= padj_cutoff, abs(log2FoldChange) >= lfc_cutoff)
}

save_venn <- function(dge_ids, dte_ids, title, filename) {
  venn_plot <- venn.diagram(
    x = list(DGE = dge_ids, DTE = dte_ids),
    filename = NULL,
    fill = c("skyblue", "lightyellow"),
    alpha = 0.5,
    cex = 1.5,
    cat.cex = 1.5,
    main = title,
    main.cex = 2
  )

  png(
    filename = file.path(figure_dir, filename),
    width = 70,
    height = 70,
    units = "mm",
    res = 600,
    bg = "white"
  )
  grid.newpage()
  grid.draw(venn_plot)
  dev.off()
}

make_volcano <- function(res, label_col, selected_labels, title, subtitle,
                         highlight_color, point_size, max_overlaps) {
  plot_df <- res %>%
    mutate(
      padj = replace_na(padj, 1),
      label = if_else(
        .data[[label_col]] %in% selected_labels,
        .data[[label_col]],
        NA_character_
      )
    )

  EnhancedVolcano(
    plot_df,
    lab = plot_df$label,
    selectLab = selected_labels,
    x = "log2FoldChange",
    y = "padj",
    pCutoff = padj_cutoff,
    FCcutoff = lfc_cutoff,
    col = c("gray50", "gray50", "gray50", highlight_color),
    boxedLabels = TRUE,
    max.overlaps = max_overlaps,
    labCol = "black",
    labSize = 1.7,
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

to_entrez <- function(symbols) {
  symbols <- unique(symbols[!is.na(symbols) & symbols != ""])

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
    parts <- str_split(ratio, "/", simplify = TRUE)
    as.numeric(parts[1]) / as.numeric(parts[2])
  })
}

# Assigns an Ensembl gene ID, via the gene symbol, to isoforms that lack one.
backfill_associated_gene <- function(res, gene_id_map) {
  res %>%
    left_join(gene_id_map, by = "gene_symbol", relationship = "many-to-many") %>%
    mutate(
      associated_gene = case_when(
        !is.na(associated_gene) & str_starts(associated_gene, "ENSMUSG") ~ associated_gene,
        !is.na(mapped_associated_gene) ~ mapped_associated_gene,
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(-mapped_associated_gene)
}

run_det_only_go <- function(target_genes, universe_genes) {
  ego <- enrichGO(
    gene = to_entrez(target_genes),
    universe = to_entrez(universe_genes),
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = go_ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = padj_cutoff,
    qvalueCutoff = padj_cutoff,
    readable = TRUE
  )

  clusterProfiler::simplify(
    ego,
    cutoff = go_simplify_cutoff,
    by = "p.adjust",
    select_fun = min
  )@result %>%
    arrange(p.adjust) %>%
    slice_head(n = go_top_n) %>%
    as_tibble()
}

plot_go_dot <- function(go_df, title, dot_color) {
  df <- go_df %>%
    mutate(
      gene_ratio_num = parse_gene_ratio(GeneRatio),
      neglog10_fdr = -log10(p.adjust),
      Description = fct_reorder(str_wrap(Description, width = 38), Count)
    )

  ggplot(df, aes(x = gene_ratio_num, y = Description)) +
    geom_point(aes(size = Count, color = neglog10_fdr)) +
    scale_size(range = c(0.8, 3), breaks = pretty_breaks(n = 4)) +
    scale_color_gradient(low = "grey70", high = dot_color, breaks = pretty_breaks(n = 4)) +
    labs(
      x = "Gene ratio",
      y = NULL,
      title = title,
      color = expression(-Log[10]~italic(adj.~P)),
      size = "Count"
    ) +
    theme_bw(base_size = 8) +
    theme(
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 8, face = "bold"),
      axis.text.x = element_text(size = 7, color = "black"),
      axis.text.y = element_text(size = 8, color = "black"),
      panel.grid.major = element_line(linewidth = 0.15, color = "grey90"),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 6),
      plot.margin = margin(t = 2, r = 4, b = 1, l = 0, unit = "mm")
    )
}

# ==============================================================================
# Gene symbol to Ensembl gene ID map
# ==============================================================================

gene_id_map <- read_tsv(count_path, show_col_types = FALSE) %>%
  filter(
    !is.na(gene_symbol),
    !is.na(associated_gene),
    str_starts(associated_gene, "ENSMUSG")
  ) %>%
  distinct(gene_symbol, mapped_associated_gene = associated_gene)

# ==============================================================================
# Run each comparison
# ==============================================================================

for (comparison in comparisons) {
  cfg <- comparison_config[[comparison]]
  message("---- ", cfg$contrast_label, " ----")

  gene_res <- read_csv(file.path(results_dir, cfg$gene_file), show_col_types = FALSE)
  transcript_raw <- read_csv(file.path(results_dir, cfg$isoform_file), show_col_types = FALSE)

  transcript_res <- backfill_associated_gene(transcript_raw, gene_id_map)

  gene_sig <- filter_significant(gene_res)
  transcript_sig <- filter_significant(transcript_res)

  # ----------------------------------------------------------------------------
  # DGE/DTE overlap (by Ensembl gene ID)
  # ----------------------------------------------------------------------------

  dge_ids <- unique(gene_sig$associated_gene)
  dte_ids <- unique(transcript_sig$associated_gene)

  shared_ids <- intersect(dge_ids, dte_ids)
  dge_only_ids <- setdiff(dge_ids, dte_ids)
  dte_only_ids <- setdiff(dte_ids, dge_ids)

  message(
    "DGE-only: ", length(dge_only_ids),
    " | shared: ", length(shared_ids),
    " | DTE-only: ", length(dte_only_ids)
  )

  dge_only_df <- gene_sig %>%
    filter(associated_gene %in% dge_only_ids)

  shared_df <- left_join(
    transcript_sig %>% filter(associated_gene %in% shared_ids),
    gene_sig %>% filter(associated_gene %in% shared_ids),
    by = "associated_gene",
    suffix = c("_transcript", "_gene")
  )

  dte_only_df <- transcript_sig %>%
    filter(associated_gene %in% dte_only_ids)

  isoform_summary <- transcript_sig %>%
    group_by(gene_symbol) %>%
    summarise(
      num_isoforms = dplyr::n(),
      transcript_symbols = paste(unique(transcript_symbol), collapse = "; "),
      .groups = "drop"
    ) %>%
    arrange(dplyr::desc(num_isoforms))

  write_csv(dge_only_df, file.path(overlap_dir, paste0(comparison, "_Injured_vs_C0_DGE_only_genes.csv")))
  write_csv(shared_df, file.path(overlap_dir, paste0(comparison, "_Injured_vs_C0_DGE_shared_genes.csv")))
  write_csv(dte_only_df, file.path(overlap_dir, paste0(comparison, "_Injured_vs_C0_DTE_only_genes.csv")))
  write_csv(isoform_summary, file.path(overlap_dir, paste0(comparison, "_DTE_isoform_number_shared_summary.csv")))

  # ----------------------------------------------------------------------------
  # Figure S3A: DGE/DTE Venn diagram
  # ----------------------------------------------------------------------------

  save_venn(
    dge_ids = dge_ids,
    dte_ids = dte_ids,
    title = cfg$contrast_label,
    filename = paste0("FigS3A_", comparison, "_DGE_DTE_Venn.png")
  )

  # ----------------------------------------------------------------------------
  # Figure S3B-C: DGE and DTE volcano plots
  # ----------------------------------------------------------------------------

  p_dge <- make_volcano(
    res = gene_res,
    label_col = "gene_symbol",
    selected_labels = cfg$gene_labels,
    title = "Differential Gene Expression",
    subtitle = cfg$contrast_label,
    highlight_color = cfg$color,
    point_size = 0.8,
    max_overlaps = 20
  )

  save_panel(p_dge, paste0("FigS3B_", comparison, "_DGE.png"), width_mm = 75, height_mm = 68)

  p_dte <- make_volcano(
    res = transcript_raw,
    label_col = "transcript_symbol",
    selected_labels = cfg$transcript_labels,
    title = "Differential Transcript Expression",
    subtitle = cfg$contrast_label,
    highlight_color = cfg$color,
    point_size = 0.7,
    max_overlaps = Inf
  )

  save_panel(p_dte, paste0("FigS3C_", comparison, "_DTE.png"), width_mm = 75, height_mm = 68)

  # ----------------------------------------------------------------------------
  # Figure S3D: GO enrichment of DET-only genes (background: multi-isoform genes)
  # ----------------------------------------------------------------------------

  # Single-isoform genes are removed from the DET-only list and the background.
  mono_isoform_genes <- transcript_res %>%
    filter(!is.na(associated_gene), !is.na(transcript_id)) %>%
    distinct(associated_gene, transcript_id) %>%
    dplyr::count(associated_gene, name = "n_isoforms") %>%
    filter(n_isoforms == 1) %>%
    pull(associated_gene)

  det_only_genes <- dte_only_df %>%
    filter(!is.na(associated_gene), !associated_gene %in% mono_isoform_genes) %>%
    pull(gene_symbol) %>%
    unique()

  universe_genes <- transcript_res %>%
    filter(!is.na(associated_gene), !associated_gene %in% mono_isoform_genes) %>%
    pull(gene_symbol) %>%
    unique()

  message(
    "DET-only genes: ", n_distinct(dte_only_df$associated_gene),
    " (", n_distinct(dte_only_df$associated_gene) -
      n_distinct(dte_only_df$associated_gene[!dte_only_df$associated_gene %in% mono_isoform_genes]),
    " mono-isoform removed) | GO background: ", length(universe_genes)
  )

  go_df <- run_det_only_go(det_only_genes, universe_genes)

  write_csv(
    go_df,
    file.path(figure_dir, paste0("FigS3D_", comparison, "_DET_only_GO_source_data.csv"))
  )

  p_go <- plot_go_dot(
    go_df,
    title = paste0("GO: ", comparison, " DET-only Genes"),
    dot_color = cfg$color
  )

  save_panel(p_go, paste0("FigS3D_", comparison, "_DET_only_GO.png"), width_mm = 120, height_mm = 80)
}

message("Supplementary Figure 3 complete.")
