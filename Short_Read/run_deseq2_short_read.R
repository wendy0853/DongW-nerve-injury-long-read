#!/usr/bin/env Rscript

# ============================================================
# Short-read RNA-seq DESeq2 differential gene expression script
# ============================================================
#
# Input:
#   - featureCounts gene-level count matrix (.xlsx or .csv)
#
# Output:
#   - DESeq2 results tables
#   - significant DEG tables
#   - normalized counts
#   - PCA, MA, and volcano plots
#
#
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(readxl)
  library(DESeq2)
  library(tidyverse)
  library(ggplot2)
  library(EnhancedVolcano)
  library(pheatmap)
})

# -----------------------------
# Command-line arguments
# -----------------------------

option_list <- list(
  make_option(c("--counts"), type = "character", help = "Path to gene count matrix (.xlsx or .csv)"),
  make_option(c("--outdir"), type = "character", default = "results", help = "Output directory [default: results]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$counts)) {
  stop("Please provide a count matrix using --counts")
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# User-adjustable settings
# -----------------------------

padj_cutoff <- 0.05
lfc_cutoff <- 1
min_count <- 10
min_replicates <- 2

# If your count matrix has annotation columns before sample columns,
# update this number if needed.
first_sample_column <- 9

# -----------------------------
# Load count matrix
# -----------------------------

message("Reading count matrix: ", opt$counts)

if (grepl("\\.xlsx$", opt$counts, ignore.case = TRUE)) {
  expression_matrix_raw <- readxl::read_excel(opt$counts)
} else if (grepl("\\.csv$", opt$counts, ignore.case = TRUE)) {
  expression_matrix_raw <- read.csv(opt$counts, check.names = FALSE)
} else {
  stop("Input count matrix must be .xlsx or .csv")
}

expression_matrix_raw <- as.data.frame(expression_matrix_raw)

# Clean sample column names if featureCounts/core output used sample.* prefix
colnames(expression_matrix_raw) <- gsub("^sample\\.", "", colnames(expression_matrix_raw))

sample_cols <- colnames(expression_matrix_raw)[first_sample_column:ncol(expression_matrix_raw)]

required_annotation_cols <- c(
  "ensembl_gene_id",
  "external_gene_name",
  "description",
  "gene_biotype"
)

missing_cols <- setdiff(required_annotation_cols, colnames(expression_matrix_raw))
if (length(missing_cols) > 0) {
  stop("Missing required annotation columns: ", paste(missing_cols, collapse = ", "))
}

# -----------------------------
# Remove duplicate gene IDs
# -----------------------------

expression_matrix_raw <- expression_matrix_raw %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

# -----------------------------
# Prepare count matrix
# -----------------------------

count_data <- expression_matrix_raw[, sample_cols]
count_data <- as.data.frame(lapply(count_data, as.integer))
rownames(count_data) <- expression_matrix_raw$ensembl_gene_id

# Remove genes with very low expression
keep <- rowSums(count_data >= min_count) >= min_replicates
count_data <- count_data[keep, ]

gene_annotation <- expression_matrix_raw[keep, c(
  "ensembl_gene_id",
  "external_gene_name",
  "description",
  "gene_biotype"
)]

# -----------------------------
# Create sample metadata
# -----------------------------

sample_info <- data.frame(
  sample = sample_cols,
  condition = sample_cols,
  stringsAsFactors = FALSE
)

sample_info$time_point <- sub("^(C[0-9]+).*", "\\1", sample_info$condition)

sample_info$treatment <- case_when(
  grepl("Injured", sample_info$condition, ignore.case = TRUE) ~ "Injured",
  grepl("Control", sample_info$condition, ignore.case = TRUE) ~ "Control",
  sample_info$time_point == "C0" ~ "Baseline",
  TRUE ~ "Unknown"
)

sample_info$group <- paste(sample_info$time_point, sample_info$treatment, sep = "_")
sample_info$group[sample_info$time_point == "C0"] <- "C0_Baseline"

sample_info$group <- factor(
  sample_info$group,
  levels = c(
    "C0_Baseline",
    "C3_Control",
    "C3_Injured",
    "C7_Control",
    "C7_Injured"
  )
)

rownames(sample_info) <- sample_info$sample

if (any(is.na(sample_info$group))) {
  stop("Some samples could not be assigned to expected groups. Check sample names.")
}

# -----------------------------
# Run DESeq2
# -----------------------------

dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = sample_info,
  design = ~ group
)

mcols(dds)$gene_symbol <- gene_annotation$external_gene_name
mcols(dds)$ensembl_id <- gene_annotation$ensembl_gene_id
mcols(dds)$description <- gene_annotation$description
mcols(dds)$gene_biotype <- gene_annotation$gene_biotype

dds <- DESeq(dds)

# Save normalized counts
normalized_counts <- counts(dds, normalized = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("ensembl_id") %>%
  left_join(
    gene_annotation %>%
      rename(gene_symbol = external_gene_name),
    by = c("ensembl_id" = "ensembl_gene_id")
  ) %>%
  relocate(gene_symbol, .after = ensembl_id)

write.csv(
  normalized_counts,
  file.path(opt$outdir, "normalized_counts.csv"),
  row.names = FALSE
)

# -----------------------------
# Helper functions
# -----------------------------

extract_results <- function(dds, numerator, denominator, gene_annotation) {
  res <- results(dds, contrast = c("group", numerator, denominator))
  res_df <- as.data.frame(res) %>%
    rownames_to_column("ensembl_id") %>%
    left_join(
      gene_annotation %>%
        rename(gene_symbol = external_gene_name),
      by = c("ensembl_id" = "ensembl_gene_id")
    )

  norm_counts <- counts(dds, normalized = TRUE)

  numerator_samples <- rownames(colData(dds))[colData(dds)$group == numerator]
  denominator_samples <- rownames(colData(dds))[colData(dds)$group == denominator]

  res_df$mean_counts_numerator <- rowMeans(norm_counts[res_df$ensembl_id, numerator_samples, drop = FALSE])
  res_df$mean_counts_denominator <- rowMeans(norm_counts[res_df$ensembl_id, denominator_samples, drop = FALSE])

  res_df %>%
    select(
      gene_symbol,
      ensembl_id,
      gene_biotype,
      description,
      baseMean,
      log2FoldChange,
      lfcSE,
      stat,
      pvalue,
      padj,
      mean_counts_numerator,
      mean_counts_denominator
    ) %>%
    arrange(padj)
}

filter_sig_degs <- function(res_df, padj_cutoff = 0.05, lfc_cutoff = 1) {
  res_df %>%
    filter(!is.na(padj)) %>%
    filter(padj <= padj_cutoff, abs(log2FoldChange) >= lfc_cutoff)
}

save_result_pair <- function(res_df, comparison_name) {
  write.csv(
    res_df,
    file.path(opt$outdir, paste0(comparison_name, "_results.csv")),
    row.names = FALSE
  )

  sig_df <- filter_sig_degs(res_df, padj_cutoff, lfc_cutoff)

  write.csv(
    sig_df,
    file.path(opt$outdir, paste0(comparison_name, "_significant_DEGs.csv")),
    row.names = FALSE
  )

  return(sig_df)
}

create_ma_plot <- function(res_df, title) {
  ggplot(res_df, aes(x = log10(baseMean + 1), y = log2FoldChange)) +
    geom_point(aes(color = !is.na(padj) & padj <= padj_cutoff & abs(log2FoldChange) >= lfc_cutoff),
               alpha = 0.6, size = 1) +
    labs(
      title = title,
      x = "log10(baseMean + 1)",
      y = "log2FoldChange",
      color = "Significant"
    ) +
    theme_minimal()
}

create_volcano_plot <- function(res_df, title) {
  EnhancedVolcano(
    res_df,
    lab = res_df$gene_symbol,
    x = "log2FoldChange",
    y = "padj",
    title = title,
    pCutoff = padj_cutoff,
    FCcutoff = lfc_cutoff,
    pointSize = 1.5,
    labSize = 3
  )
}

# -----------------------------
# Define comparisons
# -----------------------------

comparisons <- list(
  C3_Injured_vs_C3_Control = c("C3_Injured", "C3_Control"),
  C7_Injured_vs_C7_Control = c("C7_Injured", "C7_Control"),
  C3_Injured_vs_C0_Baseline = c("C3_Injured", "C0_Baseline"),
  C7_Injured_vs_C0_Baseline = c("C7_Injured", "C0_Baseline"),
  C7_Injured_vs_C3_Injured = c("C7_Injured", "C3_Injured")
)

all_results <- list()
all_sig_results <- list()

for (comparison_name in names(comparisons)) {
  numerator <- comparisons[[comparison_name]][1]
  denominator <- comparisons[[comparison_name]][2]

  message("Running comparison: ", comparison_name)

  res_df <- extract_results(dds, numerator, denominator, gene_annotation)
  sig_df <- save_result_pair(res_df, comparison_name)

  all_results[[comparison_name]] <- res_df
  all_sig_results[[comparison_name]] <- sig_df

  ma_plot <- create_ma_plot(res_df, paste("MA Plot -", comparison_name))
  ggsave(
    file.path(opt$outdir, paste0("MA_plot_", comparison_name, ".png")),
    ma_plot,
    width = 8,
    height = 6,
    dpi = 300
  )

  volcano_plot <- create_volcano_plot(res_df, paste("Volcano Plot -", comparison_name))
  ggsave(
    file.path(opt$outdir, paste0("Volcano_plot_", comparison_name, ".png")),
    volcano_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
}


# -----------------------------
# PCA plot
# -----------------------------

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = c("time_point", "treatment"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = time_point, shape = treatment)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("PCA Plot of Short-Read RNA-seq Samples") +
  theme_minimal()

ggsave(
  file.path(opt$outdir, "PCA_plot_short_read_RNAseq.png"),
  pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# -----------------------------
# Heatmap of top DEGs
# -----------------------------

get_top_genes <- function(res_df, n = 50) {
  res_df %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    slice_head(n = n) %>%
    pull(ensembl_id)
}

top_genes <- unique(c(
  get_top_genes(all_results$C3_Injured_vs_C3_Control),
  get_top_genes(all_results$C7_Injured_vs_C7_Control),
  get_top_genes(all_results$C3_Injured_vs_C0_Baseline),
  get_top_genes(all_results$C7_Injured_vs_C0_Baseline)
))

top_genes <- top_genes[top_genes %in% rownames(vsd)]

mat <- assay(vsd)[top_genes, , drop = FALSE]
mat <- mat - rowMeans(mat)

annotation_col <- data.frame(
  Time_point = colData(dds)$time_point,
  Treatment = colData(dds)$treatment,
  row.names = colnames(mat)
)

png(
  file.path(opt$outdir, "Heatmap_top_DEGs_short_read_RNAseq.png"),
  width = 10,
  height = 8,
  units = "in",
  res = 300
)

pheatmap(
  mat,
  annotation_col = annotation_col,
  show_rownames = FALSE,
  main = "Top Differentially Expressed Genes",
  clustering_method = "ward.D2",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  angle_col = 315
)

dev.off()

message("DESeq2 analysis complete. Results saved to: ", opt$outdir)
