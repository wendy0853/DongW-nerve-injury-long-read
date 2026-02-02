# ---------------------------------------------------------------
# Isoform-Level Differential Expression Analysis (DESeq2)
# Author: Wendy Dong | WUSM
# ---------------------------------------------------------------

# --- Setup -----------------------------------------------------
setwd("/Users/wendydong/Documents/WUSM PhD/Long Read RNA Seq/Analysis_Long_Short_updated/")
data_dir <- "/Users/wendydong/Documents/WUSM PhD/Long Read RNA Seq/Long_read_results_updated/"
results_dir <- file.path(getwd(), "C3_C7_vs_C0_results_FINAL")

# --- Libraries -------------------------------------------------
library(DESeq2)
library(dplyr)
library(readr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)

# --- Load Annotated Count Table --------------------------------
counts <- read_tsv(file.path(results_dir, "all_transcripts_with_associated_genes_count.tsv"), show_col_types = FALSE)
tpm <- read_tsv(file.path(results_dir, "all_transcripts_with_associated_genes_tpm.tsv"), show_col_types = FALSE)

# --- TPM Filtering (≥1 TPM in ≥2 samples) -----------------------
tpm_filtered <- tpm %>%
  mutate(samples_ge1 = rowSums(across(5:ncol(.), ~ .x >= 1))) %>%
  filter(samples_ge1 >= 2)
write_tsv(tpm_filtered, file.path(results_dir, "TPM_filtered.tsv"))

filtered_ids <- tpm_filtered$transcript_symbol
counts <- counts %>% filter(transcript_symbol %in% filtered_ids)

# --- Prepare DESeq2 Inputs -------------------------------------
sample_cols <- colnames(counts)[5:ncol(counts)]
count_data <- counts[, sample_cols] %>% mutate_all(round)
rownames(count_data) <- counts$transcript_id

# Sample metadata
sample_info <- data.frame(
  sample = sample_cols,
  condition = sample_cols
) %>%
  mutate(
    time_point = sub("^(C[0-9]+).*", "\\1", condition),
    treatment = case_when(
      grepl("Injured", condition) ~ "Injured",
      grepl("Control", condition) ~ "Control",
      TRUE ~ "Baseline"
    ),
    group = factor(paste0(time_point, "_", treatment),
                   levels = c("C0_Control", "C3_Injured", "C7_Injured"))
  )

# --- Run DESeq2 ------------------------------------------------
dds <- DESeqDataSetFromMatrix(count_data, sample_info, design = ~ group)
dds <- DESeq(dds)
#saveRDS(dds, file.path(results_dir, "dds_full_object.rds"))
#dds <- readRDS(file.path(results_dir, "dds_full_object.rds"))

# --- Helper Function for Result Extraction ---------------------
extract_results <- function(dds, num, den) {
  res <- results(dds, contrast = c("group", num, den), cooksCutoff = FALSE) %>% as.data.frame()
  norm_counts <- counts(dds, normalized = TRUE)
  
  res %>%
    mutate(
      transcript_id = rownames(res),
      transcript_symbol = counts$transcript_symbol[match(transcript_id, counts$transcript_id)],
      associated_gene = counts$associated_gene[match(transcript_id, counts$transcript_id)],
      gene_symbol = counts$gene_symbol[match(transcript_id, counts$transcript_id)],
      mean_num = rowMeans(norm_counts[, colData(dds)$group == num]),
      mean_den = rowMeans(norm_counts[, colData(dds)$group == den])
    ) %>%
    arrange(padj)
}

# --- Perform Comparisons ---------------------------------------
res_C3_vs_C0 <- extract_results(dds, "C3_Injured", "C0_Control")
res_C7_vs_C0 <- extract_results(dds, "C7_Injured", "C0_Control")
res_C7_vs_C3 <- extract_results(dds, "C7_Injured", "C3_Injured")

# --- Save Results ----------------------------------------------
result_list <- list(
  C3_vs_C0 = res_C3_vs_C0,
  C7_vs_C0 = res_C7_vs_C0,
  C7_vs_C3 = res_C7_vs_C3
)

for (nm in names(result_list)) {
  write_csv(result_list[[nm]], file.path(results_dir, paste0(nm, "_isoform_results.csv")))
}

# --- Visualization ---------------------------------------------
## MA Plot
plot_MA <- function(res, title) {
  ggplot(res, aes(x = log10(baseMean), y = log2FoldChange, color = padj < 0.05)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("grey", "red")) +
    labs(title = title, x = "log10(baseMean)", y = "log2FoldChange") +
    theme_minimal()
}

## Volcano Plot
plot_volcano <- function(res, title) {
  EnhancedVolcano(res,
                  lab = res$transcript_symbol,
                  x = "log2FoldChange", y = "padj",
                  pCutoff = 0.05, FCcutoff = 1,
                  title = title, subtitle = "Differential Transcript Expression",
                  caption = "p_adj < 0.05 & |log2FC| > 1",
                  legendPosition = "bottom")
}

for (nm in names(result_list)) {
  ggsave(paste0("MA_", nm, ".png"), plot_MA(result_list[[nm]], nm), width = 8, height = 6)
  ggsave(paste0("Volcano_", nm, ".png"), plot_volcano(result_list[[nm]], nm), width = 8, height = 6)
}

# --- PCA --------------------------------------------------------
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = c("time_point", "treatment"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = time_point, shape = treatment)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  ggtitle("PCA: Nerve Crush Injury Isoforms")
ggsave(file.path(results_dir,"PCA_nerve_injury.png"), pca_plot, width = 8, height = 6)

# --- Heatmap (Top 50 DE Transcripts per Comparison) -------------
top_tx <- unique(c(
  head(res_C3_vs_C0$transcript_id, 50),
  head(res_C7_vs_C0$transcript_id, 50),
  head(res_C7_vs_C3$transcript_id, 50)
))
mat <- assay(vsd)[top_tx, ] - rowMeans(assay(vsd)[top_tx, ])

ann_col <- data.frame(Time = dds$time_point, Treatment = dds$treatment,
                      row.names = colnames(mat))
pheatmap(mat,
         annotation_col = ann_col,
         main = "Top Differentially Expressed Isoforms",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         clustering_method = "ward.D2",
         angle_col = 315)


