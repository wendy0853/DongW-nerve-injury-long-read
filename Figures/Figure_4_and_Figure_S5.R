#!/usr/bin/env Rscript

# =============================================================================
# Figure 4 and Supplementary Figure S5 single-cell plotting
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tidyr)
  library(scales)
  library(RColorBrewer)
  library(grid)
})

# -----------------------------
# User-defined files/directories
# -----------------------------

results_dir <- "/path/to/single_cell_results"                 # <-- MODIFY HERE
figure_dir <- "/path/to/figure_output_directory"              # <-- MODIFY HERE

seurat_object_path <- file.path(
  results_dir,
  "integrated_annotated_seurat_object.rds"
)                                                             # <-- MODIFY HERE IF NEEDED

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Load object
# -----------------------------

obj <- readRDS(seurat_object_path)

# -----------------------------
# UMAP broad cell types
# -----------------------------

p_umap <- DimPlot(
  obj,
  reduction = "umap.harmony",
  group.by = "annotations",
  label = TRUE,
  label.size = 2.2,
  pt.size = 1e-10,
  raster = FALSE
) +
  ggtitle("UMAP of Injured Sciatic Nerves") +
  theme_void() +
  theme(
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

ggsave(
  file.path(figure_dir, "UMAP_broad_annotations.png"),
  p_umap,
  width = 60,
  height = 60,
  units = "mm",
  dpi = 600
)

# -----------------------------
# UMAP detailed cell types
# -----------------------------

p_umap_detailed <- DimPlot(
  obj,
  reduction = "umap.harmony",
  group.by = "detailed_annotations",
  label = TRUE,
  label.size = 2,
  pt.size = 1e-10,
  raster = FALSE
) +
  ggtitle("UMAP of Injured Sciatic Nerves") +
  theme_void() +
  theme(
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

ggsave(
  file.path(figure_dir, "UMAP_detailed_annotations.png"),
  p_umap_detailed,
  width = 90,
  height = 90,
  units = "mm",
  dpi = 600
)

# -----------------------------
# UMAP by timepoint
# -----------------------------

p_umap_timepoint <- DimPlot(
  obj,
  reduction = "umap.harmony",
  split.by = "timepoint",
  pt.size = 1e-10,
  raster = FALSE
) +
  ggtitle(NULL) +
  theme_void()

ggsave(
  file.path(figure_dir, "UMAP_split_by_timepoint.png"),
  p_umap_timepoint,
  width = 180,
  height = 60,
  units = "mm",
  dpi = 600
)

# -----------------------------
# Highlight each timepoint
# -----------------------------

tp_cols <- c(
  C0 = "#0072B2",
  C3 = "#E69F00",
  C7 = "#66BD63"
)

for (tp in names(tp_cols)) {
  obj$highlight <- factor(
    obj$timepoint == tp,
    levels = c(FALSE, TRUE)
  )

  p <- DimPlot(
    obj,
    reduction = "umap.harmony",
    group.by = "highlight",
    cols = c("FALSE" = "grey85", "TRUE" = tp_cols[[tp]]),
    pt.size = 1e-10,
    order = TRUE,
    raster = FALSE
  ) +
    ggtitle(tp) +
    theme_void() +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      legend.position = "none"
    )

  ggsave(
    file.path(figure_dir, paste0("UMAP_highlight_", tp, ".png")),
    p,
    width = 90,
    height = 90,
    units = "mm",
    dpi = 600
  )
}

# -----------------------------
# Cell composition bar plot
# -----------------------------

pt <- table(obj$annotations, obj$timepoint) %>%
  as.data.frame()

colnames(pt) <- c("cell_type", "timepoint", "Freq")

pt$timepoint <- factor(pt$timepoint, levels = c("C0", "C3", "C7"))

cluster_names <- c(
  "Endoneurial",
  "Epineurial",
  "Perineurial",
  "Schwann Cells",
  "Immune",
  "Pericytes/VSMCs",
  "Endothelial"
)

p_bar <- ggplot(
  pt,
  aes(x = timepoint, y = Freq, fill = factor(cell_type, levels = cluster_names))
) +
  geom_col(position = "fill", width = 0.85, color = "black", linewidth = 0.3) +
  xlab(NULL) +
  ylab("Proportion") +
  ggtitle("Cell Composition") +
  theme_minimal(base_size = 8) +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 9, face = "bold"),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  )

ggsave(
  file.path(figure_dir, "Bargraph_cell_composition.png"),
  p_bar,
  width = 80,
  height = 70,
  units = "mm",
  dpi = 600
)

# -----------------------------
# Canonical marker dot plot
# -----------------------------

canonical_markers <- c(
  "Pcolce2", "Cdkn2a", "Cldn1", "Pecam1", "Notch3",
  "Sox10", "Ccl8", "Cx3cr1", "Napsa"
)

p_dot <- DotPlot(
  obj,
  features = canonical_markers,
  group.by = "annotations"
) +
  xlab(NULL) +
  ylab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "Spectral")))

ggsave(
  file.path(figure_dir, "DotPlot_canonical_markers.png"),
  p_dot,
  width = 150,
  height = 100,
  units = "mm",
  dpi = 600
)

# -----------------------------
# Detailed marker dot plot
# -----------------------------

detailed_markers <- c(
  "Ms4a2", "Cpa3",
  "Cxcr2", "S100a9", "Lcn2",
  "Ly6c2", "Siglech",
  "Mki67", "Top2a",
  "Ncr1", "Klra7",
  "Cldn1", "Cxadr", "Mpzl2",
  "Cd3g", "Lat",
  "Pecam1", "Vwf", "Eng",
  "Gas1", "Angptl1", "Cilp", "Pi16", "Dpt", "Pcolce2",
  "Acta2", "Des", "Notch3",
  "Cd209a", "Napsa",
  "Chil3", "Klra2", "Trem2", "Fcgr1", "C1qa",
  "Sox10", "Mpz", "Prx", "Slc35f1", "Scn7a",
  "Sox9", "Bmp7", "Abca6"
)

p_dot_detailed <- DotPlot(
  obj,
  features = detailed_markers,
  group.by = "detailed_annotations"
) +
  xlab(NULL) +
  ylab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "Spectral")))

ggsave(
  file.path(figure_dir, "DotPlot_detailed_markers.png"),
  p_dot_detailed,
  width = 350,
  height = 170,
  units = "mm",
  dpi = 600
)

# -----------------------------
# Mbp-Golli / Mbp-Classic dot plot
# -----------------------------

mbp_features <- c("Mbp-Classic", "Mbp-Golli")

available_mbp <- mbp_features[mbp_features %in% rownames(obj)]

if (length(available_mbp) > 0) {
  p_dot_mbp <- DotPlot(
    obj,
    features = available_mbp,
    group.by = "annotations"
  ) +
    xlab(NULL) +
    ylab(NULL) +
    coord_flip() +
    scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "Spectral"))) +
    theme_minimal(base_size = 7) +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.3),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 6, face = "bold"),
      legend.key.size = unit(2.5, "mm"),
      legend.text = element_text(size = 5),
      legend.title = element_text(size = 5)
    )

  ggsave(
    file.path(figure_dir, "DotPlot_Mbp_Golli_Classic.png"),
    p_dot_mbp,
    width = 100,
    height = 40,
    units = "mm",
    dpi = 600
  )
}

# -----------------------------
# Mbp-Golli violin plot
# -----------------------------

if ("Mbp-Golli" %in% rownames(obj)) {
  p_vln <- VlnPlot(
    obj,
    features = "Mbp-Golli",
    group.by = "annotations",
    raster = FALSE
  ) +
    xlab(NULL) +
    theme(legend.position = "none")

  ggsave(
    file.path(figure_dir, "VlnPlot_Mbp_Golli.png"),
    p_vln,
    width = 200,
    height = 90,
    units = "mm",
    dpi = 600
  )
}

message("Figure plotting complete. Figures saved to: ", figure_dir)

# =============================================================================
# Multi-DET cell-type enrichment dot plot
# =============================================================================

enrichment_path <- file.path(
  results_dir,
  "multi_DET_cell_type_enrichment_results.csv"
)                                                                 # <-- MODIFY HERE IF NEEDED

if (file.exists(enrichment_path)) {
  enrichment_df <- readr::read_csv(enrichment_path, show_col_types = FALSE) %>%
    mutate(
      contrast = factor(
        contrast,
        levels = c("C3_vs_C0", "C7_vs_C0"),
        labels = c("C3 vs C0", "C7 vs C0")
      ),
      cell_type = factor(
        cell_type,
        levels = c(
          "Schwann Cells",
          "Perineurial",
          "Epineurial",
          "Endoneurial",
          "Immune",
          "Pericytes VSMCs",
          "Endothelial"
        )
      )
    )

  p_enrichment <- ggplot(
    enrichment_df,
    aes(
      x = contrast,
      y = cell_type,
      size = odds_ratio,
      color = neg_log10_FDR
    )
  ) +
    geom_point(alpha = 0.9) +
    scale_size_continuous(name = "Odds ratio", range = c(1.5, 7)) +
    scale_color_gradient(name = "-log10(FDR)") +
    xlab(NULL) +
    ylab(NULL) +
    theme_minimal(base_size = 8) +
    theme(
      axis.text.x = element_text(size = 8, face = "bold"),
      axis.text.y = element_text(size = 8),
      panel.grid.major = element_line(linewidth = 0.2),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6)
    )

  ggsave(
    file.path(figure_dir, "DotPlot_multi_DET_cell_type_enrichment.png"),
    p_enrichment,
    width = 85,
    height = 65,
    units = "mm",
    dpi = 600
  )
} else {
  warning("Enrichment file not found: ", enrichment_path)
}
