#!/usr/bin/env Rscript

# =============================================================================
# Single-nucleus RNA-seq integration and pseudobulk DESeq2 analysis
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  options(Seurat.object.assay.version = "v5")
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(DESeq2)
  library(harmony)
  library(SeuratWrappers)
  library(speckle)
})

# -----------------------------
# User-defined directories
# -----------------------------

data_base <- "/path/to/cellbender_filtered_h5_files"          # <-- MODIFY HERE
scrub_base <- "/path/to/scrublet_doublet_score_files"         # <-- MODIFY HERE
results_dir <- "/path/to/single_cell_results"                 # <-- MODIFY HERE

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# User-defined sample metadata
# -----------------------------

samples <- c(                                                  # <-- MODIFY HERE IF NEEDED
  "C0_1", "C0_2", "C0_3",
  "C3_1", "C3_2", "C3_3", "C3_4", "C3_5",
  "C7_1", "C7_2", "C7_3", "C7_4", "C7_5", "C7_6", "C7_7", "C7_8"
)

# -----------------------------
# QC parameters
# -----------------------------

min_features <- 500                                            # <-- MODIFY HERE IF NEEDED
max_features <- 8000                                           # <-- MODIFY HERE IF NEEDED
min_counts <- 10                                               # <-- MODIFY HERE IF NEEDED
max_percent_mt <- 15                                           # <-- MODIFY HERE IF NEEDED
mt_pattern <- "^MT-|^mt-"                                      # <-- MODIFY HERE IF NEEDED

# -----------------------------
# Helper functions
# -----------------------------

get_timepoint <- function(sample_id) {
  sub("_.*$", "", sample_id)
}

find_elbow <- function(seurat_obj, reduction = "pca", max_pcs = 50) {
  stdev <- seurat_obj[[reduction]]@stdev[1:max_pcs]
  pct_change <- diff(stdev) / stdev[-length(stdev)]
  elbow <- which(abs(pct_change) < 0.05)[1]
  ifelse(is.na(elbow), 30, elbow)
}

process_one_sample <- function(sample_id) {
  message("Processing sample: ", sample_id)

  h5_file <- file.path(data_base, sample_id, "cellbender_output_file_filtered_seurat.h5") # <-- MODIFY HERE IF FILE NAME DIFFERS
  scrub_file <- file.path(scrub_base, paste0(sample_id, "_doublet_scores.csv"))           # <-- MODIFY HERE IF FILE NAME DIFFERS

  if (!file.exists(h5_file)) {
    stop("Missing CellBender filtered h5 file: ", h5_file)
  }

  if (!file.exists(scrub_file)) {
    stop("Missing Scrublet score file: ", scrub_file)
  }

  counts <- Read10X_h5(filename = h5_file, use.names = TRUE)

  seu <- CreateSeuratObject(
    counts = counts,
    project = sample_id
  )

  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = mt_pattern)
  seu$timepoint <- get_timepoint(sample_id)
  seu$sample_id <- sample_id

  seu <- subset(
    seu,
    subset =
      nFeature_RNA > min_features &
      nFeature_RNA < max_features &
      nCount_RNA > min_counts &
      percent.mt < max_percent_mt
  )

  scrub <- read.csv(scrub_file, stringsAsFactors = FALSE)
  rownames(scrub) <- scrub$barcode

  scrub <- scrub[colnames(seu), , drop = FALSE]
  seu <- AddMetaData(seu, metadata = scrub)

  seu <- subset(seu, subset = Predicted_Doublets == "False")

  return(seu)
}

# =============================================================================
# Step 1: Load and QC individual samples
# =============================================================================

seu_list <- setNames(
  lapply(samples, process_one_sample),
  samples
)

cell_counts <- tibble(
  sample = names(seu_list),
  n_cells_after_QC = sapply(seu_list, ncol),
  n_genes = sapply(seu_list, nrow)
)

write_csv(cell_counts, file.path(results_dir, "cell_counts_after_QC.csv"))

# =============================================================================
# Step 2: Merge samples
# =============================================================================

obj <- merge(
  x = seu_list[[1]],
  y = seu_list[-1],
  add.cell.ids = names(seu_list),
  project = "single_cell_merged"
)

# =============================================================================
# Step 3: Initial normalization, PCA, Harmony integration, and clustering
# =============================================================================

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 50)

obj <- IntegrateLayers(
  object = obj,
  orig.reduction = "pca",
  method = HarmonyIntegration,
  new.reduction = "harmony",
  verbose = FALSE
)

optimal_pcs <- find_elbow(obj, reduction = "pca", max_pcs = 50)

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:optimal_pcs)
obj <- FindClusters(obj, resolution = 0.5, cluster.name = "harmony_clusters")
obj <- RunUMAP(
  obj,
  reduction = "harmony",
  dims = 1:optimal_pcs,
  reduction.name = "umap.harmony"
)

# =============================================================================
# Step 4 (Optional): Remove and re-integrate afterr noise-cluster removal if applicable
# =============================================================================

noise_clusters <- c("17")                                      # <-- MODIFY HERE IF NEEDED
# Noise cluster confirmed based on high doblet score, ambient RNA as top differentially
# expressed genes, and disperse clustering

if (any(obj$seurat_clusters %in% noise_clusters)) {
  obj <- subset(obj, seurat_clusters %in% setdiff(unique(obj$seurat_clusters), noise_clusters))
}

obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 50)

obj <- IntegrateLayers(
  object = obj,
  orig.reduction = "pca",
  method = HarmonyIntegration,
  new.reduction = "harmony",
  verbose = FALSE
)

optimal_pcs <- find_elbow(obj, reduction = "pca", max_pcs = 50)

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:optimal_pcs)
obj <- FindClusters(obj, resolution = 0.3, cluster.name = "harmony_clusters")
obj <- RunUMAP(
  obj,
  reduction = "harmony",
  dims = 1:optimal_pcs,
  reduction.name = "umap.harmony"
)

# =============================================================================
# Step 5: Marker genes and cluster annotation
# =============================================================================

obj <- JoinLayers(obj)

allClusterMarkers <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

allClusterMarkers$pct.diff <- allClusterMarkers$pct.1 - allClusterMarkers$pct.2

write_csv(
  allClusterMarkers,
  file.path(results_dir, "Cluster_Markers.csv")
)

# Broad annotations
Idents(obj) <- obj$harmony_clusters

obj <- RenameIdents(
  obj,
  "0" = "Endoneurial",
  "1" = "Immune",
  "2" = "Epineurial",
  "3" = "Schwann Cells",
  "4" = "Immune",
  "5" = "Immune",
  "6" = "Pericytes/VSMCs",
  "7" = "Epineurial",
  "8" = "Endothelial",
  "9" = "Immune",
  "10" = "Schwann Cells",
  "11" = "Perineurial",
  "12" = "Immune",
  "13" = "Immune",
  "14" = "Immune",
  "15" = "Immune",
  "16" = "Immune"
)

obj$annotations <- Idents(obj)

# Detailed annotations
Idents(obj) <- obj$harmony_clusters

obj <- RenameIdents(
  obj,
  "0" = "Endoneurial",
  "1" = "Macrophages-c1",
  "2" = "Epineurial-c1",
  "3" = "Schwann Cells",
  "4" = "Macrophages-c2",
  "5" = "MoDC",
  "6" = "Pericytes/VSMCs",
  "7" = "Epineurial-c2",
  "8" = "Endothelial",
  "9" = "T Cells",
  "10" = "Prolif. SC",
  "11" = "Perineurial",
  "12" = "NK Cells",
  "13" = "Prolif. Immune",
  "14" = "pDC",
  "15" = "Granulocytes",
  "16" = "Mast Cells"
)

obj$detailed_annotations <- Idents(obj)

saveRDS(
  obj,
  file.path(results_dir, "integrated_annotated_seurat_object.rds")
)

# =============================================================================
# Step 6: Cell proportion analysis using Propeller
# =============================================================================

propeller_results <- propeller(
  clusters = obj$annotations,
  sample = obj$orig.ident,
  group = obj$timepoint,
  transform = "asin"
)

write_csv(
  as.data.frame(propeller_results),
  file.path(results_dir, "propeller_cell_proportion_results.csv")
)

# =============================================================================
# Step 7: Pseudobulk aggregation
# =============================================================================

pb <- AggregateExpression(
  obj,
  group.by = c("annotations", "timepoint", "orig.ident"),
  assays = "RNA",
  slot = "counts",
  return.seurat = FALSE
)

pb_mat <- pb$RNA

parts <- strsplit(colnames(pb_mat), "_")

meta <- data.frame(
  sample = colnames(pb_mat),
  annotations = vapply(parts, `[`, character(1), 1),
  timepoint = vapply(parts, `[`, character(1), 2),
  orig.ident = vapply(parts, function(x) paste(x[3:length(x)], collapse = "_"), character(1)),
  row.names = colnames(pb_mat),
  stringsAsFactors = FALSE
)

meta$timepoint <- factor(meta$timepoint, levels = c("C0", "C3", "C7"))

stopifnot(identical(colnames(pb_mat), rownames(meta)))

write.csv(
  meta,
  file.path(results_dir, "pseudobulk_sample_metadata.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(pb_mat),
  file.path(results_dir, "pseudobulk_count_matrix.csv")
)

# =============================================================================
# Step 8: Pseudobulk DESeq2 by cell type
# =============================================================================

run_within_compartment_DE <- function(annotation_name,
                                      min_count = 10,
                                      min_n_samples = 2) {
  keep <- meta$annotations == annotation_name

  if (sum(keep) < 3) {
    warning("Not enough pseudobulk samples for: ", annotation_name)
    return(NULL)
  }

  meta_sub <- droplevels(meta[keep, , drop = FALSE])
  mat_sub <- pb_mat[, keep, drop = FALSE]

  stopifnot(identical(colnames(mat_sub), rownames(meta_sub)))

  keep_genes <- rowSums(mat_sub >= min_count) >= min_n_samples
  mat_sub <- mat_sub[keep_genes, , drop = FALSE]

  dds <- DESeqDataSetFromMatrix(
    countData = round(mat_sub),
    colData = meta_sub,
    design = ~ timepoint
  )

  dds <- DESeq(dds)

  list(
    C3_vs_C0 = results(dds, contrast = c("timepoint", "C3", "C0")),
    C7_vs_C0 = results(dds, contrast = c("timepoint", "C7", "C0"))
  )
}

to_df <- function(res) {
  df <- as.data.frame(res)
  df$gene <- rownames(df)
  df
}

cell_types <- c(
  "Schwann Cells",
  "Perineurial",
  "Epineurial",
  "Endoneurial",
  "Immune",
  "Pericytes/VSMCs",
  "Endothelial"
)

contrasts <- c("C3_vs_C0", "C7_vs_C0")

for (ct in cell_types) {
  message("Running pseudobulk DESeq2 for: ", ct)

  de <- run_within_compartment_DE(
    annotation_name = ct,
    min_count = 10,
    min_n_samples = 2
  )

  if (is.null(de)) next

  ct_safe <- str_replace_all(ct, "[^A-Za-z0-9]+", "_")

  for (con in contrasts) {
    res <- to_df(de[[con]])

    out_file <- file.path(
      results_dir,
      paste0("DE_", ct_safe, "_", con, ".csv")
    )

    write_csv(res, out_file)
  }
}


  
