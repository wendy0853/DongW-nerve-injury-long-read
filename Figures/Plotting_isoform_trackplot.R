#!/usr/bin/env Rscript

# =============================================================================
# Isoform structure track plotting
# =============================================================================
#
# Purpose:
#   Plot exon / intron structures for selected isoforms.
#
# Annotation strategy:
#   - Known annotated isoforms are drawn from Ensembl GRCm39 release 115 GTF to
#      ensure the most accurate structures.
#   - Novel isoforms are drawn from SQANTI3 corrected/renamed GTF files.
#
#
# Inputs:
#   - Ensembl reference GTF for known isoforms
#   - SQANTI3 corrected/renamed GTFs for novel isoforms
#   - SQANTI3 classification files containing ORF/CDS predictions
#
# Output:
#   - Isoform structure track plot
#
# =============================================================================

suppressPackageStartupMessages({
  library(ggtranscript)
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(readr)
  library(rtracklayer)
  library(stringr)
  library(RColorBrewer)
  library(tibble)
  library(tidyr)
  library(grid)
  library(ggtext)
})

# -----------------------------
# User-defined files/directories
# -----------------------------

figure_dir <- "/path/to/figure_output_directory"                  # <-- MODIFY HERE

ensembl_gtf_path <- "/path/to/Mus_musculus.GRCm39.115.chr.gtf"    # <-- MODIFY HERE

sqanti_gtf_dir <- "/path/to/SQANTI3_corrected_renamed_gtfs"       # <-- MODIFY HERE
sqanti_class_dir <- "/path/to/SQANTI3_classification_files"       # <-- MODIFY HERE

# -----------------------------
# User-defined isoforms to plot
# -----------------------------

gene_of_interest <- "Lama4"                                      # <-- MODIFY HERE

isoforms_to_plot <- c(                                           # <-- MODIFY HERE
  "Lama4-201",
  "Lama4-novel-6"
)

# Optional custom y-axis labels.
# If left NULL, isoform IDs will be used.
y_label_vec <- NULL                                               # <-- MODIFY HERE IF NEEDED

# Example:
# y_label_vec <- c(
#   "Lama4-201" = paste0(
#     "<b>Lama4-201</b><br>",
#     "<span style='font-size:6pt;'>(ENSMUST00000019992.6)</span>"
#   ),
#   "Lama4-novel-6" = "<b>Lama4-novel-6</b>"
# )

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Input checks
# -----------------------------

if (!file.exists(ensembl_gtf_path)) {
  stop("Ensembl GTF not found: ", ensembl_gtf_path)
}

if (!dir.exists(sqanti_gtf_dir)) {
  stop("SQANTI3 GTF directory not found: ", sqanti_gtf_dir)
}

if (!dir.exists(sqanti_class_dir)) {
  stop("SQANTI3 classification directory not found: ", sqanti_class_dir)
}

# =============================================================================
# Helper functions
# =============================================================================

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

split_exon_by_orf <- function(exon_row) {
  exon <- as_tibble(exon_row)

  cds_start <- exon$CDS_genomic_start
  cds_end <- exon$CDS_genomic_end
  s <- exon$start
  e <- exon$end
  strand <- exon$strand

  parts <- list()

  if (is.na(cds_start) || is.na(cds_end)) {
    exon$region_type <- "UTR"
    parts[[1]] <- exon
  } else {
    cds_min <- min(cds_start, cds_end)
    cds_max <- max(cds_start, cds_end)

    if (strand == "+") {
      if (s < cds_min) {
        parts[[length(parts) + 1]] <- exon %>%
          mutate(
            end = min(e, cds_min - 1),
            region_type = "5UTR"
          )
      }

      if (e > cds_min && s < cds_max) {
        parts[[length(parts) + 1]] <- exon %>%
          mutate(
            start = max(s, cds_min),
            end = min(e, cds_max),
            region_type = "CDS"
          )
      }

      if (e > cds_max) {
        parts[[length(parts) + 1]] <- exon %>%
          mutate(
            start = max(s, cds_max + 1),
            region_type = "3UTR"
          )
      }
    }

    if (strand == "-") {
      if (e > cds_max) {
        parts[[length(parts) + 1]] <- exon %>%
          mutate(
            start = max(s, cds_max + 1),
            region_type = "5UTR"
          )
      }

      if (s < cds_max && e > cds_min) {
        parts[[length(parts) + 1]] <- exon %>%
          mutate(
            start = max(s, cds_min),
            end = min(e, cds_max),
            region_type = "CDS"
          )
      }

      if (s < cds_min) {
        parts[[length(parts) + 1]] <- exon %>%
          mutate(
            end = min(e, cds_min - 1),
            region_type = "3UTR"
          )
      }
    }
  }

  bind_rows(parts) %>%
    filter(end >= start)
}

# =============================================================================
# Step 1: Separate known and novel isoforms
# =============================================================================

known_isoforms <- isoforms_to_plot[!str_detect(isoforms_to_plot, "novel")]
novel_isoforms <- isoforms_to_plot[str_detect(isoforms_to_plot, "novel")]

message("Known isoforms: ", paste(known_isoforms, collapse = ", "))
message("Novel isoforms: ", paste(novel_isoforms, collapse = ", "))

# =============================================================================
# Step 2: Load known isoform structures from Ensembl reference GTF
# =============================================================================

known_structure <- tibble()

if (length(known_isoforms) > 0) {
  message("Importing Ensembl reference GTF...")

  ensembl_gtf <- rtracklayer::import(ensembl_gtf_path) %>%
    as.data.frame() %>%
    mutate(
      transcript_name = if_else(
        is.na(transcript_name),
        transcript_id,
        transcript_name
      )
    )

  known_structure <- ensembl_gtf %>%
    filter(
      transcript_name %in% known_isoforms,
      type %in% c("CDS", "five_prime_utr", "three_prime_utr", "exon")
    ) %>%
    mutate(
      tx_id = transcript_name,
      region_type = case_when(
        type == "CDS" ~ "CDS",
        type %in% c("five_prime_utr", "three_prime_utr") ~ "UTR",
        type == "exon" ~ "UTR",
        TRUE ~ NA_character_
      ),
      source = "Ensembl_reference"
    ) %>%
    select(
      tx_id,
      seqnames,
      start,
      end,
      strand,
      region_type,
      source
    ) %>%
    distinct()
}

# =============================================================================
# Step 3: Load novel isoform ORF predictions from SQANTI3 classifications
# =============================================================================

novel_structure <- tibble()

if (length(novel_isoforms) > 0) {
  class_files <- list.files(
    sqanti_class_dir,
    pattern = "classification.*\\.txt$|flag_with_transcriptID\\.txt$",
    full.names = TRUE
  )

  if (length(class_files) == 0) {
    stop("No SQANTI3 classification files found in: ", sqanti_class_dir)
  }

  orf_annotation <- map_dfr(
    class_files,
    ~ readr::read_tsv(.x, show_col_types = FALSE)
  ) %>%
    select(
      transcript_id,
      strand,
      CDS_genomic_start,
      CDS_genomic_end
    ) %>%
    mutate(
      CDS_genomic_start = suppressWarnings(as.numeric(CDS_genomic_start)),
      CDS_genomic_end = suppressWarnings(as.numeric(CDS_genomic_end)),
      CDS_length = abs(CDS_genomic_end - CDS_genomic_start)
    ) %>%
    group_by(transcript_id) %>%
    arrange(desc(replace_na(CDS_length, -Inf)), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      start0 = CDS_genomic_start,
      end0 = CDS_genomic_end,
      flip = !is.na(start0) & !is.na(end0) & start0 > end0 & strand == "-"
    ) %>%
    transmute(
      transcript_id,
      CDS_genomic_start = if_else(flip, end0, start0),
      CDS_genomic_end = if_else(flip, start0, end0)
    ) %>%
    distinct(transcript_id, .keep_all = TRUE)

  # ---------------------------------------------------------------------------
  # Load SQANTI3 corrected/renamed GTF files for novel isoform exon structures.
  # ---------------------------------------------------------------------------

  sqanti_gtf_files <- list.files(
    sqanti_gtf_dir,
    pattern = "_corrected_renamed\\.gtf$|_corrected\\.gtf$|\\.gtf$",
    full.names = TRUE
  )

  if (length(sqanti_gtf_files) == 0) {
    stop("No SQANTI3 GTF files found in: ", sqanti_gtf_dir)
  }

  sqanti_gtf_all <- map_dfr(sqanti_gtf_files, function(f) {
    rtracklayer::import(f) %>%
      as.data.frame() %>%
      mutate(source_file = basename(f))
  })

  cds_from_gtf <- sqanti_gtf_all %>%
    filter(type == "CDS") %>%
    group_by(transcript_id) %>%
    summarise(
      CDS_genomic_start = min(start),
      CDS_genomic_end = max(end),
      .groups = "drop"
    )

  orf_annotation_combined <- orf_annotation %>%
    full_join(cds_from_gtf, by = "transcript_id") %>%
    mutate(
      CDS_genomic_start = coalesce(CDS_genomic_start.x, CDS_genomic_start.y),
      CDS_genomic_end = coalesce(CDS_genomic_end.x, CDS_genomic_end.y)
    ) %>%
    select(transcript_id, CDS_genomic_start, CDS_genomic_end)

  gtf_exons <- sqanti_gtf_all %>%
    filter(type == "exon") %>%
    left_join(orf_annotation_combined, by = "transcript_id")

  novel_selected <- gtf_exons %>%
    filter(transcript_id %in% novel_isoforms) %>%
    arrange(transcript_id, start) %>%
    distinct(
      seqnames,
      start,
      end,
      strand,
      transcript_id,
      gene_id,
      CDS_genomic_start,
      CDS_genomic_end,
      .keep_all = TRUE
    )

  novel_structure <- novel_selected %>%
    rowwise() %>%
    do(split_exon_by_orf(.)) %>%
    ungroup() %>%
    mutate(
      tx_id = transcript_id,
      region_type = case_when(
        region_type == "CDS" ~ "CDS",
        region_type %in% c("5UTR", "3UTR", "UTR") ~ "UTR",
        TRUE ~ NA_character_
      ),
      source = "SQANTI3_ORF_prediction"
    ) %>%
    select(
      tx_id,
      seqnames,
      start,
      end,
      strand,
      region_type,
      source
    ) %>%
    filter(!is.na(region_type)) %>%
    distinct()
}

# =============================================================================
# Step 4: Combine known and novel structures
# =============================================================================

df_split_annot <- bind_rows(
  known_structure,
  novel_structure
) %>%
  filter(tx_id %in% isoforms_to_plot) %>%
  mutate(
    tx_id = factor(tx_id, levels = isoforms_to_plot)
  )

if (nrow(df_split_annot) == 0) {
  stop("No transcript structures found for requested isoforms.")
}

# -----------------------------
# Compute introns
# -----------------------------

introns_selected <- ggtranscript::to_intron(
  df_split_annot,
  "tx_id"
) %>%
  mutate(
    tx_id = factor(tx_id, levels = isoforms_to_plot),
    arrow_start = if_else(strand == "-", end, start),
    arrow_end = if_else(strand == "-", start, end)
  )

# -----------------------------
# Colors
# -----------------------------

n_isoforms <- length(unique(df_split_annot$tx_id))
color_palette <- make_extended_palette(n_isoforms)
names(color_palette) <- isoforms_to_plot

# -----------------------------
# Labels
# -----------------------------

if (is.null(y_label_vec)) {
  y_label_vec <- setNames(
    paste0("<b>", isoforms_to_plot, "</b>"),
    isoforms_to_plot
  )
}

# =============================================================================
# Step 5: Plot isoform structures
# =============================================================================

p_isoforms <- ggplot() +
  geom_range(
    data = df_split_annot %>%
      filter(region_type %in% c("5UTR", "3UTR", "UTR")),
    aes(
      xstart = start,
      xend = end,
      y = tx_id
    ),
    fill = "gray60",
    height = 0.2,
    linewidth = 0.1
  ) +
  geom_range(
    data = df_split_annot %>%
      filter(region_type == "CDS"),
    aes(
      xstart = start,
      xend = end,
      y = tx_id,
      fill = tx_id
    ),
    height = 0.4,
    linewidth = 0.1,
    color = "black"
  ) +
  geom_intron(
    data = introns_selected,
    aes(
      xstart = arrow_start,
      xend = arrow_end,
      y = tx_id
    ),
    color = "gray60",
    linewidth = 0.25,
    arrow.min.intron.length = 100,
    arrow = grid::arrow(
      length = unit(0.8, "mm"),
      type = "open"
    )
  ) +
  scale_fill_manual(values = color_palette) +
  scale_y_discrete(
    limits = rev(levels(df_split_annot$tx_id)),
    labels = y_label_vec
  ) +
  labs(
    x = "Genomic position (bp)",
    y = NULL,
    title = gene_of_interest
  ) +
  theme_bw(base_size = 7) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 6),
    axis.text.y = ggtext::element_markdown(
      size = 7,
      color = "black",
      lineheight = 1.2
    ),
    plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

# -----------------------------
# Save figure
# -----------------------------

out_file <- file.path(
  figure_dir,
  paste0(gene_of_interest, "_IsoformTrackplot.png")
)

ggsave(
  out_file,
  p_isoforms,
  width = 90,
  height = 40,
  units = "mm",
  dpi = 600
)

message("Saved isoform track plot: ", out_file)
