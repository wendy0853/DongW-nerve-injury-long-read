#!/usr/bin/env Rscript

################################################################################
# Figure 1: Long-read RNA-seq reveals high-quality transcriptome profiling and 
#  isoform remodeling after peripheral nerve injury. 
#
# Panels:
#   A: Project schematic
#   B: Isoform length distribution
#   C: Exon count distribution
#   D: Isoforms per gene distribution
#   E: SQANTI3 structural category distribution
#   F: Total novel isoforms
#   G: Novel coding isoforms
#
# Inputs:
#   - SQANTI3 classification files
#   - SQANTI3 structural category summary table
#
# Outputs:
#   - Publication-ready Figure 1 panels
#
################################################################################

# ==============================================================================
# Setup
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(grid)
})

# ==============================================================================
# User-defined directories
# ==============================================================================

# Directory containing SQANTI3 classification outputs
sqanti_dir <- "/path/to/SQANTI3_results"                         # <-- MODIFY HERE

# Directory containing processed sample metric summaries
metrics_dir <- "/path/to/sample_metrics"                         # <-- MODIFY HERE

# Output directory for Figure 1 panels
figure_dir <- "/path/to/Figure_1_output"                         # <-- MODIFY HERE

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Input files
# ==============================================================================

# SQANTI3 classification files generated after filtering/rescue
classification_files <- c(                                      # <-- MODIFY HERE IF NEEDED
  "Filtered_SQANTI3_Results/C0_Sciatic_1_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt",
  "Filtered_SQANTI3_Results/C0_Sciatic_2_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt",
  "Filtered_SQANTI3_Results/C0_Sciatic_3_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt",
  "Filtered_SQANTI3_Results/C3_Injured_Sciatic_1_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt",
  "Filtered_SQANTI3_Results/C3_Injured_Sciatic_2_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt",
  "Filtered_SQANTI3_Results/C3_Injured_Sciatic_3_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt",
  "Filtered_SQANTI3_Results/C7_Injured_Sciatic_1_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt",
  "Filtered_SQANTI3_Results/C7_Injured_Sciatic_2_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt",
  "Filtered_SQANTI3_Results/C7_Injured_Sciatic_3_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt"
)

classification_paths <- file.path(
  sqanti_dir,
  classification_files
)

# SQANTI3 structural category summary table
category_summary_path <- file.path(
  metrics_dir,
  "Filtered_SQANTI3_Category_summary_FINAL.csv"
)                                                                # <-- MODIFY HERE IF NEEDED

# ==============================================================================
# Input validation
# ==============================================================================

missing_class_files <- classification_paths[!file.exists(classification_paths)]

if (length(missing_class_files) > 0) {
  stop(
    "Missing SQANTI3 classification file(s):\n",
    paste(missing_class_files, collapse = "\n"),
    call. = FALSE
  )
}

if (!file.exists(category_summary_path)) {
  stop(
    "Missing category summary file:\n",
    category_summary_path,
    call. = FALSE
  )
}

# ==============================================================================
# Colors
# ==============================================================================

condition_fill <- c(
  "C0 Control" = "#A6CEE3",
  "C3 Injured" = "#FFD580",
  "C7 Injured" = "#B2DF8A"
)

condition_line <- c(
  "C0 Control" = "#0072B2",
  "C3 Injured" = "#E69F00",
  "C7 Injured" = "#66BD63"
)

condition_fill_short <- c(
  "C0" = "#0072B2",
  "C3" = "#E69F00",
  "C7" = "#66BD63"
)

condition_levels <- c(
  "C0 Control",
  "C3 Injured",
  "C7 Injured"
)

condition_labels <- c(
  "C0 Control" = "C0",
  "C3 Injured" = "C3",
  "C7 Injured" = "C7"
)

# ==============================================================================
# Shared plotting theme
# ==============================================================================

theme_pub <- function(base_size = 7) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.3),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(0.15, "cm"),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 5.5, vjust = 0.8),
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
      plot.margin = margin(
        t = 8,
        r = 8,
        b = 4,
        l = 4,
        unit = "mm"
      )
    )
}

# ==============================================================================
# Helper functions
# ==============================================================================

save_panel <- function(
  plot,
  filename,
  width_mm = 50,
  height_mm = 50
) {

  ggsave(
    filename = file.path(figure_dir, filename),
    plot = plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 600
  )
}

get_sample_name <- function(path) {
  basename(path) %>%
    str_remove(
      "_filtered_RulesFilter_result_classification_with_rescue_flag_with_transcriptID.txt"
    )
}

get_condition_short <- function(x) {
  case_when(
    str_detect(x, "C0") ~ "C0",
    str_detect(x, "C3") ~ "C3",
    str_detect(x, "C7") ~ "C7",
    TRUE ~ NA_character_
  )
}

get_condition_long <- function(x) {
  case_when(
    str_detect(x, "C0") ~ "C0 Control",
    str_detect(x, "C3") ~ "C3 Injured",
    str_detect(x, "C7") ~ "C7 Injured",
    TRUE ~ NA_character_
  )
}

# ==============================================================================
# Read and format filtered SQANTI3 classification data
# ==============================================================================

classification_df <- map_dfr(classification_paths, function(path) {
  read_tsv(path, show_col_types = FALSE) %>%
    filter(filter_rescued_result == "Isoform") %>%
    transmute(
      sample = get_sample_name(path),
      condition_short = get_condition_short(path),
      condition = get_condition_long(path),
      gene = gene_symbol,
      transcript_id = transcript_id,
      length = length,
      exons = exons
    )
}) %>%
  mutate(
    condition = factor(condition, levels = condition_levels),
    condition_short = factor(condition_short, levels = c("C0", "C3", "C7"))
  )

# ==============================================================================
# Figure 1B: Isoform length distribution
# ==============================================================================

length_trim <- classification_df %>%
  group_by(condition_short) %>%
  mutate(q9999 = quantile(length, 0.9999, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(length <= q9999)

p_iso_len_box <- ggplot(length_trim, aes(x = condition_short, y = length / 1000, fill = condition_short)) +
  geom_boxplot(
    width = 0.6,
    linewidth = 0.2,
    outlier.size = 0.02,
    outlier.alpha = 0.1
  ) +
  scale_fill_manual(values = condition_fill_short) +
  labs(
    x = NULL,
    y = "Isoform length (kb)",
    title = "Isoform Length"
  ) +
  theme_pub() +
  theme(legend.position = "none")

p_iso_len_box
save_panel(p_iso_len_box, "Fig1B_Isoform_Length.png")

# ==============================================================================
# Figure 1C: Exon count distribution
# ==============================================================================

exon_trim <- classification_df %>%
  group_by(condition_short) %>%
  mutate(q9999 = quantile(exons, 0.9999, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(exons <= q9999)

p_iso_exon_box <- ggplot(exon_trim, aes(x = condition_short, y = exons, fill = condition_short)) +
  geom_boxplot(
    width = 0.6,
    linewidth = 0.2,
    outlier.size = 0.02,
    outlier.alpha = 0.1
  ) +
  scale_fill_manual(values = condition_fill_short) +
  labs(
    x = NULL,
    y = "Exon count",
    title = "Number of Exons"
  ) +
  theme_pub() +
  theme(legend.position = "none")

p_iso_exon_box
save_panel(p_iso_exon_box, "Fig1C_Exon_Count.png")

# ==============================================================================
# 4. Figure 1D: Distribution of isoform counts per gene
# ==============================================================================

gene_iso_summary <- classification_df %>%
  group_by(condition, sample, gene) %>%
  summarise(n_isoforms = n_distinct(transcript_id), .groups = "drop") %>%
  mutate(isoform_bin = if_else(n_isoforms > 5, "6+", as.character(n_isoforms))) %>%
  group_by(condition, sample, isoform_bin) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(condition, sample) %>%
  mutate(percent = 100 * count / sum(count)) %>%
  ungroup() %>%
  group_by(condition, isoform_bin) %>%
  summarise(
    mean_percent = mean(percent, na.rm = TRUE),
    sd_percent = sd(percent, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(isoform_bin = factor(isoform_bin, levels = c(as.character(1:5), "6+")))

p_isoform_count <- ggplot(gene_iso_summary, aes(x = isoform_bin, y = mean_percent, fill = condition)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black",
    linewidth = 0.2
  ) +
  geom_errorbar(
    aes(ymin = mean_percent - sd_percent, ymax = mean_percent + sd_percent),
    position = position_dodge(width = 0.8),
    width = 0.3,
    linewidth = 0.2
  ) +
  scale_fill_manual(values = condition_line) +
  scale_y_continuous(limits = c(0, 70), breaks = seq(0, 70, 10), expand = c(0, 0)) +
  labs(
    x = "Number of isoforms per gene",
    y = "Percentage of genes (%)",
    title = "Distribution of Isoform Counts per Gene"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(size = 5, vjust = 0.8),
    legend.position = c(0.9, 0.75),
    legend.background = element_rect(color = "black", fill = "white", linewidth = 0.3)
  )

p_isoform_count
save_panel(p_isoform_count, "Fig1D_Isoforms_Per_Gene.png", width_mm = 95, height_mm = 55)

# ==============================================================================
# Figure 1E: SQANTI3 structural category distribution
# ==============================================================================

# Read SQANTI3 structural category summary
category_df <- read_csv(category_summary_path, show_col_types = FALSE) %>%
  mutate(
    condition = case_when(
      str_detect(Sample_ID, "^C0_") ~ "C0 Control",
      str_detect(Sample_ID, "^C3_") ~ "C3 Injured",
      str_detect(Sample_ID, "^C7_") ~ "C7 Injured",
      TRUE ~ NA_character_
    ),
    condition = factor(condition, levels = condition_levels)
  )

structural_summary <- category_df %>%
  group_by(condition) %>%
  summarise(
    FSM_coding = mean(`full-splice_match_coding_pct`, na.rm = TRUE),
    FSM_non_coding = mean(`full-splice_match_non_coding_pct`, na.rm = TRUE),
    ISM_coding = mean(`incomplete-splice_match_coding_pct`, na.rm = TRUE),
    ISM_non_coding = mean(`incomplete-splice_match_non_coding_pct`, na.rm = TRUE),
    NIC_coding = mean(`novel_in_catalog_coding_pct`, na.rm = TRUE),
    NIC_non_coding = mean(`novel_in_catalog_non_coding_pct`, na.rm = TRUE),
    NNC_coding = mean(`novel_not_in_catalog_coding_pct`, na.rm = TRUE),
    NNC_non_coding = mean(`novel_not_in_catalog_non_coding_pct`, na.rm = TRUE),
    .groups = "drop"
  )

structural_long <- structural_summary %>%
  pivot_longer(
    cols = -condition,
    names_to = c("category", "biotype"),
    names_pattern = "^(FSM|ISM|NIC|NNC)_(coding|non_coding)$",
    values_to = "percentage"
  ) %>%
  mutate(
    category = factor(category, levels = c("FSM", "ISM", "NIC", "NNC")),
    biotype = recode(biotype, coding = "Coding", non_coding = "Non-coding")
  ) %>%
  pivot_wider(names_from = biotype, values_from = percentage) %>%
  arrange(category, condition) %>%
  group_by(category) %>%
  mutate(
    bar_id = row_number(),
    total = Coding + `Non-coding`
  ) %>%
  ungroup() %>%
  pivot_longer(cols = c(Coding, `Non-coding`), names_to = "biotype", values_to = "percentage") %>%
  mutate(
    biotype = factor(biotype, levels = c("Non-coding", "Coding")),
    group_id = as.numeric(category),
    x_position = (group_id - 1) * 4 + bar_id,
    fill_group = interaction(bar_id, biotype, sep = ".")
  )

structural_labels <- structural_long %>%
  filter(biotype == "Non-coding") %>%
  distinct(x_position, total)

category_labels <- structural_long %>%
  filter(bar_id == 2, biotype == "Coding") %>%
  distinct(x_position, category)

structural_fill <- c(
  "1.Coding" = "#0072B2", "1.Non-coding" = "#A6CEE3",
  "2.Coding" = "#E69F00", "2.Non-coding" = "#FFD580",
  "3.Coding" = "#66BD63", "3.Non-coding" = "#B2DF8A"
)

structural_fill_labels <- c(
  "1.Coding" = "C0 Coding", "1.Non-coding" = "C0 Non-coding",
  "2.Coding" = "C3 Coding", "2.Non-coding" = "C3 Non-coding",
  "3.Coding" = "C7 Coding", "3.Non-coding" = "C7 Non-coding"
)

p_structural <- ggplot(structural_long, aes(x = x_position, y = percentage, fill = fill_group)) +
  geom_col(width = 1, color = "black", linewidth = 0.2) +
  geom_text(
    data = structural_labels,
    aes(x = x_position, y = total + 1.5, label = paste0(round(total, 1), "%")),
    inherit.aes = FALSE,
    size = 4 * 0.39
  ) +
  scale_fill_manual(values = structural_fill, labels = structural_fill_labels) +
  scale_x_continuous(
    breaks = category_labels$x_position,
    labels = category_labels$category,
    expand = c(0, 0.5)
  ) +
  scale_y_continuous(breaks = seq(0, 100, 20), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(
    x = NULL,
    y = "Percentage of isoforms (%)",
    title = "Distribution of Isoforms Across Structural Categories",
    fill = NULL
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(size = 5),
    legend.position = c(1, 0.5),
    legend.justification = c("right", "center"),
    legend.background = element_rect(color = "black", fill = "white", linewidth = 0.3),
    legend.box.margin = margin(t = -20, b = -5, unit = "mm"),
    legend.spacing.y = unit(1, "mm"),
    legend.spacing.x = unit(1, "mm")
  )

p_structural
save_panel(p_structural, "Fig1E_SQANTI3_Structural_Categories.png", width_mm = 95, height_mm = 55)

# ==============================================================================
# Figure 1F/G: Total novel isoforms, Novel Coding Isoforms
# ==============================================================================

# Shared function for novel isoform panels
make_novel_df <- function(df, novel_type = c("total", "coding")) {
  novel_type <- match.arg(novel_type)
  
  if (novel_type == "total") {
    novel_cols <- names(df)[str_detect(names(df), "^novel_.*_pct$")]
  } else {
    novel_cols <- names(df)[str_detect(names(df), "^novel_.*catalog_coding_pct$")]
  }
  
  df %>%
    rowwise() %>%
    mutate(novel_percent = sum(c_across(all_of(novel_cols)), na.rm = TRUE)) %>%
    ungroup() %>%
    transmute(
      sample = Sample_ID,
      condition = condition,
      novel_percent = novel_percent
    )
}

plot_novel_panel <- function(df, title, y_breaks = waiver(), label_y = NULL) {
  comparisons <- list(
    c("C0 Control", "C3 Injured"),
    c("C0 Control", "C7 Injured")
  )
  
  p <- ggplot(df, aes(x = condition, y = novel_percent)) +
    geom_point(
      aes(fill = condition, color = condition),
      shape = 21,
      size = 1,
      stroke = 0.5,
      position = position_jitter(width = 0.1, seed = 1)
    ) +
    stat_summary(
      fun = mean,
      geom = "crossbar",
      aes(color = condition),
      width = 0.7,
      linewidth = 0.3
    ) +
    scale_fill_manual(values = condition_fill) +
    scale_color_manual(values = condition_line) +
    scale_x_discrete(labels = condition_labels) +
    scale_y_continuous(breaks = y_breaks) +
    labs(
      x = NULL,
      y = "Percentage of isoforms (%)",
      title = title
    ) +
    theme_pub() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 5.5, vjust = 0.8, hjust = 0.8)
    )
  
  if (!is.null(label_y)) {
    p <- p + stat_compare_means(
      comparisons = comparisons,
      method = "t.test",
      label = "p.format",
      label.y = label_y,
      size = 1.5
    )
  }
  
  p
}

# Figure 1F: Total novel coding isoforms

novel_total_df <- make_novel_df(category_df, novel_type = "total")

p_novel_total <- plot_novel_panel(
  novel_total_df,
  title = "Total Novel Isoforms",
  y_breaks = seq(9, 12, by = 0.5),
  label_y = c(11.44, 11.6)
)

p_novel_total
save_panel(p_novel_total, "Fig1F_Total_Novel_Isoforms.png")

# Figure 1G: Novel coding isoforms

novel_coding_df <- make_novel_df(category_df, novel_type = "coding")

p_novel_coding <- plot_novel_panel(
  novel_coding_df,
  title = "Novel Coding Isoforms",
  label_y = c(8, 8.16)
)

p_novel_coding
save_panel(p_novel_coding, "Fig1G_Novel_Coding_Isoforms.png")
