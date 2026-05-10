# Short-Read RNA-seq Differential Gene Expression Analysis

This directory contains code used for differential gene expression analysis of Illumina short-read RNA-seq data.

---

## Overview

Short-read RNA-seq libraries were prepared from the same bulk RNA samples used for long-read sequencing.

Library preparation and sequencing were performed by GTAC@MGI using:

- SMARTer Ultra Low RNA Kit for Illumina Sequencing
- Illumina NovaSeq X Plus
- at least 50 million paired-end reads per sample

Base calling and demultiplexing were performed using Illumina bcl2fastq.

Read alignment and gene-level quantification were performed by the sequencing core using:

- STAR v2.7.11b
- featureCounts v2.0.8

Reference files:

- Genome: GRCm39
- Annotation: GENCODE vM38

---

## Sample Design

The short-read dataset contained 15 samples across five groups:

| Group | Description |
|---|---|
| C0_Baseline | Uninjured baseline sciatic nerve |
| C3_Control | Uninjured contralateral side from C3 injury mice |
| C3_Injured | Injured sciatic nerve at 3 days post crush |
| C7_Control | Uninjured contralateral side from C7 injury mice |
| C7_Injured | Injured sciatic nerve at 7 days post crush |

---

## Differential Gene Expression Analysis

Differential gene expression analysis was performed in R using DESeq2 v1.44.0.

Genes with low counts were filtered before analysis, retaining genes with:

- ≥10 counts in at least two biological replicates of any condition

DESeq2 analysis included:

- median-of-ratios normalization
- negative binomial generalized linear modeling
- Wald tests
- Benjamini-Hochberg multiple testing correction

Significant DEGs were defined as:

- adjusted p-value ≤ 0.05
- |log2FoldChange| ≥ 1

---

## Main Comparisons

The following contrasts were tested:

| Comparison | Description |
|---|---|
| C3_Injured vs C3_Control | Injury effect at C3 relative to matched uninjured side |
| C7_Injured vs C7_Control | Injury effect at C7 relative to matched uninjured side |
| C3_Injured vs C0_Baseline | Injury effect at C3 relative to baseline |
| C7_Injured vs C0_Baseline | Injury effect at C7 relative to baseline |
| C7_Injured vs C3_Injured | Temporal injury progression |

Injury-specific DEGs were additionally identified as genes significant in both:

- injured vs matched contralateral control
- injured vs C0 baseline

---

## Input Files

Required inputs:

| File | Description |
|---|---|
| `all.gene_counts.xlsx` or `.csv` | featureCounts gene-level count matrix |
| sample metadata | sample group, time point, and treatment labels |

The script assumes the count matrix contains gene annotation columns followed by sample count columns.

---

## Running the Analysis

Example:

```bash
Rscript run_deseq2_short_read.R \
  --counts all.gene_counts.xlsx \
  --outdir results
```

---

## Outputs

The script generates:

| Output | Description |
|---|---|
| `*_results.csv` | DESeq2 results for each comparison |
| `*_significant_DEGs.csv` | Significant DEGs passing padj and log2FC thresholds |
| `normalized_counts.csv` | DESeq2 normalized counts |
| `PCA_plot_nerve_injury.png` | PCA plot of short-read samples |
| `MA_plot_*.png` | MA plots for each comparison |
| `Volcano_plot_*.png` | Volcano plots for each comparison |

---

## Software

Primary R packages:

- DESeq2
- readxl
- tidyverse
- ggplot2
- EnhancedVolcano
- pheatmap


---

## Reference

DESeq2:
https://bioconductor.org/packages/release/bioc/html/DESeq2.html
