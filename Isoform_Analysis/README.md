# Isoform-Level Transcriptomic Analysis

This directory contains scripts used for isoform-level transcriptomic analyses of long-read RNA-seq data, including:

- differential transcript expression (DTE)
- differential transcript usage (DTU)
- isoform switching
- transcript emergence classification
- structural remodeling analysis

These analyses were performed using transcript models generated from IsoQuant and refined through SQANTI3 filtering and rescue workflows.

---

# Workflow Overview

The isoform analysis workflow was performed in the following order:

| Step | Script | Purpose |
|---|---|---|
| 1 | `Unifying_Novel_Transcript_Names.R` | Harmonize novel transcript identifiers across samples |
| 2 | `Generating_Count_Matrix.R` | Generate transcript-level count matrices |
| 3 | `DTE_Analysis_DESeq2.R` | Differential transcript expression analysis |
| 4 | `DET_Remodeling_Classification.R` | Transcript emergence and remodeling classification |
| 5 | `DTU_Analysis_IsoformSwitchAnalyzeR_DEXSeq.R` | Differential transcript usage and isoform switching analysis |

---

# Differential Transcript Expression (DTE)

Differential transcript expression analysis was performed using DESeq2 v1.44.0 in R.

Input:
- isoform-level count matrices derived from long-read RNA-seq

Isoforms with low expression were filtered prior to analysis:
- retained isoforms required ≥10 counts in at least two biological replicates of any condition

DESeq2 analysis included:
- median-of-ratios normalization
- negative binomial generalized linear modeling
- Wald tests
- Benjamini-Hochberg multiple testing correction

Significance thresholds:
- adjusted p-value ≤ 0.05
- |log2FoldChange| ≥ 1

Because transcript-level analyses exhibit increased sparsity and variance relative to gene-level analyses, Cook’s distance filtering was disabled (cooksCutoff = FALSE) to avoid exclusion of biologically meaningful transcript-specific changes.

---

# Differential Transcript Usage (DTU)

DTU analysis was performed using:

- IsoformSwitchAnalyzeR v2.4.0
- DEXSeq v1.50.0

Input:
- filtered transcript annotations generated from IsoQuant and SQANTI3

Analyses included:
- isoform switching
- differential isoform fraction usage
- coding sequence alterations
- domain gain/loss
- transcript structure remodeling

Significance thresholds:
- adjusted p-value ≤ 0.05
- differential isoform fraction (dIF) ≥ 0.1

---

# DET Emergence and Remodeling Classification

Differentially expressed transcripts (DETs) were classified to distinguish:

- transcript-level emergence
- gene-level activation
- isoform remodeling within actively transcribed loci

An isoform was considered expressed at baseline if:
- ≥10 counts were detected in at least two C0 biological replicates

DETs were assigned into three categories:

| Category | Description |
|---|---|
| Expressed Baseline | DET expressed at baseline |
| Emerged: Activation | DET absent at baseline and no isoforms from same gene expressed |
| Emerged: Remodeling | DET absent at baseline but alternative isoforms from same gene expressed |

---

# Structural Remodeling Classification

For DETs classified as:

```text
Emerged: Remodeling
```

each transcript was compared against the dominant baseline isoform of the same gene.

The dominant baseline isoform was defined as:
- the highest mean expressed isoform across C0 samples

Structural consequences were classified as:

| Category | Description |
|---|---|
| UTR Only | Changes restricted to untranslated regions |
| ORF Modification | Coding sequence structural changes |
| ORF Gain | Baseline isoform noncoding → DET coding |
| ORF Loss | Baseline isoform coding → DET noncoding |
| Noncoding | Both transcripts noncoding |

Known transcript CDS annotations were obtained from:
- Ensembl GRCm39 release 115

Novel transcript ORFs were inferred using:
- SQANTI3 ORF predictions

---

# Functional Enrichment Analysis

Gene Ontology enrichment analyses were performed using:

- clusterProfiler v4.12.6

GO analyses were performed for:
- multi-isoform regulated genes
- DET-only genes

---

# Main Outputs

Typical outputs include:

| Output | Description |
|---|---|
| transcript count matrices | isoform-level quantification |
| DTE result tables | transcript-level differential expression |
| DTU result tables | differential isoform usage |
| transcript emergence annotations | DET category assignments |
| structural remodeling annotations | ORF/UTR remodeling classifications |
| GO enrichment tables | functional enrichment analyses |

---

# Software

Primary R packages used:

- DESeq2
- IsoformSwitchAnalyzeR
- DEXSeq
- clusterProfiler
- tidyverse
- GenomicRanges
- biomaRt
- AnnotationDbi
