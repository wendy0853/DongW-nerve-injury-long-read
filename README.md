# Isoform Remodeling Shapes Peripheral Nerve Response to Injury

Companion repository for the submitted manuscript **Dong W et al. Peripheral Nerve Injury Engages Transcript Isoform Remodeling to Coordinate Cellular and Extracellular Response. (2026)**

This repository contains the computational workflows, analysis scripts, and supporting resources used to generate and analyze long-read and short-read transcriptomic datasets from mouse sciatic nerve crush injury.

The study integrates:

- Oxford Nanopore long-read RNA sequencing
- Illumina short-read RNA sequencing
- Single-cell and pseudobulk transcriptomic analyses

to characterize isoform-level remodeling during Wallerian degeneration.

The background and tested gene lists for every enrichment analysis in the manuscript are provided in [`Enrichment_Analysis/`](Enrichment_Analysis/).

An interactive companion browser for exploring transcript- and gene-level results is available at **IsoNerve** (https://isonerve.pages.dev)

If you use this repository or analysis framework, please cite:

> Dong W et al. *Peripheral Nerve Injury Engages Transcript Isoform Remodeling to Coordinate Cellular and Extracellular Response.* (2026)

## Study Design

Adult mouse sciatic nerves were profiled at:

| Condition | Description |
|---|---|
| C0 | Uninjured control |
| C3 | 3 days post crush injury |
| C7 | 7 days post crush injury |

Long-read and short-read sequencing were generated from the same bulk RNA samples collected from distal sciatic nerve segments following injury.

---

## Repository Structure

```text
Figures/               Manuscript figures and supplementary figures
IsoQuant/              Long-read alignment and transcript reconstruction workflows
SQANTI3/               Isoform QC, filtering, annotation, and quantification workflows
Isoform_Analysis/      Differential transcript expression and usage analyses
Gene_Analysis/         Differential gene expression analyses
Enrichment_Analysis/   Background and tested gene lists for all enrichment analyses
single-cell/           Single-cell and pseudobulk transcriptomic analyses
```

---

## Brief Methods

Long-read cDNA libraries were prepared using the Oxford Nanopore PCR-cDNA Barcoding Kit (SQK-PCB114.24) and sequenced on the PromethION platform (R10.4.1, Oxford Nanopore, FLO-PRO114). Demultiplexing and base calling was performed using Dorado (v0.9.1) SUP (super high accuracy) model with default parameters. Reference files used were using GRCm39 mouse genome and GENCODE vM38 annotation. Short-read cDNA libraries were prepared using the SMARTer Ultra Low RNA Kit for Illumina Sequencing (Takara-Clontech) and sequenced on the Illumina NovaSeq X Plus platform with at least 50 million paired-end reads per sample. Base calling and demultiplexing were performed using Illumina's bcl2fastq software.

---

## Isoform Reconstruction, Filtering, and Quantification

### IsoQuant

Long-read reads were aligned and reconstructed using [IsoQuant](https://github.com/ablab/IsoQuant).

Main outputs include:
- transcript annotations
- isoform abundance estimates
- reconstructed transcript models
- splice junction information

### SQANTI3

[SQANTI3](https://github.com/ConesaLab/SQANTI3) was used for transcript classification, quality filtering, and structural annotation.

Analyses included:
- FSM/NIC/NNC classification
- ORF prediction
- artifact filtering
- short-read splice support integration
- transcript rescue and confidence filtering

Low-confidence isoforms identified by SQANTI3 filtering were excluded from downstream analyses.

### Quantification

Isoform-level abundance was quantified from the matched short-read data using kallisto, as implemented in SQANTI3, with short reads from each sample pseudoaligned to that sample's long-read-derived transcript models. Estimated counts for isoforms retained after filtering and rescue were merged across samples into a single isoform count matrix after harmonizing novel transcript identifiers.

Gene-level counts were obtained by summing the estimated counts of all isoforms assigned to the same Ensembl gene, so that gene- and isoform-level analyses were derived from the same quantification.

---

## Analysis Summary

### Differential Transcript Expression (DTE) and Differential Gene Expression (DGE)

**DESeq2**

Differential transcript expression analyses were performed using the isoform count matrix described above. Differential gene expression analyses were performed using:

- gene-level counts summed from the isoform count matrix, for direct comparison with DTE
- gene-level counts from standard short-read alignment (STAR + featureCounts), used to assess concordance

(Threshold: adjusted p-value ≤ 0.05 and |Log2FoldChange| ≥ 1)

### Differential Transcript Usage (DTU)

**IsoformSwitchAnalyzeR + DEXSeq**

Differential transcript usage analyses were performed using IsoformSwitchAnalyzeR and DEXSeq. (Threshold: adjusted p-value ≤ 0.05 and differential isoform fraction (dIF) ≥ 0.1)

Analyses included:
- isoform switching
- isoform fraction changes
- coding sequence alterations
- domain gain/loss
- transcript emergence
- structural remodeling


### Single-Cell RNA-seq analysis

Single-cell datasets were integrated to infer cell-type-specific contributions to transcript remodeling.

Processing steps included:
- Custom Mbp-Golli / Mbp-Classic reference generation
- Cell Ranger alignment
- CellBender ambient RNA removal
- Scrublet doublet filtering
- Harmony integration
- Seurat clustering
- Pseudobulk DESeq2 analysis

### Enrichment Analyses

Three enrichment analyses were performed:

| Analysis | Method | Figure |
|---|---|---|
| GO enrichment of multi-DET genes | clusterProfiler | Figure 2C |
| GO enrichment of DET-only genes | clusterProfiler | Supplementary Figure 3D |
| Enrichment of multi-DET genes among cell-type DEGs | Logistic regression | Figure 4C |

Gene Ontology analyses used genes represented by at least two isoforms in the filtered isoform set as background. The cell-type analysis used, for each cell type and comparison, genes with at least two tested isoforms that were also tested in that cell type's pseudobulk analysis.

---

## Data Availability

**Bulk RNA Sequencing**

Raw ONT long-read and Illumina short-read datasets are available through SRA **BioProject:** PRJNA1462824. Data will become publicly available upon publication.


**Single-Cell Datasets**

Previously published datasets used in this study are from GSE291435 and GSE198582

**Enrichment Analysis Gene Lists**

The background (universe) and tested gene lists for all enrichment analyses are available in this repository under [`Enrichment_Analysis/`](Enrichment_Analysis/):

| File | Contents |
|---|---|
| `GO_background_multi_isoform_genes.csv` | Background gene list for both GO analyses |
| `GO_tested_multi_DET_genes.csv` | Multi-DET genes tested in Figure 2C |
| `GO_tested_DET_only_genes.csv` | DET-only genes tested in Supplementary Figure 3D |
| `Celltype_enrichment_background_genes.csv` | Background genes and model variables for each cell type and comparison in Figure 4C |
| `Celltype_enrichment_background_summary.csv` | Background size and gene counts per cell type and comparison |

See [`Enrichment_Analysis/README.md`](Enrichment_Analysis/README.md) for column descriptions.


## Software and Packages

Primary software used in this study includes:

- [IsoQuant](https://github.com/ablab/IsoQuant)
- [SQANTI3](https://github.com/ConesaLab/SQANTI3)
- DESeq2
- IsoformSwitchAnalyzeR
- DEXSeq
- kallisto
- clusterProfiler
- Seurat
- Harmony
- CellBender
- Scrublet
- STAR
- featureCounts
- Dorado

## Contact

**Wendy Dong**  
MD-PhD Candidate  
Washington University School of Medicine in St. Louis
wendy.dong@wustl.edu

## License

This repository is distributed under the MIT License. See `LICENSE` for details.
