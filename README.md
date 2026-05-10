# Isoform Remodeling Shapes Peripheral Nerve Response to Injury

Companion repository for the submitted manuscript **Dong W et al. Isoform Remodeling Shapes Peripheral Nerve Response to Injury. (2026)**

This repository contains the computational workflows, analysis scripts, and supporting resources used to generate and analyze long-read and short-read transcriptomic datasets from mouse sciatic nerve crush injury.

The study integrates:

- Oxford Nanopore long-read RNA sequencing
- Illumina short-read RNA sequencing
- Single-cell and pseudobulk transcriptomic analyses

to characterize isoform-level remodeling during Wallerian degeneration.

An interactive companion browser for exploring transcript- and gene-level results is available at **IsoNerve** (https://isonerve.pages.dev)

If you use this repository or analysis framework, please cite:

> Dong W et al. *Isoform Remodeling Shapes Peripheral Nerve Response to Injury.* (2026)

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
Figures/            Manuscript figures and supplementary figures
IsoQuant/           Long-read alignment and transcript reconstruction workflows
SQANTI3/            Isoform QC, filtering, and annotation workflows
Isoform_Analysis/   Differential transcript expression and usage analyses
Short_Read/         Differential gene expression analysis
single-cell/        Single-cell and pseudobulk transcriptomic analyses
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

---

## Analysis Summary

### Differential Transcript Expression (DTE) and Differential Gene Expression (DGE)

**DESeq2**

Differential transcript expression analyses were performed using isoform-level count matrices derived from long-read RNA sequencing. Differential gene expression analyses were performed using gene-level count matrices derived from short-read RNA sequencing. (Threshold: adjusted p-value ≤ 0.05 and Log2FoldChange ≥ 1)

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

---

## Data Availability

**Bulk RNA Sequencing**

Raw ONT long-read and Illumina short-read datasets are available through SRA **BioProject:** PRJNA1462824. Data will become publicly available upon publication.


**Single-Cell Datasets**

Previously published datasets used in this study are from GSE291435 and GSE198582


## Software and Packages

Primary software used in this study includes:

- [IsoQuant](https://github.com/ablab/IsoQuant)
- [SQANTI3](https://github.com/ConesaLab/SQANTI3)
- DESeq2
- IsoformSwitchAnalyzeR
- DEXSeq
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
