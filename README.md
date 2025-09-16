# Long-Read RNA-seq Analysis of Mice Crush Nerve Injury

This repository contains the code and scripts used for the analysis of long-read RNA-seq data in a mouse peripheral nerve crush injury model (C0, C3, C7 days post-injury). The analysis integrates Oxford Nanopore long-read RNA sequencing with short-read RNA-seq to identify, quality-control, and quantify isoforms, followed by downstream differential expression and isoform-level analyses.


---

## Analysis Workflow

1. **IsoQuant (long-read alignment & isoform reconstruction)**  
   - Tool: [IsoQuant](https://github.com/ablab/IsoQuant)  
   - Input: Nanopore `fastq.gz` files, GRCm39 genome reference, GENCODE annotation  
   - Output: Sample-specific isoform GTFs and counts  

2. **SQANTI3 (QC & filtering)**  
   - Tool: [SQANTI3](https://github.com/ConesaLab/SQANTI3)  
   - Functions:  
     - Classification of isoforms (FSM, NIC, NNC, antisense, etc.)  
     - Quality filtering (junction support, ORF validation, artifact removal)  
     - Integration with short-read RNA-seq via *kallisto* quantification  

3. **Downstream R Analysis**  
   - **DESeq2**: Differential gene and isoform expression across conditions (C0, C3, C7)  
   - **Novel isoform analysis**: Identification and characterization of injury-induced novel transcripts  
   - **DTU (Differential Transcript Usage)**: Performed with [IsoformSwitchAnalyzeR](https://bioconductor.org/packages/release/bioc/html/IsoformSwitchAnalyzeR.html) to detect isoform switching events and functional consequences  

---

## Requirements

- **IsoQuant** (docker(etycksen/isoquant:latest))  
- **SQANTI3** (Python 3 + dependencies)  
- **Kallisto** (≥ v0.46)  
- **R (≥ 4.2.0)** with the following packages:  
  - `DESeq2`  
  - `IsoformSwitchAnalyzeR`  
  - `dplyr`, `ggplot2`, `tidyr`  
  - `Biostrings`, `GenomicFeatures`, `tximport`  


---

## Citation

If you use this repository or workflow, please cite our manuscript (in preparation):  

*Dong W, et al. Long-read RNA sequencing reveals novel isoform dynamics in peripheral nerve injury.*  

---

## Contact

For questions about the code or analysis, please contact:  
**Wendy Dong** – MD-PhD Student, Washington University in St. Louis


