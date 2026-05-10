# IsoQuant Workflow

This directory contains scripts used for long-read alignment, transcript reconstruction, and isoform quantification using [IsoQuant](https://github.com/ablab/IsoQuant).

---

# Overview

Oxford Nanopore long-read RNA sequencing reads were aligned to the mouse GRCm39 genome reference and processed using:
- IsoQuant v3.6.0
- minimap2 v2.28-r1209

Reference annotation:
- GENCODE vM38

Across all samples:
- Mean mapped reads per sample: 30.11 million
- Mean mapping rate: 91.64%

IsoQuant outputs were used for:
- transcript reconstruction
- isoform quantification
- novel isoform discovery
- downstream SQANTI3 filtering and annotation
- differential transcript expression analyses

---

# Computational Environment

IsoQuant was executed within a Docker container environment to ensure reproducibility.

Docker image used:

```bash
etycksen/isoquant:latest
```

Containerized execution included:
- IsoQuant v3.6.0
- minimap2 v2.28-r1209
- required Python dependencies

---

# Requirements

## Software

- [Docker](https://www.docker.com/)
- [IsoQuant](https://github.com/ablab/IsoQuant)


---

# Input Files

Required inputs:
- Oxford Nanopore FASTQ files
- GRCm39 genome FASTA
- GENCODE vM38 GTF annotation

---

# Running IsoQuant

The workflow was executed inside the IsoQuant Docker container.

Example Docker execution:

```bash
docker run \
-v /path/to/project:/workspace \
-w /workspace \
etycksen/isoquant:latest \
bash run_isoquant.sh
```

The script:
1. merges FASTQ files
2. compresses merged reads
3. aligns reads to the reference genome
4. reconstructs transcript models
5. generates SQANTI3-compatible outputs

---

# Important Parameters

| Parameter | Description |
|---|---|
| `-d nanopore` | Specifies Oxford Nanopore long-read data |
| `--stranded forward` | Forward-stranded library preparation |
| `--complete_genedb` | Uses complete annotation database |
| `--sqanti_output` | Generates SQANTI3-compatible outputs |
| `--check_canonical` | Checks canonical splice junctions |
| `--count_exons` | Reports exon-level counts |
| `--report_novel_unspliced true` | Retains novel unspliced transcript models |

---

# Main Output Files

Typical outputs include:

| File Type | Description |
|---|---|
| `.gtf` | Reconstructed transcript annotations |
| `.tsv` | Transcript abundance estimates |
| splice junction summaries | Canonical/noncanonical splice information |

These outputs were subsequently used for:
- SQANTI3 filtering
- differential transcript expression (DTE)
- differential transcript usage (DTU)
- novel isoform characterization

---

# Reference

IsoQuant GitHub repository:
https://github.com/ablab/IsoQuant
