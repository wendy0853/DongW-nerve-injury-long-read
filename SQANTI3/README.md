# SQANTI3 Workflow

This directory contains scripts used for isoform annotation, quality control, filtering, and rescue using [SQANTI3](https://github.com/ConesaLab/SQANTI3).

---

## Overview

Isoform annotation and quality control were performed using SQANTI3 v5.5.3 with default parameters.

SQANTI3 was used to:

- classify reconstructed isoforms
- assess splice junction support
- annotate FSM, NIC, NNC, and other transcript classes
- identify low-confidence isoforms
- perform rules-based filtering and rescue
- generate filtered transcript models for downstream analysis

Short-read RNA-seq alignments from matched samples were incorporated to provide splice junction support. Reference CAGE peak, poly(A) motif, and poly(A) site files supplied with SQANTI3 were used for transcript classification.

Following SQANTI3 default rules-based filtering and rescue, low-confidence isoforms were excluded from downstream analyses.

---

## Computational Environment

SQANTI3 was executed in a Docker container.

Docker image used:

```bash
jinlab/sqanti3:vs1
```

Example Docker execution:

```bash
docker run \
-v /path/to/project:/workspace \
-v /path/to/reference:/reference \
-w /workspace \
jinlab/sqanti3:vs1 \
bash run_sqanti3_qc.sh sample_name
```

---

## Workflow

The SQANTI3 workflow was performed in three main steps:

1. `sqanti3_qc.py`
2. `sqanti3_filter.py`
3. `sqanti3_rescue.py rules --mode automatic`

---

## Step 1: SQANTI3 QC

Run:

```bash
bash run_sqanti3_qc.sh sample_name
```

This step takes IsoQuant transcript models as input and generates SQANTI3 classification and corrected transcript files.

Main outputs include:

| Output | Description |
|---|---|
| `*_classification.txt` | SQANTI3 transcript classification file |
| `*_corrected.gtf` | Corrected transcript annotation |
| `*_corrected.fasta` | Corrected transcript sequences |
| `*_junctions.txt` | Splice junction annotation |

---

## Step 2: SQANTI3 Rules-Based Filtering

Run:

```bash
bash run_sqanti3_filter.sh sample_name
```

This step applies SQANTI3 rules-based filtering to remove low-confidence isoforms.

Main outputs include:

| Output | Description |
|---|---|
| `*_filtered_RulesFilter_result_classification.txt` | Filtered classification file |
| `*_filtered_RulesFilter_result.gtf` | Filtered transcript annotation |

---

## Step 3: SQANTI3 Automatic Rescue

Run:

```bash
bash run_sqanti3_rescue_rules_automatic.sh sample_name
```

This step applies SQANTI3 automatic rescue rules to recover supported transcripts after filtering.

Main outputs include:

| Output | Description |
|---|---|
| `*_rescue_rescue_rules_classification.txt` | Rescued transcript classification |
| rescued GTF outputs | Final transcript models after filtering and rescue |

---

## Required Inputs

| Input | Description |
|---|---|
| IsoQuant transcript GTF | Transcript models from IsoQuant |
| Reference GTF | GENCODE vM38 annotation |
| Reference FASTA | GRCm39 genome FASTA |
| Short-read BAM/FOFN | Matched short-read RNA-seq alignments |
| Poly(A) motif file | SQANTI3-supplied mouse/human poly(A) motif list |
| CAGE peak file | SQANTI3-supplied mouse CAGE/refTSS annotation |

---

## Reference

SQANTI3 GitHub repository:  
https://github.com/ConesaLab/SQANTI3
