# Manuscript Figure Generation

This directory contains scripts used to generate manuscript main and supplementary figures.

Most figure panels were generated using outputs from:

- `Isoform_Analysis/`
- `single-cell/`
- `Short_Read/`

---

## Figure Scripts

| Script | Description |
|---|---|
| Figure_1.R | Long-read RNA-seq reveals high-quality transcriptome profiling and isoform remodeling after peripheral nerve injury. |
| Figure_2.R | Peripheral nerve injury induces widespread isoform remodeling with distinct modes of transcript emergence and functional consequences. |
| Figure_3.R | Isoform-level analysis reveals differential transcript expression independent of gene-level changes after peripheral nerve injury. |
| Figure_4.R | Heterogeneous isoform regulation across nerve-resident cell types. |
| Figure_5.R | Characterization of novel isoforms and validation of Lama4-novel-6 after peripheral nerve injury. |
| Figure_6.R | Differential transcript usage (DTU) analysis identifies Postn isoform switching, exemplifying dynamic extracellular matrix remodeling during nerve repair. |
| Supplementary_Figure_*.R | Supplementary figure generation scripts |

---

## Shared Plotting Utilities

Several figures use shared helper plotting scripts for transcript visualization and isoform-level quantification.

| Script | Purpose |
|---|---|
| Plotting_isoform_expression.R | Isoform expression line plot across conditions |
| Plotting_isoform_proportion.R | Isoform proportion bar plot |
| Plotting_isoform_trackplot.R | Transcript structure and exon architecture visualization |

These scripts were reused across multiple figures and panels to maintain consistent formatting and visualization styles throughout the manuscript.

---

## Notes

Figure scripts are organized by final manuscript figure number rather than analysis workflow.

Many scripts assume that upstream analyses from:

- `Isoform_Analysis/`
- `single-cell/`
- `Short_Read/`

have already been completed.

Paths marked with:

```r
# <-- MODIFY HERE
```

should be updated prior to execution.
