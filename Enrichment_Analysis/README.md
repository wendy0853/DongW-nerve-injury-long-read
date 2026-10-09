# Enrichment Analysis Gene Lists

This directory contains the background (universe) and tested gene lists used in the enrichment analyses of this study, so that each analysis can be reproduced.

---

## Gene Ontology Enrichment

Used for:

- GO enrichment of multi-DET genes (Figure 2C)
- GO enrichment of DET-only genes (Supplementary Figure 3D)

Both analyses use the same background: genes represented by at least two isoforms in the filtered isoform set tested for differential transcript expression. Genes represented by a single isoform were excluded.

| File | Description |
|---|---|
| `GO_background_multi_isoform_genes.csv` | Background gene list shared by both GO analyses |
| `GO_tested_multi_DET_genes.csv` | Genes with two or more differentially expressed transcripts, by comparison |
| `GO_tested_DET_only_genes.csv` | Genes with differential transcript expression but no gene-level differential expression, by comparison |

Columns:

| Column | Description |
|---|---|
| `comparison` | C3 vs. C0 or C7 vs. C0 (tested gene lists only) |
| `gene_symbol` | Gene symbol |
| `n_isoforms` | Number of isoforms in the filtered isoform set (background only) |
| `n_DETs` | Number of differentially expressed transcripts (multi-DET list only) |
| `entrez_id` | Entrez gene ID used by clusterProfiler; empty if the symbol could not be mapped |

Genes without an Entrez ID were not included in the GO test.

---

## Cell-Type Enrichment of Multi-DET Genes

Used for:

- Multi-DET enrichment across cell types (Figure 4C)

For each cell type and comparison, the background is the set of genes with at least two tested isoforms in the bulk isoform analysis that were also tested in that cell type's pseudobulk differential expression analysis.

| File | Description |
|---|---|
| `Celltype_enrichment_background_genes.csv` | One row per gene, cell type and comparison, with the variables used in the logistic regression |
| `Celltype_enrichment_background_summary.csv` | Background size and gene counts per cell type and comparison |

Columns in `Celltype_enrichment_background_genes.csv`:

| Column | Description |
|---|---|
| `comparison` | C3_vs_C0 or C7_vs_C0 |
| `cell_type` | Cell type |
| `gene_symbol` | Gene symbol |
| `n_iso` | Number of tested isoforms for the gene |
| `multiDET` | 1 if the gene has two or more differentially expressed transcripts, otherwise 0 |
| `is_DEG` | 1 if the gene is differentially expressed in that cell type's pseudobulk analysis, otherwise 0 |
| `log_expr` | log2 bulk gene expression (sum of isoform baseMean + 1) |
| `log_expr_ct` | log2 pseudobulk expression in that cell type (baseMean + 1) |

The enrichment model was:

```text
is_DEG ~ multiDET + log_expr + log_expr_ct + log2(n_iso)
```

---

## Related Scripts

| Script | Analysis |
|---|---|
| `Figures/Figure_2.R` | GO enrichment of multi-DET genes |
| `Figures/Supplementary_Figure_3.R` | GO enrichment of DET-only genes |
| `single-cell/Multi_DET_Celltype_Enrichment_Analysis.R` | Cell-type enrichment of multi-DET genes |
