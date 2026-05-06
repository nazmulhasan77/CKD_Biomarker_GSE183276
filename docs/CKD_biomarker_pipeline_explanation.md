# CKD Candidate Biomarker Pipeline: GSE183276 scCv3 Kidney Atlas

## Goal

Thesis topic: **Identification of Candidate Biomarkers for Chronic Kidney Disease Using Single-Cell Spatial Transcriptomic Data**.

Dataset used here:

`GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Seurat_03282022.h5ad`

The AnnData object contains 109,741 cells and 37,080 genes from kidney samples labelled as CKD, AKI, and healthy reference (`Ref`). The primary comparison is:

**CKD vs Ref**

AKI is used as an additional specificity comparison:

**CKD vs AKI**

Important note: in this `.h5ad`, `adata.X` contains scaled values and includes negative numbers. Therefore, the pipeline uses `adata.raw.X`, which contains sparse log-normalized expression.

## Why Patient-Level Pseudobulk?

Single-cell data has many cells, but the biological replicate is the **patient**, not each cell. If every cell is treated as an independent sample, p-values can become overconfident because thousands of cells from the same patient are correlated.

So the pipeline:

1. Groups cells by patient and condition.
2. Computes average log-normalized expression for each gene per patient.
3. Runs CKD vs Ref differential testing across patient-level profiles.
4. Ranks genes using effect size, adjusted p-value, expression prevalence, and CKD-vs-AKI specificity.

Bangla summary: ekhane protita cell ke independent sample dhora hoy nai. Patient-level pseudobulk use kora hoyeche, tai result thesis-er jonno statistically beshi defensible.

## Pipeline Steps

### 1. Load Data

The script opens the `.h5ad` in backed mode:

```python
adata = ad.read_h5ad(input_file, backed="r")
```

The matrix selection logic is:

```python
matrix = adata.raw.X
genes = adata.raw.var["_index"]
```

### 2. Metadata Summary

The code exports:

- `patient_summary.csv`
- `patient_cell_composition.csv`
- `patient_subclass_counts.csv`

These are useful for describing cohort composition and checking whether CKD samples contain different proportions of epithelial, immune, endothelial, and stromal cells.

### 3. All-Cell CKD vs Ref Biomarker Ranking

For each gene, the pipeline calculates:

- mean expression in CKD patients
- mean expression in Ref patients
- CKD minus Ref expression difference
- Welch t-test p-value
- Benjamini-Hochberg FDR
- Cohen's d effect size
- percentage of cells expressing the gene in CKD and Ref
- composite biomarker rank score

Main output:

`results/ckd_vs_ref_all_cells_de.csv`

### 4. CKD vs AKI Specificity

A gene may be a general kidney injury marker, not a CKD-specific marker. To reduce this problem, the pipeline also compares CKD vs AKI.

Main output:

`results/ckd_vs_aki_all_cells_de.csv`

The final candidate biomarker table combines CKD-vs-Ref strength with CKD-vs-AKI specificity:

`results/ckd_candidate_biomarkers_ranked.csv`

### 5. Cell-Subclass-Specific Markers

CKD biomarkers can be cell-type-specific. The pipeline repeats pseudobulk analysis inside `subclass.l1`, such as:

- PT
- TAL
- DCT
- CNT
- IC
- PC
- IMM
- EC
- VSM/P
- FIB
- POD

Main output:

`results/subclass_l1_ckd_vs_ref_top_markers.csv`

This table helps answer: **which kidney cell population is driving the biomarker signal?**

### 6. Elastic-Net Feature Importance

The script also fits a patient-level elastic-net logistic regression model using the top differential genes. This is not the primary statistical evidence because the number of patients is small, but it gives an additional machine-learning feature-importance view.

Main output:

`results/elasticnet_gene_importance.csv`

## How To Run

From the dataset folder:

```powershell
python scripts/ckd_biomarker_pipeline.py
```

Faster first run without subclass analysis:

```powershell
python scripts/ckd_biomarker_pipeline.py --skip-subclass
```

Lower-memory run:

```powershell
python scripts/ckd_biomarker_pipeline.py --chunk-size 500
```

## How To Interpret Top Genes

Prioritize genes with:

- high `combined_biomarker_score`
- `adj_p_value < 0.05`
- positive `diff_mean_log_expr` for CKD-upregulated biomarkers
- high `cohen_d`
- higher `pct_expr_case` than `pct_expr_control`
- plausible kidney disease biology from literature
- strong signal in relevant cell classes/subclasses

Example interpretation sentence:

> The candidate gene X showed significantly higher patient-level expression in CKD compared with reference kidneys, with positive effect size, higher CKD detection frequency, and enrichment in proximal tubule/injury-associated epithelial states, suggesting its potential role as a CKD-associated biomarker.

## Thesis Caution

This analysis identifies **candidate biomarkers**, not clinically validated diagnostic biomarkers. Final validation should include:

- independent CKD cohort validation
- external bulk RNA-seq or scRNA-seq validation
- protein-level evidence if available
- ROC/AUC analysis in independent patient samples
- biological interpretation using CKD literature

