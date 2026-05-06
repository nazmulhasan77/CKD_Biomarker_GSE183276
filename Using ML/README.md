# Using ML: Random Forest and KNN Biomarker Discovery

This folder contains a machine-learning workflow for candidate CKD biomarker discovery using the GSE183276 kidney atlas.

## Why Patient-Level ML?

The dataset has many cells, but only a limited number of patients. For biomarker discovery, the patient is the biological replicate. Therefore, this workflow first aggregates cells into patient-level pseudobulk expression profiles, then trains ML models on patients.

Bangla note: cell-level ML korle same patient-er thousands of cells model-e chole jabe, tai performance artificially high dekhabe. Patient-level pseudobulk use kora thesis-er jonno safer.

## Models Used

1. **Random Forest Classifier**
   - Trains CKD vs Ref classification.
   - Uses native `feature_importances_` to rank genes.

2. **K-Nearest Neighbors Classifier**
   - Trains CKD vs Ref classification.
   - KNN does not have built-in gene coefficients.
   - The workflow estimates KNN evidence in two ways:
     - single-gene KNN cross-validation performance
     - permutation importance after fitting the final KNN model

## Files

- `01_build_patient_pseudobulk.py`
  - Reads the `.h5ad` file.
  - Uses `adata.raw.X`.
  - Builds patient-level average log-expression.
  - Saves `data/patient_pseudobulk_logexpr.csv`.

- `02_random_forest_knn_biomarkers.py`
  - Trains Random Forest and KNN.
  - Performs cross-validation.
  - Computes RF importance, KNN permutation importance, and combined ML biomarker score.

- `run_using_ml_pipeline.ps1`
  - Runs both scripts in order.

## How To Run

From the main dataset folder:

```powershell
python "Using ML\01_build_patient_pseudobulk.py"
python "Using ML\02_random_forest_knn_biomarkers.py"
```

Or run both together:

```powershell
powershell -ExecutionPolicy Bypass -File "Using ML\run_using_ml_pipeline.ps1"
```

## Main Outputs

- `results/random_forest_gene_importance.csv`
- `results/knn_permutation_importance.csv`
- `results/knn_single_gene_cv_importance.csv`
- `results/ml_combined_candidate_biomarkers.csv`
- `results/ml_upregulated_candidate_biomarkers.csv`
- `results/ml_downregulated_reference_loss_genes.csv`
- `results/model_cv_metrics.csv`
- `figures/top20_random_forest_importance.png`
- `figures/top20_knn_permutation_importance.png`
- `figures/top20_knn_single_gene_cv_importance.png`
- `figures/top20_combined_ml_biomarkers.png`
- `figures/top20_ml_upregulated_biomarkers.png`

## Recommended Thesis Interpretation

Use `ml_upregulated_candidate_biomarkers.csv` for CKD-upregulated biomarker candidates. Use `ml_combined_candidate_biomarkers.csv` when discussing the complete ML ranking, including genes that are lower in CKD.

Suggested wording:

> Patient-level pseudobulk expression profiles were generated from single-cell transcriptomic data to avoid pseudo-replication. Random Forest and K-Nearest Neighbors classifiers were trained to distinguish CKD from healthy reference samples. Random Forest feature importance, KNN single-gene cross-validation importance, and patient-level differential expression statistics were combined to prioritize candidate CKD biomarkers.

## Caution

The patient count is small for high-dimensional ML. Treat these as candidate biomarkers. Independent validation is needed before claiming diagnostic biomarkers.
