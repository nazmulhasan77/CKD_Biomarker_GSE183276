# CKD Biomarker Discovery Project README

## Thesis Topic

**Identification of Candidate Biomarkers for Chronic Kidney Disease Using Single-Cell Spatial Transcriptomic Data**

Ei project-e `GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Seurat_03282022.h5ad` dataset use kore CKD candidate biomarker identify kora hoyeche.

Dataset-ti Lake et al. kidney healthy-injury atlas-er scCv3 single-cell dataset. Ei `.h5ad` file-e single-cell transcriptomic expression, cell annotation, patient annotation, disease condition, cell class/subclass, state annotation ebong reference UMAP ache. Ei file-e direct tissue spatial x-y coordinate nei, tai ei analysis mainly single-cell transcriptomic atlas metadata based. Spatial validation korte hole corresponding spatial transcriptomics file/data lagbe.

## Dataset Summary

Input file:

```text
GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Seurat_03282022.h5ad
```

AnnData size:

```text
109,741 cells x 37,080 genes
```

Condition distribution:

```text
CKD    52,314 cells
AKI    35,777 cells
Ref    21,650 cells
```

Patient distribution:

```text
CKD    15 patients
AKI    12 patients
Ref    18 patients
```

Cell class distribution:

```text
epithelial cells      80,625
immune cells          17,176
endothelial cells      8,550
stroma cells           3,383
neural cells               7
```

Major cell subclasses:

```text
PT, TAL, IMM, EC, IC, PC, CNT, DCT, VSM/P, DTL, PEC, ATL, FIB, POD, PapE, NEU
```

## Important Data Handling Decision

Ei `.h5ad` file-e `adata.X` scaled expression, tai negative value ache. Differential expression ba ML analysis-er jonno scaled matrix direct use kora uchit na.

Tai pipeline-e expression source hisebe use kora hoyeche:

```python
adata.raw.X
```

Gene symbol source:

```python
adata.raw.var["_index"]
```

Eta kora hoyeche karon `adata.raw.X` sparse log-normalized expression matrix, ja biomarker analysis-er jonno better.

## Main Question

Ei project-e duita main question answer kora hoyeche:

1. **Full dataset theke overall CKD candidate biomarker konta?**
2. **Cell type/subclass wise CKD biomarker konta?**

Tar sathe extra machine learning analysis kora hoyeche:

3. **Random Forest and KNN model diye candidate biomarker ranking kora jay kina?**

## Why Patient-Level Pseudobulk Used

Single-cell dataset-e cell onek, kintu biological replicate holo **patient**, individual cell na.

Jodi 109,741 cell-ke independent sample dhora hoy, tahole pseudo-replication hobe. Same patient-er thousands of cells model/test-e chole gele p-value ebong model accuracy artificially high hote pare.

Tai ami patient-level pseudobulk approach use korechi:

1. Same patient-er cells group kora hoy.
2. Protita patient-er jonno gene-wise average log-expression calculate kora hoy.
3. Tarpor CKD patient vs Ref patient compare kora hoy.

Bangla summary:

```text
Cell-level na, patient-level analysis kora hoyeche jate thesis-er result statistically defensible hoy.
```

## Analysis 1: Full Dataset CKD vs Reference Biomarker

Script:

```text
scripts/ckd_biomarker_pipeline.py
```

Main comparison:

```text
CKD vs Ref
```

Additional specificity check:

```text
CKD vs AKI
```

Ki kora hoyeche:

- All cells use kore patient-level pseudobulk matrix banano hoyeche.
- CKD vs Ref differential expression kora hoyeche.
- Welch t-test use kora hoyeche.
- Benjamini-Hochberg FDR correction kora hoyeche.
- Cohen's d effect size calculate kora hoyeche.
- CKD and Ref-e gene detection percentage calculate kora hoyeche.
- Final gene ranking-er jonno composite biomarker score calculate kora hoyeche.

Main output files:

```text
results/ckd_vs_ref_all_cells_de.csv
results/ckd_vs_aki_all_cells_de.csv
results/ckd_candidate_biomarkers_ranked.csv
results/ckd_upregulated_candidate_biomarkers.csv
results/ckd_downregulated_reference_loss_genes.csv
```

## Full Dataset Top CKD-Upregulated Candidate Biomarkers

Output file:

```text
results/ckd_upregulated_candidate_biomarkers.csv
```

Top CKD-upregulated candidate genes:

```text
WFDC2
ITGB8
MAP3K1
TBL1XR1
LRRK2
TNFRSF11B
TSPAN1
PIGR
PROM1
MMP7
CYBA
TNFSF10
MARCKS
SPP1
CLU
CDH6
SESN3
ITGB6
CFI
VCAM1
```

Interpretation:

Positive `diff_mean_log_expr` mane gene expression CKD patient-e Ref-er cheye beshi.

Example:

```text
WFDC2, ITGB8, MAP3K1, PROM1, MMP7, VCAM1 ei genes-gulo CKD vs Ref comparison-e upregulated candidate biomarker hisebe rank peyechhe.
```

## Full Dataset CKD-Downregulated / Reference-Loss Genes

Output file:

```text
results/ckd_downregulated_reference_loss_genes.csv
```

Top CKD-downregulated genes:

```text
PDK4
ZFAND5
CEBPD
KLF9
GADD45A
PER1
NFKBIA
HERPUD1
TSC22D3
ERRFI1
NR4A1
SLC25A33
FKBP5
RASD1
SGK1
```

Interpretation:

Negative `diff_mean_log_expr` mane gene expression CKD patient-e Ref-er cheye kom. Ederke direct CKD-up biomarker bola thik na. Ederke **reference-loss**, **CKD-associated downregulated gene**, ba **loss-of-normal-state marker** bola better.

## CKD vs AKI Specificity Check

CKD biomarker jodi AKI-teo high hoy, tahole seta general kidney injury marker hote pare, CKD-specific marker na.

Tai ami CKD vs AKI comparison-o korechi.

Output:

```text
results/ckd_vs_aki_all_cells_de.csv
results/ckd_upregulated_and_aki_specific_biomarkers.csv
```

Important result:

```text
Strict all-cell CKD-up gene jeta CKD vs Ref FDR < 0.05 and CKD vs AKI FDR < 0.05 pass kore, emon gene ei run-e paoa jayni.
```

Tai final conclusion:

```text
Upregulated genes-gulo CKD-vs-healthy candidate biomarker. CKD-specific diagnostic biomarker claim korte independent validation lagbe.
```

## Analysis 2: Cell-Type / Subclass-Wise Biomarker

Full dataset analysis charao ami cell subclass wise analysis korechi.

Used annotation:

```text
obs["subclass.l1"]
```

Included subclasses:

```text
CNT
DCT
EC
IC
IMM
PC
PEC
PT
TAL
VSM/P
```

Output:

```text
results/subclass_l1_ckd_vs_ref_top_markers.csv
```

Ki kora hoyeche:

- Protita subclass-er moddhe CKD vs Ref compare kora hoyeche.
- Patient-level pseudobulk rakha hoyeche.
- Protita subclass-er jonno top overall, CKD-up, and CKD-down marker save kora hoyeche.

Important column:

```text
subclass.l1
marker_set
gene
rank_score
adj_p_value
diff_mean_log_expr
cohen_d
```

`marker_set` column-er meaning:

```text
top_overall  = strongest subclass marker, up or down both included
up_in_CKD    = CKD-te expression beshi
down_in_CKD  = CKD-te expression kom
```

Example subclass-wise signals:

```text
PT:   ERRFI1, GADD45A, PDK4, USP2, TSC22D3 down in CKD; CLU, S100A11 up signal
TAL:  ZFAND5, USP2, CEBPD, PDK4 down in CKD; WFDC2 up signal
EC:   MT2A, KLF9, MT1E down in CKD; SOX18, MARCKS up signal
IMM:  CXCR4, SRGN, TSC22D3 down in CKD; ZFP36L1 up signal
```

Cell-type-wise result thesis-e use korle bhalo sentence:

```text
The biomarker signal was further localized to kidney cell subclasses using patient-level pseudobulk differential expression within each annotated subclass.
```

## Analysis 3: Machine Learning Biomarker Discovery

ML folder:

```text
Using ML/
```

ML scripts:

```text
Using ML/01_build_patient_pseudobulk.py
Using ML/02_random_forest_knn_biomarkers.py
```

ML README:

```text
Using ML/README.md
```

ML-e ki kora hoyeche:

1. CKD and Ref patient-level pseudobulk matrix banano hoyeche.
2. Random Forest classifier train kora hoyeche.
3. KNN classifier train kora hoyeche.
4. Random Forest feature importance calculate kora hoyeche.
5. KNN-er jonno single-gene cross-validation importance calculate kora hoyeche.
6. KNN permutation importance-o export kora hoyeche.
7. RF importance, KNN evidence, and statistical score combine kore ML biomarker rank banano hoyeche.

ML dataset:

```text
CKD vs Ref only
15 CKD patients
18 Ref patients
33 total patient profiles
29,010 non-constant genes
```

ML model performance:

```text
Random Forest CV balanced accuracy: 1.00
KNN CV balanced accuracy: about 0.94
```

Important caution:

Patient count small, tai ML result candidate biomarker hisebe use korte hobe. Clinical diagnostic biomarker claim kora jabe na without external validation.

## ML Output Files

```text
Using ML/data/patient_pseudobulk_logexpr.csv
Using ML/data/patient_metadata.csv
Using ML/results/model_cv_metrics.csv
Using ML/results/random_forest_gene_importance.csv
Using ML/results/knn_single_gene_cv_importance.csv
Using ML/results/knn_permutation_importance.csv
Using ML/results/ml_combined_candidate_biomarkers.csv
Using ML/results/ml_upregulated_candidate_biomarkers.csv
Using ML/results/ml_downregulated_reference_loss_genes.csv
```

## Top ML CKD-Upregulated Candidate Biomarkers

Output:

```text
Using ML/results/ml_upregulated_candidate_biomarkers.csv
```

Top genes:

```text
TLR1
PROM1
DYRK2
GALNT3
C9orf139
IKBKE
GRHL3
MAP3K1
RASGRP1
TLR7
PLEKHH1
SAMD12
CDH6
TBL1XR1
MLIP
VCAM1
SESN3
CXADR
ITGB8
SYT16
CLDN1
ENPP5
AP1G2
DTX4
GPX8
```

ML and DE overlap-er important genes:

```text
PROM1
MAP3K1
CDH6
TBL1XR1
VCAM1
SESN3
ITGB8
```

Ei overlap genes-gulo thesis discussion-e priority candidate hisebe highlight kora jay, karon differential expression and ML duita angle thekei support pachhe.

## Top Combined ML Genes

Output:

```text
Using ML/results/ml_combined_candidate_biomarkers.csv
```

Top combined ML genes:

```text
KLF9
ZFAND5
ELL2
PDK4
CEBPD
NFIL3
PER1
SLC25A33
TIPARP
APOLD1
CSRNP1
F2RL3
USP2
TLR1
ERRFI1
```

Note:

Ei combined list-e upregulated and downregulated both gene ache. Biomarker candidate bolte jodi CKD-te high expression marker bujhao, tahole `ml_upregulated_candidate_biomarkers.csv` use koro.

## Figures

Main figures:

```text
figures/umap_condition_l1.png
figures/umap_class.png
figures/cell_class_composition_by_condition.png
figures/volcano_ckd_vs_ref_all_cells.png
figures/top20_ckd_upregulated_candidate_biomarkers.png
figures/top20_ckd_downregulated_reference_loss_genes.png
```

ML figures:

```text
Using ML/figures/top20_random_forest_importance.png
Using ML/figures/top20_knn_single_gene_cv_importance.png
Using ML/figures/top20_knn_permutation_importance.png
Using ML/figures/top20_combined_ml_biomarkers.png
Using ML/figures/top20_ml_upregulated_biomarkers.png
```

## How To Run Full DE Pipeline

From project root:

```powershell
python scripts/ckd_biomarker_pipeline.py
```

Faster run without subclass analysis:

```powershell
python scripts/ckd_biomarker_pipeline.py --skip-subclass
```

Lower memory run:

```powershell
python scripts/ckd_biomarker_pipeline.py --chunk-size 500
```

## How To Run ML Pipeline

From project root:

```powershell
python "Using ML\01_build_patient_pseudobulk.py"
python "Using ML\02_random_forest_knn_biomarkers.py"
```

Or run both with PowerShell:

```powershell
powershell -ExecutionPolicy Bypass -File "Using ML\run_using_ml_pipeline.ps1"
```

## Important Column Meaning

`gene`

```text
Gene symbol
```

`diff_mean_log_expr`

```text
CKD mean log-expression minus Ref mean log-expression
Positive = CKD-te beshi
Negative = CKD-te kom
```

`adj_p_value`

```text
FDR-adjusted p-value
Usually < 0.05 significant dhora hoy
```

`cohen_d`

```text
Effect size
Positive = CKD up
Negative = CKD down
Absolute value joto beshi, effect toto strong
```

`rank_score`

```text
DE-based composite ranking score
```

`combined_ml_biomarker_score`

```text
Random Forest importance, KNN single-gene importance, and statistical score combined ML ranking
```

`random_forest_importance`

```text
Random Forest model-er gene importance
```

`knn_single_gene_balanced_accuracy_mean`

```text
Single gene diye KNN cross-validation-e CKD vs Ref distinguish korar performance
```

## Recommended Final Candidate Gene Strategy

Thesis-er jonno ekta balanced candidate list banate chaile ei rule follow kora jay:

1. Gene CKD vs Ref-e upregulated hote hobe.
2. `adj_p_value < 0.05` thaka bhalo.
3. Positive `cohen_d` thaka bhalo.
4. RF/KNN ML ranking-e support thaka bhalo.
5. Cell-type/subclass localization thaka bhalo.
6. Literature support thaka bhalo.

Strong candidate discussion-er jonno suggested overlap genes:

```text
PROM1
MAP3K1
CDH6
TBL1XR1
VCAM1
SESN3
ITGB8
WFDC2
MMP7
TNFRSF11B
PIGR
```

## Suggested Thesis Methods Text

```text
Single-cell transcriptomic data from GSE183276 were analyzed to identify candidate biomarkers associated with chronic kidney disease. Since individual cells from the same patient are not independent biological replicates, expression values were aggregated into patient-level pseudobulk profiles. Differential expression analysis was performed between CKD and healthy reference samples using patient-level profiles, followed by multiple-testing correction and effect-size estimation. Candidate genes were ranked using adjusted p-value, expression difference, detection frequency, and Cohen's d. Cell-subclass-specific analyses were also performed to localize biomarker signals to kidney cell populations. In addition, Random Forest and K-Nearest Neighbors classifiers were trained on patient-level pseudobulk expression profiles to provide complementary machine-learning-based feature prioritization.
```

## Suggested Thesis Result Text

```text
The all-cell patient-level analysis identified WFDC2, ITGB8, MAP3K1, TBL1XR1, LRRK2, TNFRSF11B, PROM1, MMP7, and VCAM1 as CKD-upregulated candidate biomarkers compared with healthy reference kidneys. Machine-learning-based prioritization using Random Forest and KNN supported several overlapping candidates including PROM1, MAP3K1, CDH6, TBL1XR1, VCAM1, SESN3, and ITGB8. Subclass-wise analysis further localized disease-associated transcriptional changes across nephron and non-epithelial compartments, including PT, TAL, EC, IMM, and VSM/P populations.
```

## Limitations

Important limitations:

```text
1. Patient count small for high-dimensional ML.
2. CKD group includes DKD and hypertensive CKD conditions.
3. The h5ad file analyzed here is scCv3 single-cell data, not direct spatial coordinate matrix.
4. Strict CKD-vs-AKI-specific upregulated biomarkers were not detected at all-cell level.
5. Candidate biomarkers require external validation.
```

Validation suggestions:

```text
Independent CKD cohort validation
Bulk RNA-seq validation
External single-cell validation
Protein-level validation
ROC/AUC analysis in independent samples
Literature-based biological interpretation
```

## Final Short Summary

Ei project-e ami:

```text
1. Full dataset theke CKD vs Ref biomarker identify korechi.
2. Cell subclass wise biomarker identify korechi.
3. CKD vs AKI specificity check korechi.
4. Random Forest and KNN diye ML-based biomarker ranking korechi.
5. Sob output CSV and figure folder-e save korechi.
6. Patient-level pseudobulk use kore pseudo-replication avoid korechi.
```

Main final files:

```text
results/ckd_upregulated_candidate_biomarkers.csv
results/subclass_l1_ckd_vs_ref_top_markers.csv
Using ML/results/ml_upregulated_candidate_biomarkers.csv
Using ML/results/ml_combined_candidate_biomarkers.csv
```

