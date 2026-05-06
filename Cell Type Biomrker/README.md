# কোষ-ধরনভিত্তিক CKD বায়োমার্কার বিশ্লেষণ

## উদ্দেশ্য

এই ফোল্ডারে বড় চারটি কোষ শ্রেণির জন্য আলাদা CKD বনাম Ref বায়োমার্কার বিশ্লেষণ করা হয়েছে।

বিশ্লেষিত কোষ শ্রেণি:

```text
epithelial cells
immune cells
endothelial cells
stroma cells
```

প্রতিটি কোষ শ্রেণির জন্য তিন ধরনের evidence ব্যবহার করা হয়েছে:

১. Traditional differential expression

২. Random Forest feature importance

৩. KNN single-gene cross-validation importance

এরপর তিন পদ্ধতিতে সমর্থিত জিনগুলোকে consensus biomarker হিসেবে ranking করা হয়েছে।

## কেন আলাদা cell type analysis করা হয়েছে

সম্পূর্ণ dataset থেকে overall CKD biomarker বের করলে বোঝা যায় কোন জিন CKD অবস্থায় পরিবর্তিত। কিন্তু সেই জিন কোন ধরনের কিডনি কোষে বেশি পরিবর্তিত হচ্ছে, সেটা বোঝার জন্য cell type-wise analysis দরকার।

এই বিশ্লেষণের মাধ্যমে দেখা যায়:

```text
কোন biomarker epithelial cell-এ বেশি গুরুত্বপূর্ণ
কোন biomarker immune cell-এ বেশি গুরুত্বপূর্ণ
কোন biomarker endothelial cell-এ বেশি গুরুত্বপূর্ণ
কোন biomarker stroma cell-এ বেশি গুরুত্বপূর্ণ
```

## ব্যবহৃত পদ্ধতি

প্রতিটি কোষ শ্রেণির জন্য একই পদ্ধতি অনুসরণ করা হয়েছে।

### ধাপ ১: কোষ নির্বাচন

`obs["class"]` annotation ব্যবহার করে নির্দিষ্ট cell class নির্বাচন করা হয়েছে।

### ধাপ ২: রোগী-স্তরের ছদ্ম-বাল্ক তৈরি

একই রোগীর একই cell class-এর কোষগুলো একত্র করে প্রতিটি রোগীর জন্য gene-wise average log-expression গণনা করা হয়েছে।

এটি করা হয়েছে কারণ biological replicate হলো রোগী, individual cell নয়।

### ধাপ ৩: Traditional differential expression

প্রতিটি cell class-এর মধ্যে CKD বনাম Ref তুলনা করা হয়েছে।

ব্যবহৃত পরিসংখ্যান:

```text
Welch t-test
Benjamini-Hochberg FDR correction
Cohen's d effect size
expression prevalence
```

### ধাপ ৪: Random Forest

প্রতিটি cell class-এর রোগী-স্তরের profile ব্যবহার করে Random Forest classifier চালানো হয়েছে।

Random Forest থেকে gene feature importance নেওয়া হয়েছে।

### ধাপ ৫: KNN

KNN classifier ব্যবহার করে CKD বনাম Ref আলাদা করার ক্ষমতা পরীক্ষা করা হয়েছে।

KNN-এর নিজস্ব gene coefficient নেই, তাই single-gene KNN cross-validation importance গণনা করা হয়েছে।

### ধাপ ৬: Consensus biomarker

Consensus biomarker বলতে এখানে বোঝানো হয়েছে:

```text
জিনটি traditional significant CKD-up DE table-এ আছে
জিনটি ML upregulated table-এ আছে
জিনটি Random Forest selected importance table-এ আছে
```

Traditional side-এ শুধু CKD-up হওয়া যথেষ্ট ধরা হয়নি; `adj_p_value < 0.05` ব্যবহার করা হয়েছে।

## ব্যবহৃত কোড

প্রধান স্ক্রিপ্ট:

```text
Cell Type Biomrker/cell_type_biomarker_pipeline.py
```

চালানোর কমান্ড:

```powershell
python "Cell Type Biomrker\cell_type_biomarker_pipeline.py"
```

এই run-এ ব্যবহৃত কমান্ড:

```powershell
python "Cell Type Biomrker\cell_type_biomarker_pipeline.py" --top-n-filter 300 --rf-trees 1000 --permutation-repeats 30 --chunk-size 2000
```

## রোগী ও কোষ সংখ্যা

`min-cells-per-patient = 10` filter ব্যবহার করা হয়েছে।

সারাংশ:

```text
epithelial cells:   33 patient profiles, 56,451 cells
immune cells:       32 patient profiles,  8,349 cells
endothelial cells:  32 patient profiles,  6,634 cells
stroma cells:       32 patient profiles,  2,513 cells
```

## প্রধান summary file

```text
Cell Type Biomrker/results/cell_type_biomarker_summary.csv
```

এখানে প্রতিটি cell class-এর patient সংখ্যা, cell সংখ্যা, traditional significant up gene সংখ্যা, ML up gene সংখ্যা এবং consensus gene সংখ্যা আছে।

## কোষ-ধরনভিত্তিক শীর্ষ consensus biomarker

### Epithelial cells

শীর্ষ consensus genes:

```text
ORC5
TBL1XR1
WFDC2
MAP3K1
HOXA3
DYRK2
PROM1
TNFAIP8
ITGB8
SAMD12
```

Consensus gene সংখ্যা:

```text
118
```

প্রধান result file:

```text
Cell Type Biomrker/results/epithelial_cells_consensus_traditional_ml_rf.csv
```

### Immune cells

শীর্ষ consensus genes:

```text
TLR7
ZFP36L1
CD180
TNF
MAP3K1
MAPK8IP3
ZMIZ2
SP1
MYO1C
TP53
```

Consensus gene সংখ্যা:

```text
16
```

প্রধান result file:

```text
Cell Type Biomrker/results/immune_cells_consensus_traditional_ml_rf.csv
```

### Endothelial cells

শীর্ষ consensus genes:

```text
ZNF503
SOX18
LPAR6
FOXP1
MARCKS
SFT2D2
RUNX1T1
VWF
CALHM2
SOX4
```

Consensus gene সংখ্যা:

```text
58
```

প্রধান result file:

```text
Cell Type Biomrker/results/endothelial_cells_consensus_traditional_ml_rf.csv
```

### Stroma cells

শীর্ষ consensus genes:

```text
NEDD4L
NR2F2
C4orf48
RUNX1T1
PROM1
EPS8L2
CDH6
SOX4
FLRT3
MUC1
```

Consensus gene সংখ্যা:

```text
38
```

প্রধান result file:

```text
Cell Type Biomrker/results/stroma_cells_consensus_traditional_ml_rf.csv
```

## প্রতিটি cell type-এর output file pattern

প্রতিটি cell class-এর জন্য একই ধরনের ফাইল তৈরি হয়েছে।

উদাহরণ, epithelial cells-এর জন্য:

```text
epithelial_cells_traditional_de.csv
epithelial_cells_traditional_up_biomarkers.csv
epithelial_cells_traditional_significant_up_biomarkers.csv
epithelial_cells_traditional_down_genes.csv
epithelial_cells_random_forest_importance.csv
epithelial_cells_knn_single_gene_importance.csv
epithelial_cells_knn_permutation_importance.csv
epithelial_cells_ml_combined_biomarkers.csv
epithelial_cells_ml_upregulated_biomarkers.csv
epithelial_cells_consensus_traditional_ml_rf.csv
epithelial_cells_model_cv_metrics.csv
```

অন্য cell class-এর জন্য একই pattern:

```text
immune_cells_...
endothelial_cells_...
stroma_cells_...
```

## তৈরি করা plot

প্রতিটি cell class-এর জন্য নিচের plot তৈরি হয়েছে:

```text
traditional_volcano
traditional_top20_up
random_forest_top20
knn_single_gene_top20
ml_top20_up
consensus_score_comparison
consensus_score_heatmap
consensus_overlap_counts
model_cv_balanced_accuracy
```

উদাহরণ:

```text
Cell Type Biomrker/figures/epithelial_cells_traditional_volcano.png
Cell Type Biomrker/figures/epithelial_cells_random_forest_top20.png
Cell Type Biomrker/figures/epithelial_cells_knn_single_gene_top20.png
Cell Type Biomrker/figures/epithelial_cells_consensus_score_comparison.png
```

সমন্বিত plot:

```text
Cell Type Biomrker/figures/cell_type_consensus_counts.png
Cell Type Biomrker/figures/cell_type_top_consensus_heatmap.png
```

## গুরুত্বপূর্ণ কলামের অর্থ

`gene`

```text
জিনের নাম
```

`diff_mean_log_expr`

```text
CKD mean expression minus Ref mean expression
ধনাত্মক মান = CKD-তে বেশি প্রকাশ
ঋণাত্মক মান = CKD-তে কম প্রকাশ
```

`adj_p_value`

```text
FDR-corrected p-value
0.05-এর কম হলে traditional significant ধরা হয়েছে
```

`cohen_d`

```text
effect size
ধনাত্মক মান = CKD-তে বেশি প্রকাশ
```

`rank_score`

```text
traditional differential expression ranking score
```

`random_forest_importance`

```text
Random Forest মডেলের feature importance
```

`knn_single_gene_balanced_accuracy_mean`

```text
শুধু একটি gene ব্যবহার করে KNN CKD বনাম Ref আলাদা করতে কতটা সক্ষম
```

`consensus_score`

```text
traditional score, ML score এবং Random Forest score মিলিয়ে তৈরি consensus score
```

## থিসিসে ব্যবহারযোগ্য সারাংশ

```text
কোষ-ধরনভিত্তিক বিশ্লেষণে epithelial, immune, endothelial এবং stroma cell population-এর মধ্যে CKD বনাম Ref তুলনা করা হয়েছে। প্রতিটি cell class-এর জন্য রোগী-স্তরের ছদ্ম-বাল্ক expression profile তৈরি করে traditional differential expression, Random Forest feature importance এবং KNN single-gene cross-validation importance ব্যবহার করা হয়েছে। তিন পদ্ধতিতে সমর্থিত জিনগুলোকে consensus candidate biomarker হিসেবে অগ্রাধিকার দেওয়া হয়েছে।
```

## গুরুত্বপূর্ণ ব্যাখ্যা

এই বিশ্লেষণে cell-type-specific biomarkers পাওয়া গেছে। তবে এগুলো এখনো candidate biomarker। চূড়ান্ত clinical biomarker দাবি করার আগে স্বাধীন dataset, protein-level validation এবং biological literature support দরকার।

## চিত্র গ্যালারি

এই ফোল্ডারে তৈরি করা সব cell-type plot নিচে যুক্ত করা হলো।

### সামগ্রিক cell-type plot

#### Consensus gene count

![Cell type consensus count](figures/cell_type_consensus_counts.png)

#### Top consensus heatmap

![Cell type top consensus heatmap](figures/cell_type_top_consensus_heatmap.png)

## Epithelial Cells

#### Traditional volcano plot

![Epithelial traditional volcano](figures/epithelial_cells_traditional_volcano.png)

#### Traditional top CKD-up biomarkers

![Epithelial traditional top up](figures/epithelial_cells_traditional_top20_up.png)

#### Random Forest top genes

![Epithelial Random Forest](figures/epithelial_cells_random_forest_top20.png)

#### KNN top genes

![Epithelial KNN](figures/epithelial_cells_knn_single_gene_top20.png)

#### ML top CKD-up biomarkers

![Epithelial ML up](figures/epithelial_cells_ml_top20_up.png)

#### Consensus score comparison

![Epithelial consensus comparison](figures/epithelial_cells_consensus_score_comparison.png)

#### Consensus heatmap

![Epithelial consensus heatmap](figures/epithelial_cells_consensus_score_heatmap.png)

#### Consensus overlap counts

![Epithelial consensus overlap](figures/epithelial_cells_consensus_overlap_counts.png)

#### Model cross-validation balanced accuracy

![Epithelial model CV](figures/epithelial_cells_model_cv_balanced_accuracy.png)

## Immune Cells

#### Traditional volcano plot

![Immune traditional volcano](figures/immune_cells_traditional_volcano.png)

#### Traditional top CKD-up biomarkers

![Immune traditional top up](figures/immune_cells_traditional_top20_up.png)

#### Random Forest top genes

![Immune Random Forest](figures/immune_cells_random_forest_top20.png)

#### KNN top genes

![Immune KNN](figures/immune_cells_knn_single_gene_top20.png)

#### ML top CKD-up biomarkers

![Immune ML up](figures/immune_cells_ml_top20_up.png)

#### Consensus score comparison

![Immune consensus comparison](figures/immune_cells_consensus_score_comparison.png)

#### Consensus heatmap

![Immune consensus heatmap](figures/immune_cells_consensus_score_heatmap.png)

#### Consensus overlap counts

![Immune consensus overlap](figures/immune_cells_consensus_overlap_counts.png)

#### Model cross-validation balanced accuracy

![Immune model CV](figures/immune_cells_model_cv_balanced_accuracy.png)

## Endothelial Cells

#### Traditional volcano plot

![Endothelial traditional volcano](figures/endothelial_cells_traditional_volcano.png)

#### Traditional top CKD-up biomarkers

![Endothelial traditional top up](figures/endothelial_cells_traditional_top20_up.png)

#### Random Forest top genes

![Endothelial Random Forest](figures/endothelial_cells_random_forest_top20.png)

#### KNN top genes

![Endothelial KNN](figures/endothelial_cells_knn_single_gene_top20.png)

#### ML top CKD-up biomarkers

![Endothelial ML up](figures/endothelial_cells_ml_top20_up.png)

#### Consensus score comparison

![Endothelial consensus comparison](figures/endothelial_cells_consensus_score_comparison.png)

#### Consensus heatmap

![Endothelial consensus heatmap](figures/endothelial_cells_consensus_score_heatmap.png)

#### Consensus overlap counts

![Endothelial consensus overlap](figures/endothelial_cells_consensus_overlap_counts.png)

#### Model cross-validation balanced accuracy

![Endothelial model CV](figures/endothelial_cells_model_cv_balanced_accuracy.png)

## Stroma Cells

#### Traditional volcano plot

![Stroma traditional volcano](figures/stroma_cells_traditional_volcano.png)

#### Traditional top CKD-up biomarkers

![Stroma traditional top up](figures/stroma_cells_traditional_top20_up.png)

#### Random Forest top genes

![Stroma Random Forest](figures/stroma_cells_random_forest_top20.png)

#### KNN top genes

![Stroma KNN](figures/stroma_cells_knn_single_gene_top20.png)

#### ML top CKD-up biomarkers

![Stroma ML up](figures/stroma_cells_ml_top20_up.png)

#### Consensus score comparison

![Stroma consensus comparison](figures/stroma_cells_consensus_score_comparison.png)

#### Consensus heatmap

![Stroma consensus heatmap](figures/stroma_cells_consensus_score_heatmap.png)

#### Consensus overlap counts

![Stroma consensus overlap](figures/stroma_cells_consensus_overlap_counts.png)

#### Model cross-validation balanced accuracy

![Stroma model CV](figures/stroma_cells_model_cv_balanced_accuracy.png)
