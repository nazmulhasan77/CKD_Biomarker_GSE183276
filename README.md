# দীর্ঘস্থায়ী কিডনি রোগের সম্ভাব্য বায়োমার্কার শনাক্তকরণ

## গবেষণার শিরোনাম

**একক-কোষ স্থানিক ট্রান্সক্রিপ্টোমিক উপাত্ত ব্যবহার করে দীর্ঘস্থায়ী কিডনি রোগের সম্ভাব্য বায়োমার্কার শনাক্তকরণ**

এই প্রকল্পে `GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Seurat_03282022.h5ad` উপাত্ত ব্যবহার করে দীর্ঘস্থায়ী কিডনি রোগ, অর্থাৎ CKD-সংশ্লিষ্ট সম্ভাব্য জিন বায়োমার্কার শনাক্ত করা হয়েছে।

এই উপাত্তটি মানব কিডনির সুস্থ ও আঘাতপ্রাপ্ত অবস্থার একক-কোষ অ্যাটলাস থেকে নেওয়া। এতে কোষভিত্তিক জিন প্রকাশ, রোগাবস্থা, রোগী পরিচয়, কোষ শ্রেণি, কোষ উপশ্রেণি, কোষীয় অবস্থা এবং UMAP অ্যানোটেশন আছে।

গুরুত্বপূর্ণ বিষয় হলো, এই নির্দিষ্ট `.h5ad` ফাইলে সরাসরি টিস্যুর স্থানিক x-y স্থানাঙ্ক নেই। তাই এখানে করা বিশ্লেষণ মূলত একক-কোষ ট্রান্সক্রিপ্টোমিক ও অ্যানোটেশনভিত্তিক। সরাসরি স্থানিক অবস্থান যাচাই করতে হলে সংশ্লিষ্ট spatial transcriptomics উপাত্ত প্রয়োজন হবে।

## ব্যবহৃত উপাত্ত

ইনপুট ফাইল:

```text
GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Seurat_03282022.h5ad
```

উপাত্তের আকার:

```text
109,741 কোষ x 37,080 জিন
```

রোগাবস্থাভিত্তিক কোষ সংখ্যা:

```text
CKD    52,314 কোষ
AKI    35,777 কোষ
Ref    21,650 কোষ
```

রোগী সংখ্যা:

```text
CKD    15 জন রোগী
AKI    12 জন রোগী
Ref    18 জন রেফারেন্স দাতা
```

কোষ শ্রেণিভিত্তিক সংখ্যা:

```text
epithelial cells      80,625
immune cells          17,176
endothelial cells      8,550
stroma cells           3,383
neural cells               7
```

প্রধান কোষ উপশ্রেণি:

```text
PT, TAL, IMM, EC, IC, PC, CNT, DCT, VSM/P, DTL, PEC, ATL, FIB, POD, PapE, NEU
```

## কোন expression matrix ব্যবহার করা হয়েছে

এই `.h5ad` ফাইলে `adata.X`-এ scaled expression আছে। সেখানে ঋণাত্মক মানও আছে। তাই বায়োমার্কার বিশ্লেষণের জন্য `adata.X` সরাসরি ব্যবহার করা হয়নি।

বিশ্লেষণের জন্য ব্যবহার করা হয়েছে:

```python
adata.raw.X
```

জিনের নাম নেওয়া হয়েছে:

```python
adata.raw.var["_index"]
```

কারণ `adata.raw.X` sparse log-normalized expression matrix হিসেবে আছে, যা differential expression এবং যন্ত্রশিক্ষা বিশ্লেষণের জন্য বেশি উপযুক্ত।

## কেন রোগী-স্তরের ছদ্ম-বাল্ক পদ্ধতি ব্যবহার করা হয়েছে

একক-কোষ উপাত্তে কোষের সংখ্যা অনেক বেশি হলেও প্রকৃত জৈবিক পুনরাবৃত্তি হলো রোগী। একই রোগীর অনেক কোষকে আলাদা আলাদা স্বাধীন নমুনা হিসেবে ধরলে pseudo-replication তৈরি হয়। এতে p-value খুব ছোট দেখাতে পারে এবং মডেলের accuracy অতিরিক্ত ভালো মনে হতে পারে।

এই সমস্যা এড়াতে রোগী-স্তরের ছদ্ম-বাল্ক পদ্ধতি ব্যবহার করা হয়েছে।

পদ্ধতি:

১. একই রোগীর সব কোষ একত্র করা হয়েছে।

২. প্রতিটি রোগীর জন্য প্রতিটি জিনের গড় log-normalized expression গণনা করা হয়েছে।

৩. এরপর CKD রোগী এবং Ref দাতার মধ্যে তুলনা করা হয়েছে।

সারকথা:

```text
কোষ-স্তরের নয়, রোগী-স্তরের বিশ্লেষণ করা হয়েছে যাতে ফলাফল পরিসংখ্যানগতভাবে বেশি গ্রহণযোগ্য হয়।
```

## এই প্রকল্পে কী কী করা হয়েছে

এই প্রকল্পে চার ধরনের প্রধান কাজ করা হয়েছে:

১. সম্পূর্ণ উপাত্ত ব্যবহার করে CKD বনাম Ref differential expression বিশ্লেষণ।

২. CKD বনাম AKI তুলনা করে CKD-নির্দিষ্টতার প্রাথমিক যাচাই।

৩. কোষ উপশ্রেণিভিত্তিক CKD বনাম Ref marker বিশ্লেষণ।

৪. Random Forest এবং KNN ব্যবহার করে যন্ত্রশিক্ষাভিত্তিক বায়োমার্কার ranking।

## প্রকল্পের ফোল্ডার বিন্যাস

```text
docs/
figures/
results/
scripts/
Using ML/
GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Seurat_03282022.h5ad
README.md
```

`scripts/` ফোল্ডারে differential expression pipeline আছে।

`results/` ফোল্ডারে differential expression ফলাফল আছে।

`figures/` ফোল্ডারে differential expression এবং UMAP-সংক্রান্ত চিত্র আছে।

`Using ML/` ফোল্ডারে Random Forest ও KNN বিশ্লেষণের কোড, ফলাফল এবং চিত্র আছে।

`docs/` ফোল্ডারে অতিরিক্ত ব্যাখ্যামূলক নথি আছে।

## বিশ্লেষণ ১: সম্পূর্ণ উপাত্ত থেকে CKD বনাম Ref বায়োমার্কার

ব্যবহৃত স্ক্রিপ্ট:

```text
scripts/ckd_biomarker_pipeline.py
```

প্রধান তুলনা:

```text
CKD বনাম Ref
```

এই বিশ্লেষণে যা করা হয়েছে:

১. সব কোষ থেকে রোগী-স্তরের ছদ্ম-বাল্ক expression profile তৈরি করা হয়েছে।

২. CKD রোগী এবং Ref দাতার মধ্যে gene-wise differential expression করা হয়েছে।

৩. Welch t-test ব্যবহার করা হয়েছে।

৪. Benjamini-Hochberg পদ্ধতিতে FDR সংশোধন করা হয়েছে।

৫. Cohen's d effect size গণনা করা হয়েছে।

৬. প্রতিটি জিন CKD ও Ref অবস্থায় কত শতাংশ কোষে প্রকাশিত হয়েছে তা গণনা করা হয়েছে।

৭. p-value, effect size, expression difference এবং expression prevalence মিলিয়ে ranking score তৈরি করা হয়েছে।

প্রধান আউটপুট:

```text
results/ckd_vs_ref_all_cells_de.csv
results/ckd_candidate_biomarkers_ranked.csv
results/ckd_upregulated_candidate_biomarkers.csv
results/ckd_downregulated_reference_loss_genes.csv
```

## সম্পূর্ণ উপাত্ত থেকে CKD-তে বেশি প্রকাশিত সম্ভাব্য বায়োমার্কার

আউটপুট ফাইল:

```text
results/ckd_upregulated_candidate_biomarkers.csv
```

শীর্ষ CKD-upregulated সম্ভাব্য বায়োমার্কার:

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

ব্যাখ্যা:

`diff_mean_log_expr` ধনাত্মক হলে জিনটি CKD-তে Ref-এর তুলনায় বেশি প্রকাশিত। তাই ওপরের জিনগুলো CKD বনাম সুস্থ রেফারেন্স অবস্থায় সম্ভাব্য upregulated বায়োমার্কার।

## CKD-তে কম প্রকাশিত জিন

আউটপুট ফাইল:

```text
results/ckd_downregulated_reference_loss_genes.csv
```

শীর্ষ CKD-downregulated জিন:

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

ব্যাখ্যা:

`diff_mean_log_expr` ঋণাত্মক হলে জিনটি CKD-তে Ref-এর তুলনায় কম প্রকাশিত। এগুলোকে সরাসরি CKD-upregulated biomarker বলা ঠিক নয়। এগুলোকে CKD-associated downregulated gene বা reference-loss marker বলা বেশি উপযুক্ত।

## CKD বনাম AKI নির্দিষ্টতা যাচাই

কোনো জিন CKD-তে বেশি প্রকাশিত হলেও সেটি AKI-তেও বেশি প্রকাশিত হতে পারে। সে ক্ষেত্রে জিনটি CKD-নির্দিষ্ট না হয়ে সাধারণ kidney injury marker হতে পারে।

এই কারণে CKD বনাম AKI তুলনাও করা হয়েছে।

আউটপুট:

```text
results/ckd_vs_aki_all_cells_de.csv
results/ckd_upregulated_and_aki_specific_biomarkers.csv
```

গুরুত্বপূর্ণ ফলাফল:

```text
এই বিশ্লেষণে এমন কোনো all-cell CKD-upregulated জিন পাওয়া যায়নি যা একই সঙ্গে CKD বনাম Ref FDR < 0.05 এবং CKD বনাম AKI FDR < 0.05 পূরণ করে।
```

তাই সিদ্ধান্ত:

```text
উল্লিখিত upregulated জিনগুলো CKD বনাম healthy reference candidate biomarker। এগুলোকে CKD-specific diagnostic biomarker বলতে হলে স্বাধীন validation প্রয়োজন।
```

## বিশ্লেষণ ২: কোষ উপশ্রেণিভিত্তিক বায়োমার্কার

সম্পূর্ণ উপাত্তভিত্তিক বিশ্লেষণের পাশাপাশি কোষ উপশ্রেণি অনুযায়ী আলাদা বিশ্লেষণ করা হয়েছে।

ব্যবহৃত অ্যানোটেশন:

```text
obs["subclass.l1"]
```

অন্তর্ভুক্ত কোষ উপশ্রেণি:

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

আউটপুট:

```text
results/subclass_l1_ckd_vs_ref_top_markers.csv
```

এই বিশ্লেষণে যা করা হয়েছে:

১. প্রতিটি কোষ উপশ্রেণির মধ্যে CKD বনাম Ref তুলনা করা হয়েছে।

২. প্রতিটি উপশ্রেণির জন্য রোগী-স্তরের ছদ্ম-বাল্ক profile তৈরি করা হয়েছে।

৩. প্রতিটি উপশ্রেণির top overall, CKD-up এবং CKD-down marker সংরক্ষণ করা হয়েছে।

গুরুত্বপূর্ণ কলাম:

```text
subclass.l1
marker_set
gene
rank_score
adj_p_value
diff_mean_log_expr
cohen_d
```

`marker_set` কলামের অর্থ:

```text
top_overall  = সবচেয়ে শক্তিশালী marker; up এবং down দুই ধরনের জিনই থাকতে পারে
up_in_CKD    = CKD-তে বেশি প্রকাশিত
down_in_CKD  = CKD-তে কম প্রকাশিত
```

উদাহরণ:

```text
PT:   ERRFI1, GADD45A, PDK4, USP2, TSC22D3 CKD-তে কম; CLU ও S100A11 CKD-up signal দেখায়
TAL:  ZFAND5, USP2, CEBPD, PDK4 CKD-তে কম; WFDC2 CKD-up signal দেখায়
EC:   MT2A, KLF9, MT1E CKD-তে কম; SOX18 ও MARCKS CKD-up signal দেখায়
IMM:  CXCR4, SRGN, TSC22D3 CKD-তে কম; ZFP36L1 CKD-up signal দেখায়
```

## বিশ্লেষণ ৩: Random Forest এবং KNN ব্যবহার করে যন্ত্রশিক্ষাভিত্তিক বায়োমার্কার

যন্ত্রশিক্ষা বিশ্লেষণের ফোল্ডার:

```text
Using ML/
```

ব্যবহৃত স্ক্রিপ্ট:

```text
Using ML/01_build_patient_pseudobulk.py
Using ML/02_random_forest_knn_biomarkers.py
```

এই ধাপে যা করা হয়েছে:

১. CKD ও Ref রোগীদের জন্য রোগী-স্তরের ছদ্ম-বাল্ক expression matrix তৈরি করা হয়েছে।

২. Random Forest classifier প্রশিক্ষণ করা হয়েছে।

৩. K-Nearest Neighbors, অর্থাৎ KNN classifier প্রশিক্ষণ করা হয়েছে।

৪. Random Forest থেকে gene feature importance বের করা হয়েছে।

৫. KNN-এর নিজস্ব coefficient নেই, তাই KNN-এর জন্য single-gene cross-validation importance গণনা করা হয়েছে।

৬. KNN permutation importance আলাদাভাবে সংরক্ষণ করা হয়েছে।

৭. Random Forest importance, KNN evidence এবং statistical score একত্র করে combined ML biomarker score তৈরি করা হয়েছে।

যন্ত্রশিক্ষায় ব্যবহৃত রোগী-স্তরের উপাত্ত:

```text
CKD বনাম Ref
15 জন CKD রোগী
18 জন Ref দাতা
মোট 33টি রোগী-স্তরের profile
29,010টি non-constant gene
```

মডেলের cross-validation ফলাফল:

```text
Random Forest CV balanced accuracy: 1.00
KNN CV balanced accuracy: প্রায় 0.94
```

সতর্কতা:

রোগীর সংখ্যা কম হওয়ায় যন্ত্রশিক্ষার ফলাফলকে সম্ভাব্য বায়োমার্কার শনাক্তকরণের সহায়ক প্রমাণ হিসেবে ব্যবহার করতে হবে। স্বাধীন validation ছাড়া clinical diagnostic biomarker দাবি করা যাবে না।

## যন্ত্রশিক্ষা আউটপুট ফাইল

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

## যন্ত্রশিক্ষা থেকে CKD-তে বেশি প্রকাশিত শীর্ষ সম্ভাব্য বায়োমার্কার

আউটপুট:

```text
Using ML/results/ml_upregulated_candidate_biomarkers.csv
```

শীর্ষ জিন:

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

Differential expression এবং যন্ত্রশিক্ষা দুই দিক থেকেই সমর্থন পাওয়া গুরুত্বপূর্ণ overlapping জিন:

```text
PROM1
MAP3K1
CDH6
TBL1XR1
VCAM1
SESN3
ITGB8
```

এই overlapping জিনগুলো থিসিসের discussion অংশে অগ্রাধিকারপ্রাপ্ত candidate biomarker হিসেবে ব্যবহার করা যেতে পারে।

## Combined ML Ranking

আউটপুট:

```text
Using ML/results/ml_combined_candidate_biomarkers.csv
```

শীর্ষ combined ML জিন:

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

গুরুত্বপূর্ণ ব্যাখ্যা:

এই combined তালিকায় CKD-তে বেশি এবং কম প্রকাশিত দুই ধরনের জিনই আছে। যদি উদ্দেশ্য হয় CKD-তে বেশি প্রকাশিত বায়োমার্কার খোঁজা, তাহলে `ml_upregulated_candidate_biomarkers.csv` ব্যবহার করা উচিত।

## তৈরি করা চিত্র

প্রধান চিত্র:

```text
figures/umap_condition_l1.png
figures/umap_class.png
figures/cell_class_composition_by_condition.png
figures/volcano_ckd_vs_ref_all_cells.png
figures/top20_ckd_upregulated_candidate_biomarkers.png
figures/top20_ckd_downregulated_reference_loss_genes.png
```

যন্ত্রশিক্ষা চিত্র:

```text
Using ML/figures/top20_random_forest_importance.png
Using ML/figures/top20_knn_single_gene_cv_importance.png
Using ML/figures/top20_knn_permutation_importance.png
Using ML/figures/top20_combined_ml_biomarkers.png
Using ML/figures/top20_ml_upregulated_biomarkers.png
```

## Differential expression pipeline চালানোর নিয়ম

Project root থেকে চালাতে হবে:

```powershell
python scripts/ckd_biomarker_pipeline.py
```

কোষ উপশ্রেণিভিত্তিক বিশ্লেষণ বাদ দিয়ে দ্রুত চালাতে চাইলে:

```powershell
python scripts/ckd_biomarker_pipeline.py --skip-subclass
```

কম memory ব্যবহার করে চালাতে চাইলে:

```powershell
python scripts/ckd_biomarker_pipeline.py --chunk-size 500
```

## যন্ত্রশিক্ষা pipeline চালানোর নিয়ম

Project root থেকে চালাতে হবে:

```powershell
python "Using ML\01_build_patient_pseudobulk.py"
python "Using ML\02_random_forest_knn_biomarkers.py"
```

অথবা একসাথে চালাতে:

```powershell
powershell -ExecutionPolicy Bypass -File "Using ML\run_using_ml_pipeline.ps1"
```

## গুরুত্বপূর্ণ কলামের অর্থ

`gene`

```text
জিনের নাম
```

`diff_mean_log_expr`

```text
CKD mean log-expression minus Ref mean log-expression
ধনাত্মক মান = CKD-তে বেশি প্রকাশ
ঋণাত্মক মান = CKD-তে কম প্রকাশ
```

`adj_p_value`

```text
FDR-adjusted p-value
সাধারণত 0.05-এর কম হলে তা পরিসংখ্যানগতভাবে গুরুত্বপূর্ণ ধরা হয়
```

`cohen_d`

```text
effect size
ধনাত্মক মান = CKD-তে বেশি প্রকাশ
ঋণাত্মক মান = CKD-তে কম প্রকাশ
মানের absolute value যত বেশি, পরিবর্তন তত শক্তিশালী
```

`rank_score`

```text
differential expression ভিত্তিক composite ranking score
```

`combined_ml_biomarker_score`

```text
Random Forest importance, KNN single-gene importance এবং statistical score মিলিয়ে তৈরি যন্ত্রশিক্ষাভিত্তিক ranking score
```

`random_forest_importance`

```text
Random Forest মডেলে জিনটির গুরুত্ব
```

`knn_single_gene_balanced_accuracy_mean`

```text
একটি জিন একা ব্যবহার করে KNN cross-validation-এ CKD বনাম Ref আলাদা করার গড় balanced accuracy
```

## চূড়ান্ত candidate gene বাছাইয়ের প্রস্তাবিত কৌশল

থিসিসের জন্য শক্তিশালী candidate biomarker তালিকা বানাতে চাইলে নিচের বিষয়গুলো বিবেচনা করা যেতে পারে:

১. জিনটি CKD বনাম Ref বিশ্লেষণে CKD-তে বেশি প্রকাশিত হওয়া উচিত।

২. `adj_p_value < 0.05` হওয়া ভালো।

৩. `cohen_d` ধনাত্মক এবং তুলনামূলকভাবে বড় হওয়া ভালো।

৪. Random Forest বা KNN ranking-এ সমর্থন থাকা ভালো।

৫. কোষ উপশ্রেণিভিত্তিক localization থাকলে আরও ভালো।

৬. পূর্ববর্তী সাহিত্য থেকে biological support থাকা প্রয়োজন।

অগ্রাধিকার দিয়ে আলোচনা করা যেতে পারে এমন candidate gene:

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

## থিসিসের Methods অংশে ব্যবহারযোগ্য লেখা

```text
GSE183276 একক-কোষ ট্রান্সক্রিপ্টোমিক উপাত্ত ব্যবহার করে দীর্ঘস্থায়ী কিডনি রোগ-সংশ্লিষ্ট সম্ভাব্য বায়োমার্কার শনাক্ত করা হয়েছে। একই রোগীর কোষগুলো স্বাধীন জৈবিক পুনরাবৃত্তি নয় বলে pseudo-replication এড়াতে কোষভিত্তিক expression রোগী-স্তরের ছদ্ম-বাল্ক profile-এ রূপান্তর করা হয়েছে। এরপর CKD এবং healthy reference রোগী profile-এর মধ্যে differential expression বিশ্লেষণ করা হয়েছে। Multiple testing correction, effect size estimation এবং expression prevalence ব্যবহার করে candidate gene ranking করা হয়েছে। কোষ উপশ্রেণিভিত্তিক বিশ্লেষণের মাধ্যমে disease-associated expression signal কোন কিডনি কোষ population-এ বেশি তা নির্ণয় করা হয়েছে। অতিরিক্তভাবে Random Forest এবং K-Nearest Neighbors মডেল ব্যবহার করে রোগী-স্তরের expression profile থেকে যন্ত্রশিক্ষাভিত্তিক feature prioritization করা হয়েছে।
```

## থিসিসের Results অংশে ব্যবহারযোগ্য লেখা

```text
রোগী-স্তরের all-cell বিশ্লেষণে WFDC2, ITGB8, MAP3K1, TBL1XR1, LRRK2, TNFRSF11B, PROM1, MMP7 এবং VCAM1 CKD-তে বেশি প্রকাশিত সম্ভাব্য বায়োমার্কার হিসেবে শনাক্ত হয়েছে। Random Forest এবং KNN-ভিত্তিক যন্ত্রশিক্ষা বিশ্লেষণে PROM1, MAP3K1, CDH6, TBL1XR1, VCAM1, SESN3 এবং ITGB8-এর মতো কয়েকটি জিনের জন্য অতিরিক্ত সমর্থন পাওয়া গেছে। কোষ উপশ্রেণিভিত্তিক বিশ্লেষণে PT, TAL, EC, IMM এবং VSM/P সহ একাধিক কিডনি কোষ population-এ CKD-সংশ্লিষ্ট transcriptional পরিবর্তন দেখা গেছে।
```

## সীমাবদ্ধতা

গুরুত্বপূর্ণ সীমাবদ্ধতা:

```text
১. উচ্চমাত্রিক যন্ত্রশিক্ষা বিশ্লেষণের তুলনায় রোগীর সংখ্যা কম।
২. CKD group-এর মধ্যে DKD এবং hypertensive CKD দুটোই আছে।
৩. এই h5ad ফাইলটি scCv3 single-cell data; সরাসরি spatial coordinate matrix নয়।
৪. all-cell level-এ কঠোর CKD বনাম AKI-specific upregulated biomarker পাওয়া যায়নি।
৫. চূড়ান্ত diagnostic biomarker দাবি করার আগে স্বাধীন validation প্রয়োজন।
```

Validation-এর জন্য প্রস্তাবনা:

```text
স্বাধীন CKD cohort validation
bulk RNA-seq validation
external single-cell validation
protein-level validation
স্বাধীন নমুনায় ROC/AUC analysis
সাহিত্যভিত্তিক biological interpretation
```

## সংক্ষিপ্ত সারাংশ

এই প্রকল্পে যা করা হয়েছে:

```text
১. সম্পূর্ণ উপাত্ত থেকে CKD বনাম Ref candidate biomarker শনাক্ত করা হয়েছে।
২. কোষ উপশ্রেণিভিত্তিক CKD marker শনাক্ত করা হয়েছে।
৩. CKD বনাম AKI নির্দিষ্টতা যাচাই করা হয়েছে।
৪. Random Forest এবং KNN ব্যবহার করে যন্ত্রশিক্ষাভিত্তিক biomarker ranking করা হয়েছে।
৫. সব ফলাফল CSV ও figure আকারে সংরক্ষণ করা হয়েছে।
৬. pseudo-replication এড়াতে রোগী-স্তরের ছদ্ম-বাল্ক পদ্ধতি ব্যবহার করা হয়েছে।
```

সবচেয়ে গুরুত্বপূর্ণ ফলাফল ফাইল:

```text
results/ckd_upregulated_candidate_biomarkers.csv
results/subclass_l1_ckd_vs_ref_top_markers.csv
Using ML/results/ml_upregulated_candidate_biomarkers.csv
Using ML/results/ml_combined_candidate_biomarkers.csv
```

