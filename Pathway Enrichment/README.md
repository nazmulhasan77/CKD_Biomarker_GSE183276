# পথওয়ে এনরিচমেন্ট বিশ্লেষণ

এই ফোল্ডারে CKD candidate biomarker gene list ব্যবহার করে pathway enrichment analysis করা হয়েছে। আগের traditional differential expression, Random Forest, KNN এবং consensus biomarker result থেকে gene list নেওয়া হয়েছে। এরপর Enrichr ব্যবহার করে GO Biological Process, KEGG, Reactome, WikiPathways এবং MSigDB Hallmark database-এ enrichment পরীক্ষা করা হয়েছে।

## কেন এই বিশ্লেষণ করা হয়েছে

শুধু biomarker gene বের করলে বোঝা যায় কোন gene CKD-তে গুরুত্বপূর্ণ হতে পারে। কিন্তু pathway enrichment করলে বোঝা যায় সেই gene-গুলো কোন biological process, signaling pathway, inflammation, fibrosis, cell death, mitochondrial metabolism বা immune response-এর সাথে যুক্ত। তাই thesis interpretation-এর জন্য এই অংশটি খুব গুরুত্বপূর্ণ।

## কী কী input ব্যবহার করা হয়েছে

| Gene list | উৎস | Gene সংখ্যা |
|---|---:|---:|
| Whole Dataset Consensus ML RF | `Using ML/results/consensus_traditional_ml_rf_biomarkers.csv` | 99 |
| Whole Dataset CKD Up Traditional | `results/ckd_upregulated_candidate_biomarkers.csv` | 200 |
| Whole Dataset CKD Down Traditional | `results/ckd_downregulated_reference_loss_genes.csv` | 200 |
| Epithelial Cells Consensus | `Cell Type Biomrker/results/epithelial_cells_consensus_traditional_ml_rf.csv` | 118 |
| Immune Cells Consensus | `Cell Type Biomrker/results/immune_cells_consensus_traditional_ml_rf.csv` | 16 |
| Endothelial Cells Consensus | `Cell Type Biomrker/results/endothelial_cells_consensus_traditional_ml_rf.csv` | 58 |
| Stroma Cells Consensus | `Cell Type Biomrker/results/stroma_cells_consensus_traditional_ml_rf.csv` | 38 |

Whole CKD up এবং down list-এ top 200 gene নেওয়া হয়েছে, কারণ অনেক বেশি gene দিলে enrichment খুব broad হয়ে যায়। Consensus list-এ যত gene পাওয়া গেছে সব ব্যবহার করা হয়েছে।

## কোন database ব্যবহার করা হয়েছে

| Short name | Enrichr library |
|---|---|
| GO Biological Process | `GO_Biological_Process_2026` |
| KEGG | `KEGG_2026` |
| Reactome | `Reactome_Pathways_2024` |
| WikiPathways | `WikiPathways_2024_Human` |
| MSigDB Hallmark | `MSigDB_Hallmark_2020` |

Significant pathway ধরা হয়েছে `adjusted_p_value < 0.05` হলে। Overlap gene কমপক্ষে 2টি হতে হয়েছে।

## কীভাবে run করা হয়েছে

PowerShell থেকে নিচের command চালানো হয়েছে:

```powershell
powershell -ExecutionPolicy Bypass -File "Pathway Enrichment\run_pathway_enrichment.ps1"
```

মূল Python code:

```text
Pathway Enrichment/pathway_enrichment_pipeline.py
```

এই code আগের result CSV পড়ে gene list বানায়, Enrichr API-তে submit করে, enrichment result save করে এবং plot তৈরি করে।

## প্রধান ফলাফল

| Gene list | সবচেয়ে শক্তিশালী biological signal |
|---|---|
| CKD downregulated genes | Mitochondrial respiration, oxidative phosphorylation, electron transport chain |
| CKD upregulated genes | KRAS signaling, nephric duct/developmental repair signal, p53-related stress signal |
| Whole consensus biomarkers | WNT signaling এবং inflammatory response |
| Epithelial consensus biomarkers | Interferon gamma response, apoptosis, regulated necrosis |
| Immune consensus biomarkers | Toll-like receptor cascade, MAPK signaling, IL18 signaling, p53 pathway |
| Endothelial consensus biomarkers | TGF-beta signaling, transcriptional regulation, vascular/developmental pathway |
| Stroma consensus biomarkers | Developmental/transcription-related GO terms; KEGG/Reactome/Hallmark-এ strong significant result কম |

সবচেয়ে শক্তিশালী overall signal এসেছে CKD downregulated gene list থেকে। সেখানে mitochondrial oxidative phosphorylation এবং respiratory electron transport pathway খুব বেশি enriched হয়েছে। এর মানে CKD tissue/cell population-এ energy metabolism loss বা mitochondrial dysfunction একটি বড় biological theme হতে পারে।

## Significant pathway count

![Significant pathway count](figures/enrichment_significant_term_counts.png)

এই heatmap-এ প্রতিটি gene list এবং database অনুযায়ী significant enriched term সংখ্যা দেখানো হয়েছে। CKD downregulated gene list-এ সবচেয়ে বেশি significant pathway পাওয়া গেছে, বিশেষ করে Reactome, GO Biological Process, WikiPathways এবং KEGG-তে।

## Combined pathway heatmap

![Combined top enriched pathways heatmap](figures/combined_top_enriched_pathways_heatmap.png)

এই heatmap-এ top enriched pathway-গুলো বিভিন্ন gene list-এ কতটা strong তা দেখানো হয়েছে। রঙ যত গাঢ়, adjusted p-value তত ছোট, অর্থাৎ enrichment তত শক্তিশালী।

## Whole dataset consensus biomarker pathway

![Whole consensus top enriched terms](figures/whole_consensus_ml_rf_top_enriched_terms.png)

![Whole consensus library dotplot](figures/whole_consensus_ml_rf_library_dotplot.png)

Whole dataset consensus biomarker list-এ WNT signaling pathway এবং inflammatory response significant এসেছে। এগুলো CKD injury, repair এবং inflammation-এর সাথে biologically relevant।

## Whole dataset CKD upregulated pathway

![Whole CKD up top enriched terms](figures/whole_ckd_up_traditional_top_enriched_terms.png)

![Whole CKD up library dotplot](figures/whole_ckd_up_traditional_library_dotplot.png)

CKD upregulated gene list-এ KRAS signaling, nephric duct formation এবং stress/development-related pathway signal পাওয়া গেছে। এগুলো injury repair এবং maladaptive epithelial response-এর সাথে ব্যাখ্যা করা যায়।

## Whole dataset CKD downregulated pathway

![Whole CKD down top enriched terms](figures/whole_ckd_down_traditional_top_enriched_terms.png)

![Whole CKD down library dotplot](figures/whole_ckd_down_traditional_library_dotplot.png)

CKD downregulated gene list-এ সবচেয়ে strong signal হলো oxidative phosphorylation, mitochondrial ATP synthesis এবং respiratory electron transport। এই result CKD-তে mitochondrial energy metabolism কমে যাওয়ার ধারণাকে support করে।

## Epithelial cell pathway

![Epithelial top enriched terms](figures/epithelial_consensus_top_enriched_terms.png)

![Epithelial library dotplot](figures/epithelial_consensus_library_dotplot.png)

Epithelial consensus biomarker list-এ interferon gamma response, apoptosis এবং regulated necrosis pathway enriched হয়েছে। CKD-তে epithelial injury, inflammatory stress এবং cell death response ব্যাখ্যার জন্য এই result গুরুত্বপূর্ণ।

## Immune cell pathway

![Immune top enriched terms](figures/immune_consensus_top_enriched_terms.png)

![Immune library dotplot](figures/immune_consensus_library_dotplot.png)

Immune consensus biomarker list-এ Toll-like receptor cascade, MAPK signaling, IL18 signaling এবং p53 pathway enriched হয়েছে। এই result CKD immune activation এবং inflammatory signaling-এর সাথে সম্পর্কিত।

## Endothelial cell pathway

![Endothelial top enriched terms](figures/endothelial_consensus_top_enriched_terms.png)

![Endothelial library dotplot](figures/endothelial_consensus_library_dotplot.png)

Endothelial consensus biomarker list-এ TGF-beta signaling এবং transcription regulation pathway এসেছে। CKD-তে endothelial dysfunction, fibrosis-associated signaling এবং vascular remodeling ব্যাখ্যার জন্য এটি ব্যবহার করা যায়।

## Stroma cell pathway

![Stroma top enriched terms](figures/stroma_consensus_top_enriched_terms.png)

![Stroma library dotplot](figures/stroma_consensus_library_dotplot.png)

Stroma consensus biomarker list-এ GO Biological Process-এ developmental এবং transcription-related signal পাওয়া গেছে। তবে KEGG, Reactome, WikiPathways এবং Hallmark database-এ strong FDR-significant result কম, তাই stroma result interpret করার সময় সতর্ক থাকতে হবে।

## Output files

| File | কাজ |
|---|---|
| `results/input_gene_sets.csv` | কোন gene list-এ কোন gene ব্যবহার হয়েছে |
| `results/input_gene_set_summary.csv` | gene list summary |
| `results/all_enrichment_results.csv` | সব database এবং সব gene list-এর full enrichment result |
| `results/significant_enriched_terms_fdr_lt_0_05.csv` | শুধু significant pathway |
| `results/top10_terms_per_gene_set_and_library.csv` | প্রতি gene list ও database-এর top 10 term |
| `results/enrichment_count_summary.csv` | significant pathway count summary |
| `logs/resolved_enrichr_libraries.json` | কোন Enrichr library version ব্যবহার হয়েছে |

## Thesis-এ কীভাবে লিখতে পারো

এই বিশ্লেষণে candidate biomarker gene-গুলোর biological interpretation করার জন্য over-representation based pathway enrichment করা হয়েছে। Enrichr database ব্যবহার করে প্রতিটি biomarker gene set-এর enriched pathway বের করা হয়। Multiple testing correction-এর জন্য adjusted p-value ব্যবহার করা হয়েছে এবং `adjusted_p_value < 0.05` pathway-কে significant ধরা হয়েছে। Result অনুযায়ী CKD downregulated genes সবচেয়ে বেশি mitochondrial oxidative phosphorylation এবং respiratory electron transport pathway-তে enriched ছিল, যা CKD-তে mitochondrial dysfunction ও energy metabolism loss নির্দেশ করে। অন্যদিকে CKD upregulated এবং consensus biomarkers inflammatory response, WNT signaling, KRAS signaling, apoptosis, Toll-like receptor cascade এবং TGF-beta signaling pathway-এর সাথে যুক্ত, যা injury repair, immune activation এবং fibrosis-associated signaling-এর সাথে সম্পর্কিত।

## সীমাবদ্ধতা

Pathway enrichment gene list-এর উপর নির্ভর করে, তাই gene selection threshold বদলালে result কিছুটা বদলাতে পারে। Enrichr একটি over-representation analysis করে; এটি pathway activity সরাসরি মাপে না। Pathway activity confirm করতে module score, GSEA, ligand-receptor analysis বা external validation dataset ব্যবহার করা যেতে পারে।
