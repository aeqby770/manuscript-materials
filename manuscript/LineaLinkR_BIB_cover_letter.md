# Cover Letter — Briefings in Bioinformatics

Dear Editor,

We submit **"LineaLinkR: lineage-resolved temporal ligand–receptor mapping along branched single-cell trajectories"** for consideration as a **Problem solving protocol** in *Briefings in Bioinformatics*.

Trajectory inference (TI) and cell–cell communication (CCC) are commonly run as separate analyses. Static CCC scores do not locate an interaction along a lineage or relative to a fate bifurcation, while TI methods order cells without locating LR programs relative to lineage bifurcations. LineaLinkR combines these operations in one R workflow and reports each ligand–receptor record by temporal bin, lineage and bifurcation-relative timing pattern.

LineaLinkR evaluates four expression-only attractor endpoint candidates, selects one by directionality, topology and stability evidence, and uses a confidence-regularized ensemble below the prespecified 0.45 threshold. It then scores LR activity in equal-frequency bins along each soft lineage, detects bifurcation bins from entropy decline and assigns each LR trajectory to a five-class operational taxonomy used as an exploratory timing summary. Benchmarking uses a complete 19-dataset ordering core shared by LineaLinkR and 11 comparators. Five scale-extension datasets and 13 native-endpoint or application datasets are analysed under separate designs. Across four two-terminal simulations scored after the shared pre-branch segment, one-to-one branch macro F1 is **0.869** versus **0.564** for slingshot; lineage ARI is **0.800** versus **0.105**. A separate 600-replicate module-level LR-score simulation gives fork-evidence AUROC **0.975**, sensitivity **0.913** and specificity **1.000**. Pancreas pseudotime correlates with stochastic RNA velocity (**ρ = 0.936**). In ten matched colorectal-cancer patients, all six independently prespecified fibroblast–macrophage LR expression proxies increase in tumour (BH **q ≤ 0.037**), and the myofibroblast fraction within stromal cells increases in every patient. Blinded gastruloid and CellTag-multi evaluations show that endpoint selection remains important for temporal direction, while the CellTag data retain non-random clonal structure.

The 19-dataset matrix is the primary ordering comparison. Scale extensions and native-endpoint analyses are reported separately. Ordering and simulated-branch metrics assess trajectory performance; CCC overlap, pathway-prior agreement and lag analyses assess concordance and timing structure.

Source code, the benchmark data registry, analysis scripts and result tables will accompany the LineaLinkR v0.6.4 release (R 4.5.3; `set.seed(12345)`). The software is released under the MIT license, and a Zenodo archive DOI will be added before publication.

This manuscript is original, is not under consideration elsewhere, and will be approved by all authors before submission. The authors declare no competing interests. [Add related-manuscript disclosure if applicable — AUTHOR PLACEHOLDER.]

Sincerely,

Yanshuo Han  
School of Chemical Engineering, Ocean and Life Sciences, Dalian University of Technology  
Department of Vascular Surgery, Central Hospital of Dalian University of Technology  
yanshuohan@dlut.edu.cn  

Yanxiao Han  
School of Nursing, Liaoning University of Traditional Chinese Medicine  
yanxiao.han@aliyun.com
