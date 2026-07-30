# LineaLinkR: lineage-resolved temporal ligand–receptor mapping along branched single-cell trajectories

**Article type:** Problem solving protocol  
**Short title:** LineaLinkR temporal ligand–receptor inference  
**Authors:** Hongjie Jiang¹, Xinyu Zhang²†, Hong Xin Wang³†, Yanxiao Han⁴∗, Yanshuo Han⁵˒⁶∗  
**Affiliations:** ^1^ School of Life and Health Sciences, Hainan University, Haikou, Hainan, China; ^2^ Yingkou Central Hospital, Yingkou, Liaoning, China; ^3^ Liaoning Cancer Hospital & Institute, Shenyang, Liaoning, China; ^4^ School of Nursing, Liaoning University of Traditional Chinese Medicine, Shenyang, Liaoning, China; ^5^ School of Chemical Engineering, Ocean and Life Sciences, Dalian University of Technology, Panjin, Liaoning, China; ^6^ Department of Vascular Surgery, Central Hospital of Dalian University of Technology, Dalian, China  
**Corresponding authors (in order):** Yanshuo Han, School of Chemical Engineering, Ocean and Life Sciences, Dalian University of Technology, Panjin, Liaoning, China; Department of Vascular Surgery, Central Hospital of Dalian University of Technology, Dalian, China; Tel: +86 (427) 2631889; Fax: +86 (427) 2631889; Email: yanshuohan@dlut.edu.cn; ORCID: https://orcid.org/0000-0002-4897-2998. Yanxiao Han, School of Nursing, Liaoning University of Traditional Chinese Medicine, Shenyang, Liaoning, China; Email: yanxiao.han@aliyun.com  
**Equal contribution:** † Xinyu Zhang and Hong Xin Wang contributed equally to this work.  
**ORCID iDs:** Yanshuo Han: https://orcid.org/0000-0002-4897-2998  
**Word count (main text, excluding references, tables, and figure legends):** [BUILDER-CALCULATED]  
**Figure count:** [BUILDER-CALCULATED]  
**Table count:** [BUILDER-CALCULATED]  
**Keywords:** single-cell RNA-seq; trajectory inference; ligand–receptor communication; pseudotime; bifurcation; cell–cell signaling

---

## Abstract

Trajectory inference (TI) orders cells but does not place signaling on temporal or branch coordinates, whereas static cell–cell communication (CCC) lacks those coordinates. LineaLinkR combines expression-only rooting, branch-aware ordering, lineage×bin ligand–receptor scoring and bifurcation-relative timing patterns to turn a snapshot into a temporal CCC map. Across the 19-dataset complete ordering core, mean Spearman was 0.615 overall and 0.690 across 16 simulated ordering datasets. In four two-terminal simulations scored after the shared pre-branch segment, one-to-one branch macro F1 was 0.869 versus 0.564 for slingshot and lineage ARI was 0.800 versus 0.105. A separate 600-replicate module-level LR-score simulation yielded fork-evidence AUROC 0.975, sensitivity 0.913 and specificity 1.000. Pancreas blind pseudotime correlated with stochastic RNA velocity (ρ = 0.936), whereas GSE274389 and GSE216606 had the opposite orientation to recorded time. Operational role labels showed low agreement with prespecified FGF, Wnt and Notch pathway priors and are retained as exploratory timing summaries. In ten matched colorectal-cancer patients, all six prespecified expression axes increased from normal to tumour and the myofibroblast fraction increased in every patient. LineaLinkR reconstructs lineage-resolved temporal communication programs from snapshot scRNA-seq and prioritizes them for downstream mechanistic testing.

---

## Introduction

Differentiation couples gene-expression change to paracrine and juxtacrine signalling. Single-cell RNA sequencing (scRNA-seq) measures both processes [1, 2], but trajectory and communication analyses are usually separated.

Trajectory-inference (TI) methods estimate pseudotime and sometimes lineage membership from snapshot data [3–16]. Their inputs and endpoints differ. BayesTraj uses topology priors [17]; scRDEN couples reconstruction to gene-rank networks [18]; InterVelo, TIVelo and NeuroVelo use velocity [19–21]; Chronocell requires spliced and unspliced counts plus a trajectory structure [22]; PHLOWER uses multimodal or spatial data [23]; STORIES requires spatially resolved timepoints [24]; and scLANE tests expression dynamics on an existing trajectory [25]. The count-only benchmark includes methods whose native inputs are available in these datasets.

Cell–cell communication (CCC) methods score ligand–receptor (LR) pairs between clusters [26–34], usually without within-cluster time or branch context. TICCI integrates interaction scoring with trajectory graphs [35]; TimeTalk uses repeated timepoints and Granger tests [36, 37]; and TraSig scores profiles along an imported pseudotime [38]. These tools do not return the same bifurcation-relative operational labels within one snapshot workflow.

For example, static CCC detects Mdk–Ncl in mouse pancreas but cannot place its peak before or after endocrine commitment, or on a particular lineage. LineaLinkR retains lineage and temporal-bin coordinates for this distinction.

Its five-class taxonomy is operational: instructive denotes sender-leading pre-fork activity; reinforcing and feedback denote sender- and receiver-leading post-fork activity; lateral denotes pre-fork opposing-lineage specificity; and bystander is residual. An auxiliary transition score records activity at the fork. These labels describe timing patterns, not causal mechanisms.

LineaLinkR combines expression-only attractor endpoint selection, lineage-weighted temporal LR activity, uncertainty-derived fork bins and role scores. Benchmarking uses a complete 19-dataset ordering core shared by LineaLinkR and 11 comparators. Five scale-extension datasets and 13 native-endpoint or application datasets are analysed under separate designs.

---

## Methods

### Benchmark design

The 37-dataset registry contains 24 ordering datasets: a 19-dataset complete core (three real and 16 simulated datasets) and five scale extensions. The remaining 13 datasets support native-endpoint comparisons and biological applications. Pancreas [39] and dentategyrus [40] support developmental CCC analyses; five GEO cohorts support disease applications. Dataset roles are recorded in Supplementary Table S1.

Eleven TI comparators were evaluated through method-specific adapters in the 19-dataset ordering core. The ordering rank comprises the 12 methods listed in Supplementary Table S2. Methods that require velocity, repeated-time or spatial inputs were outside this count-only comparison.

Adapters preserve native inputs, roots and direction policies. Root-free comparator outputs may be orientation-aligned only under the prespecified rule; LineaLinkR retains its blind direction. The primary rank comparison contains 12 methods evaluated on the same 19 complete-core datasets (Supplementary Tables S2-S4). The five scale extensions are reported separately from this ranking.

Scored outputs use matched cell identifiers without imputation. The Supplementary Methods specify returned-cell, unmatched-group and orientation rules.

### LineaLinkR workflow

**Input and graph.** The runner normalizes each cell to the dataset median library size, applies `log1p`, selects the top 2000 dispersion-ranked genes and computes 30 principal components with a scanpy-compatible preprocessing replicate [41]. A 50-neighbour PCA graph supplies trajectory geometry. State-field validation rejects sample, patient, diagnosis and identifier columns; adaptive Leiden or k-means supplies missing groups. Ground-truth time and branch labels are withheld (Supplementary Methods). Default analysis parameters are listed in Supplementary Table S5.

**Expression-only root selection.** The workflow identifies terminal-attractor basins as peripheral, low-potency states on the expression-derived graph. A double-precision Rust kernel estimates absorption probabilities from each cell to these basins. Five min–max-scaled cell features are combined with fixed weights: detected-gene potency (0.25), fate entropy (0.20), distance to the nearest attractor (0.30), local source flow (0.15) and one minus maximum fate commitment (0.10). Source flow is the inverse-distance-weighted potency decrease from a cell to its neighbours. Neither reference time nor branch labels enter these calculations.

The reported `detected_genes_ordering` protocol evaluates four expression-only endpoint candidates: high- or low-detected-gene potency crossed with a state-aggregate or cellwise-component anchor. Each candidate receives fuzzy-graph pseudotime, oriented only by high detected-gene potency, and is scored by blind embedding alignment and ordering structure. The prespecified selector retains the default high-potency endpoint unless an alternative exceeds the alignment, ordering and representation gain thresholds (0.20, 0.10 and 0.00, respectively). If the best blind ordering score is below 0.45, endpoint selection uses a confidence-regularized PC/k perturbation ensemble. An independently justified user root changes the graph origin and is reported separately.

State anchors use the highest state-median attractor score and a directed-degree representative; cellwise anchors use the highest-scoring component candidate. The final start cell is the minimum of the selected oriented pseudotime, with high-potency and cell-index tie-breakers. Below the prespecified confidence threshold, the PC and neighbour perturbation grid supplies the endpoint ensemble. User-root analyses reuse the same connected weighted graph and recompute distances and downstream quantities, isolating the effect of the starting state without changing expression preprocessing.

**Pseudotime and lineages.** Locally scaled neighbour distance is converted to fuzzy membership and edge cost,

$$FUZZY_EDGE$$

where ρᵢ and σᵢ are local smooth-kNN parameters. The graph is symmetrized; disconnected components are joined by recorded geometric bridges, while cluster-MST edge augmentation is disabled in the reported analysis. Dijkstra distance from the root is scaled to [0,1] without rank normalization. Root-to-leaf cluster-tree paths define lineages. For cell *i* and path *l*,

$$LINEAGE_WEIGHT$$

with σ = 1 and row normalization. Path-specific cells receive unit unnormalized weight on that path; trunk cells retain Gaussian weights across paths. Rust log-sum-exp normalization prevents high-dimensional underflow.

Cluster-centroid paths define candidates, whereas cell-level distances remain on the weighted graph; added bridges are recorded. Pseudotime is a graph distance, not a fitted clock, and lineage weights are path-proximity scores rather than calibrated fate probabilities. Cell-level lineage assignments and uncertainty are shown in Supplementary Figures S1 and S2.

**Temporal ligand–receptor activity.** Counts are rescaled to 10,000 per cell and log-transformed. Each lineage has 12 weighted pseudotime-quantile bins. Weighted expression and detection are calculated for each ligand, receptor, source and target state; multi-subunit complexes take the minimum subunit value and require detection ≥0.05. Activity is

$$LR_ACTIVITY$$

where *L* and *R* are complex means and *S* is the dataset-wide mean positive gene-by-state expression, not a fixed constant. This is a scoring transform, not a kinetic model. Effective size is the smaller source/target Kish value; sensitivity sweeps are in Supplementary Tables S6 and S7 and Supplementary Figures S3 and S4.

Source and target states are evaluated separately. The grid retains zeros, non-evaluable combinations, detection fractions and effective sizes. Matched static scoring omits time. “Temporal-only” means absent from this static output, not delayed signalling.

The binning contract is fixed before LR ranking. Quantile boundaries are estimated from pseudotime with lineage weights, and the same boundaries are used for ligands, receptors and sample-level recomputation. A bin is structurally valid when both source and target Kish effective sizes are at least 5; otherwise its score is `NA`. Within a structurally valid bin, insufficient ligand or receptor detection sets the score to zero and `expression_supported` to false. Pair summaries are formed from the exported lineage-by-bin grid and retain these two states separately. When biological samples are available, pooled bins define the estimand and resampling occurs over samples; cells are not treated as independent replicates. In datasets without sample identifiers, fork contrasts and temporal profiles are descriptive. Static support from CellChat, LIANA or CellPhoneDB is an external concordance field and does not alter LineaLinkR scores or the temporal-only definition.

**Fork and role evidence.** Normalized Shannon entropy [42] and its lineage-bin mean are

$$ENTROPY$$

where *L* is the lineage count. For the three-bin-smoothed curve H̃ₗ(*b*), the fork bin is

$$BIFURCATION$$

Short or near-flat curves are uninformative. For pair *p* on lineage *l*, the sibling contrast is

$$FORK_CONTRAST$$

where pre/post windows use two bins per side; focal and sibling changes use 95th-percentile-capped effective-size weights, and the sibling contribution is normalized by the sum of its retained weights. With biological samples, activity is recomputed in fixed pooled bins and samples are bootstrapped; otherwise the contrast is pooled and descriptive. Bootstrap uncertainty for leading pancreas programs is shown in Supplementary Figure S5.

For a structurally valid row, let *s*ₗ be focal-lineage specificity, *s*₋ₗ the maximum opposing-lineage specificity, Δτ the target-minus-source mean pseudotime and τ₀ the lineage timing scale; the indicators mark bins before, at or after the fork. The exported role components are

$$ROLE_SCORES$$

The five-class calibration assigns the largest of instructive, reinforcing, feedback, lateral and bystander; the transition component is reported separately and is not included in that five-class label. These are dimensionless operational scores, and the pancreas pathway-prior comparison is exploratory.

Fork eligibility is lineage-specific; insufficient bins, no sustained entropy decline or no sibling yields no call. Role assignment follows expression, size and boundary gates. Pathway-prior agreement is evaluated at the trajectory level.

### Metrics and uncertainty

Predicted and reference pseudotime are cell-ID aligned and scaled to [0,1]. Metrics are Spearman ρ, Pearson *r*, Kendall τ, RMSE and MAE. Branch structure is evaluated in four two-terminal simulations using cells with reference pseudotime ≥0.4, after the shared pre-branch segment, and summarized by ARI, NMI, macro Jaccard and one-to-one macro F1. The other 12 simulated stress datasets contribute ordering endpoints only. Datasets are replicates: mean differences use dataset bootstraps, branch comparisons add exact sign-flip tests, and paired analyses use common datasets [43]. CCC overlap, recurrence, prior agreement and lag are concordance endpoints, not causal accuracy; detailed coverage and unmatched-group rules are in the Supplement.

A separate module-level simulation starts from synthetic lineage×bin LR-score arrays with known positive or null fork evidence. Six scenarios cover a shared trend, a focal branch, a donor-confounded null, a low-expression spike, 30% dropout and a nested fork. One hundred replicates per scenario are evaluated by fork-evidence AUROC, sensitivity, specificity and scenario-level detection rate. This simulation tests the fork-evidence layer from LR scores rather than end-to-end LR recovery from a simulated expression matrix.

The primary ordering endpoint is dataset-level mean Spearman, reported overall and by real/simulated stratum. Common-core ranks are computed within dataset; Friedman tests assess the overall method effect, and pairwise bootstrap intervals describe LineaLinkR-minus-comparator differences. One-to-one Hungarian branch F1 assigns zero to unmatched groups; many-to-one mapping is reported as a secondary sensitivity analysis.

### Ablation and validation datasets

Ablations remove LR scoring, temporal bins or fork detection; graph variants rerun the trajectory. Sensitivity analyses vary neighbours, bins and Hill scale. GSE274389 supplies an external time-course analysis with sampling time and Wnt-recorder labels withheld from inference [44].

The external holdout is the GSE216606 CellTag-multi mouse LSK experiment [45]. Nine GEO RNA matrices were matched to the clone table (37,380 traced cells) and assigned neutral IDs. Day, clone, fate, author state and technical batch were withheld; inference used unsupervised clusters. Prespecified endpoints are signed day concordance, state/fate agreement and within-clone lineage concentration. The four-candidate selector retained `detected_genes_high_state`, and all reported endpoints use this single returned root. Direction-agnostic AUC measures separability only.

The annotation-free input was fixed before inference and withheld fields were joined only after fitting. Clone concentration compares within-clone lineage-pair agreement with the global lineage-pair probability and reports normalized excess with a clone-bootstrap interval.

Validation units follow the data-generating design. Dataset is the replicate for the trajectory benchmark, patient is the replicate for matched CRC contrasts, biological sample is the replicate for muscle pseudobulk tests, and clone is the resampling unit for the CellTag concentration interval. Cell-level associations are reported as descriptive separability or ordering summaries. Prespecified endpoints were written before withheld metadata were joined: signed day concordance tests orientation, direction-agnostic AUC tests separation, ARI and NMI compare partitions, and clone concentration tests whether cells sharing a barcode occupy the same inferred lineage more often than expected from global lineage frequencies. Temporal direction and clonal structure are reported as separate endpoints.

Muscle analysis combines GSE143437 day 0/2/5/7 [46] with two-replicate GSE291091 FAP depletion [47]. Blind direction is compared with recorded day before day anchors the CCC analysis. CRC tests six prespecified fibroblast–macrophage axes in ten matched GSE132465 patients [48, 49]. GSE210152 is analysed as a pooled Leiden dataset with exact-barcode plaque-atlas mapping and patient-rescoring sensitivity analyses. Preprocessing, analysis units and endpoints for these applications are specified in the Supplementary Methods.

### Reproducibility

LineaLinkR 0.6.4 is an R package with a Rust core; scanpy performs preprocessing. Analyses use R 4.5.3, Python 3 and seed 12345. Configuration, input identifiers and root choice are recorded for each analysis. Runtime measurements use fresh processes and exclude optional taxonomy calibration.

The runtime boundary includes matrix loading, normalization, variable-gene selection, PCA, kNN construction, graph connection, root inference, pseudotime, lineage weights, static LR scoring, temporal aggregation, fork evidence and serialization. Optional UMAP rendering and taxonomy calibration are excluded. Wallclock and process-tree peak memory are measured in fresh processes; stage timers separate Python preprocessing, the Rust pseudotime kernel and R-level LR aggregation.

---

## Results

### Simulated branch recovery and pseudotime ordering

Across the 19 complete-core datasets, mean Spearman is 0.743 for TSCAN, 0.700 for MultistageOT, 0.615 for LineaLinkR, 0.524 for scTour, 0.504 for scFates and 0.504 for monocle3; all other means are 0.456 or lower. Across the 16 simulated ordering datasets, LineaLinkR mean Spearman is 0.690. Terminal-branch recovery is evaluated separately on four two-terminal simulations: one-to-one macro F1 is 0.869 for LineaLinkR and 0.564 for slingshot, while lineage ARI is 0.800 and 0.105 (Fig. 2–3; Supplementary Tables S8 and S9). Method-level correlations are shown in Supplementary Figure S6.

Within the four terminal-cell comparisons, mean Jaccard is 0.888 for LineaLinkR and 0.535 for slingshot. The paired one-to-one F1 difference is 0.304 (95% dataset-bootstrap CI −0.010–0.595; exact sign-flip *P* = 0.250), with positive differences in three of four datasets. The lineage-ARI difference is 0.695 (0.294–1.000; exact sign-flip *P* = 0.125), positive in all four datasets. Many-to-one sensitivity scores are 0.928 and 0.686 (Supplementary Table S9). The remaining simulated batch and dropout datasets are analysed as pseudotime-ordering stress tests rather than branch-label benchmarks.

A Friedman test across the 19-dataset, 12-method core detects a method effect (χ² = 74.13, 11 df, *P* = 1.99 × 10⁻¹¹; Supplementary Figure S7). Direction policy and returned-cell coverage are reported alongside the common-core comparison.

### Module-level simulation evaluates fork evidence

Across 600 module-level LR-score replicates, integrated fork-evidence AUROC was 0.975, sensitivity was 0.913 and specificity was 1.000 (Fig. 3; Supplementary Table S10). Detection was 100% for focal-branch and nested-fork signals, 74% under 30% dropout and 0% for shared-trend, donor-confounded-null and low-expression controls. This simulation isolates the fork-evidence calculation; it does not measure end-to-end LR-pair recovery from expression. The separate module-output dependence check confirms which result layers remain after removing LR scoring, temporal binning or fork detection (Supplementary Table S11 and Fig. S8).

### Exploratory role-label agreement with pathway priors

LineaLinkR assigns each LR trajectory a dominant operational role. Across five pancreas lineages, 1,247 expression-supported pathway trajectories entered the prespecified timing-class comparison: FGF and Wnt instructive, and Notch feedback. Expected-role agreement was low: 0.043 for FGF (19/441), 0.121 for Wnt (31/257) and 0.013 for Notch (7/549; Fig. 4; Supplementary Tables S12 and S13). The labels are therefore used as exploratory timing summaries rather than calibrated biological roles.

### Temporal LR programs overlap static CCC partially and add temporal-only pairs

Pancreas yielded 345 pair-level temporal summaries and 31 top-decile broad-lineage candidates (Fig. 5). CellChat supported 27 candidates, LIANA supported 14 and both supported 14; four met the temporal-only definition relative to the reported external call sets (Supplementary Table S14). Mdk–Ncl ranked first with a mean temporal score of 0.867. Other leading candidates included Mpzl1–Mpzl1, Mdk–Sdc4, Ptn–Ncl, F11r–F11r and Dlk1–Notch2. The ranking summarizes lineage-resolved temporal expression patterns and cross-method support.

Dentategyrus (2,930 cells, 14 types and five soft lineages; Fig. 6) yielded 292 pair-level summaries and 26 candidates. CellChat supported 24, LIANA supported 20 and both supported 20; two met the temporal-only definition. Psap–Gpr37l1 ranked first and was supported by both external methods. Prosaposin signalling through GPR37L1/GPR37 has been demonstrated in native astrocytes, supporting the biological identity of this glial signalling axis without establishing its temporal direction here [56]. Other leading programs included Nrxn1–Nlgn3, Efnb3–Epha4, Nrxn1–Nlgn2 and Apln–Aplnr. Developmental-case output counts are summarized in Supplementary Table S15.

### Expression timing of selected temporal-only pairs

For six selected temporal-only pairs we compared context-matched log-normalized ligand and receptor profiles without reusing the MM-LR score (Supplementary Table S16 and Fig. S9). Later receptor peaks occurred for Iapp–Calcrl (+10 bins), Pyy–Npy1r (+10), Pdgfa–Pdgfrb (+3) and Mif–Ackr3 (+1), whereas Ccl7–Ccr1 and Mdk–Lrp1 peaked in the same bin. Cross-correlation supported receptor-later shifts for Iapp–Calcrl (+5), Pdgfa–Pdgfrb (+4) and Pyy–Npy1r (+3). Together, peak and cross-correlation summaries distinguish synchronous and lagged receptor profiles. A complementary all-cell lag sensitivity analysis is reported in Supplementary Table S17.

### Patient-level colorectal-cancer data reproduce a prespecified fibroblast–macrophage program

All six prespecified broad-cell-type LR expression proxies increased from matched normal to tumour samples in GSE132465 (ten patients; BH *q* ≤ 0.037; Fig. 8 and Supplementary Table S18). The mean paired symmetric change ranged from 0.275 for RARRES2–CMKLR1 to 0.912 for SPP1–CD44. The myofibroblast fraction within stromal cells increased in all ten patients (median 0.002 to 0.379; mean paired difference 0.404, 95% bootstrap CI 0.287–0.534; Wilcoxon *P* = 0.00195). The SPP1-positive fraction within myeloid cells was not significantly increased (median 0.436 to 0.496; *P* = 0.275), but across all 23 tumour patients it correlated with the myofibroblast fraction (Spearman ρ = 0.491, *P* = 0.0174; patient as the unit). Fine-subtype expression proxies were evaluable in 17/23 tumours (Supplementary Figure S10). The matched-patient analysis reproduced the published fibroblast–SPP1-positive-macrophage expression program.

### Blind-root and user-root analyses

Blind-root mode remains the primary benchmark mode, whereas user-root mode is intended for applications with an independently justified start state. The user-root analysis (Supplementary Table S19) shows how root choice changes ordering and is reported separately from the blind benchmark.

### Orthogonal dynamic axes withheld from inference

Two independent axes were withheld from LineaLinkR inference (Supplementary Table S20 / Supplementary Figure S11). In pancreas, stochastic scVelo pseudotime from spliced/unspliced layers correlated with blind-root LineaLinkR pseudotime (Spearman 0.936; Kendall 0.779) [50, 51]. In gastrulation erythroid, blind-root and user-root pseudotime correlated with withheld embryo stage at 0.238 and 0.849, respectively [52].

Three kNN-inclusive Windows runs were performed on paul15, dentategyrus and pancreas. Mean total wallclock is 10.64, 15.36 and 18.80 s, with mean process-tree peak RSS 291, 335 and 519 MB. The LineaLinkR analysis stage, including static LR and temporal/fork scoring, averages 2.27, 6.77 and 6.66 s; Rust pseudotime alone averages 0.21–0.33 s (Supplementary Table S21 / Supplementary Figure S12). Temporal and fork aggregation is the largest measured LineaLinkR stage on these inputs.

### External gastruloid direction and Wnt-history association

GSE274389 retained 7942 mouse gastruloid cells after QC; sampling time and Wnt-recorder state were withheld. The selected root reverses the 96-h to 120-h transition (Spearman −0.788; day-120 AUC 0.027), and no alternative root met the retention criteria. The selected trajectory therefore has the opposite orientation to recorded time (Fig. 7; Supplementary Table S22).

Among 120-h cells, recorder-positive AUC is 0.453 for pseudotime, 0.604 for the Wnt-target control and 0.507 for the LineaLinkR Wnt score. Cluster-adjusted cell-level odds ratios are 0.633 per pseudotime IQR, 1.467 per Wnt-target IQR and 1.023 per LineaLinkR Wnt-score IQR; cells are not biological replicates. Lineage-bin concordance is 0.419 for the Wnt-target control and −0.010 for the LineaLinkR score. The Wnt-history association is weak and the inferred temporal direction is reversed.

### CellTag-multi holdout: direction and clonal structure

GSE216606 was converted to an annotation-free H5AD before analysis; day, author state, technical batch, clone and fate were withheld. On 37 380 traced cells, the selected high-state endpoint reverses sampling day (Spearman −0.383; signed day-5 AUC 0.050); direction-agnostic AUC is 0.950 (Fig. 7). Published-state agreement is ARI 0.257 and NMI 0.318, and day-5 fate agreement is ARI 0.232 and NMI 0.270. Among 5066 clones, median dominant-lineage fraction is 0.800; within-clone pair agreement is 0.532 versus 0.255 globally, giving normalized excess 0.371 (95% CI 0.354–0.390). The result separates an inverse temporal orientation from non-random clonal concentration (Supplementary Table S23).

### Muscle-regeneration timing and FAP-depletion response

GSE143437 blind pseudotime does not follow injury day (Spearman −0.060), so the application uses recorded day. The sample-level C3–C3AR1 means at days 0/2/5/7 are 0.844/0.661/0.709/0.625. The 2–3 samples per day support a phase-dependent expression pattern; a directional handoff was not tested (Fig. 8; Supplementary Table S24).

In GSE291091 FAP depletion (two replicates per group), pseudobulk analysis finds 317 FDR-significant genes in mesenchymal/FAP-like cells and 1360 in broad phagocytes, including reduced interferon, TNF and IL6 programs. C3 is not FDR-significant and phagocyte estimates are composition-sensitive (Supplementary Figure S13). The perturbation identifies a system-level FAP-dependent response; a C3-specific perturbation was not performed.

### Atherosclerosis scRNA-seq vignette (GSE210152)

LineaLinkR analysis of 34 456 GSE210152 cells yields seven Leiden clusters and pooled temporal LR rankings (Supplementary Figures S14 and S15). Subset and patient-level analyses assess the stability of these rankings; static CellPhoneDB means and LineaLinkR temporal scores remain method-specific endpoints.

The deposited H5AD has patient identifiers but no cell type. Exact barcode joining to the processed 2025 plaque atlas [53], without expression-based imputation, maps 25 193 cells (73.1%). Basic atlas QC is passed by 99.95% of mapped cells versus 37.36% of 9263 non-members. Leiden 0 is B-cell enriched, Leiden 4 T-cell enriched and Leiden 5 mast-cell enriched; Leiden 1 is mixed macrophage enriched. Leiden 2 has only 18 mapped cells, Leiden 3 has partial plasma-cell support, and Leiden 6 is unmapped (Supplementary Tables S25-S27 and Supplementary Figures S16 and S17).

Removing Leiden 3 gives rank ρ = 0.958 and top-30 Jaccard 0.714; removing Leiden 2/3 gives ρ = 0.901 and Jaccard 0.500. Restricting to exact-atlas cells gives ρ = 0.612 and Jaccard 0.071. Individual-patient rankings are heterogeneous (median pooled-rank ρ = −0.085; top-15 Jaccard 0.364), while leave-one-patient-out rankings have median ρ = 0.963 and top-15 Jaccard 0.938; median *I*² is 0.863 (Supplementary Tables S28 and S29 / Supplementary Figure S18). GSE210152 is reported as a pooled Leiden-resolved analysis with partial atlas context; cell-type mechanistic and patient-stratified inference was not performed.

### Comparator coverage differs by native input and endpoint

The primary ordering comparison contains complete, cell-aligned outputs for LineaLinkR and 11 comparators on the same 19 datasets (Supplementary Table S2). Five registered scale extensions are reported separately, while the other 13 registry datasets support native-endpoint analyses and biological case studies, including the static-versus-temporal CCC comparison summarized in Table 2 and Supplementary Table S14.

---

## Discussion

LineaLinkR links branch-aware ordering to lineage-resolved temporal CCC in snapshot scRNA-seq. Across the complete ordering core, TSCAN and MultistageOT had higher mean Spearman, whereas LineaLinkR returned the temporal and lineage coordinates used for downstream LR analysis. Across the four evaluated two-terminal simulations, branch-recovery point estimates were higher than for slingshot under the one-to-one metric. CRC, muscle and clone analyses identified recurrent expression programs and non-random structure.

The operational role scores separate pre-fork sender-leading, post-fork sender-leading, post-fork receiver-leading, opposing-lineage pre-fork and residual timing patterns. Their low agreement with prespecified FGF, Wnt and Notch priors does not support calibrated biological role assignment; they are retained as exploratory summaries of bifurcation-relative expression timing.

Temporal LR scores measure coordinated ligand and receptor expression along inferred trajectories, not signal transmission or causality. Root direction can reverse recorded time, as observed in GSE274389 and CellTag, so analyses without independent root information require direction-aware interpretation. Branch comparisons are limited to four two-terminal simulations, and the fork-evidence simulation begins from LR-score arrays rather than simulated expression. Datasets without biological replicates support descriptive profiles rather than sample-level inference. GSE210152 therefore provides pooled Leiden-level rankings.

LineaLinkR returns lineage×bin activity and bifurcation-relative labels (Table 3; Supplementary Table S30). The comparison methods address different designs: TICCI returns trajectory graphs [35], TimeTalk models LR–TF direction in time courses [36], and TraSig tests profiles on an imported trajectory [38].

---

## Conclusion

LineaLinkR connects branch-aware ordering with lineage-resolved temporal CCC in snapshot scRNA-seq. Across four two-terminal simulations, branch-recovery point estimates were higher than for slingshot under one-to-one matching, while root direction varied across datasets. The external analyses identified recurrent expression programs and clonal structure.

## Key Points

- LineaLinkR connects branch-aware ordering with lineage×bin LR activity and bifurcation-relative roles in snapshot scRNA-seq.
- Across four two-terminal simulations, one-to-one branch macro F1 was 0.869 versus 0.564 for slingshot and lineage ARI was 0.800 versus 0.105; a separate module-level fork-evidence simulation achieved AUROC 0.975.
- Withheld axes separate concordant ordering from direction failure: pancreas agrees strongly with RNA velocity; gastrulation shows weak blind-root and strong user-root stage concordance, whereas GSE274389 and CellTag reveal blind-direction reversals.
- Matched colorectal-cancer data reproduce six prespecified expression axes, and myofibroblasts increase within the stromal compartment in all ten matched patients.

---

## Acknowledgements

None

---

## Funding

This work was supported by the Fundamental Research Funds for the Central Universities (DUT22YG107), the National Natural Science Foundation of China (81600370), and the China Postdoctoral Science Foundation (2018M640270).

---

## Conflict of interest

The authors declare no competing interests.

---

## Author contributions

H.J., Y.X.H. and Y.S.H.: Conceptualization, Methodology and Software; H.J., X.Z. and H.W.: Data curation and Writing - Original Draft; H.J. and Y.X.H.: Visualization and Investigation; Y.X.H. and Y.S.H.: Supervision; H.J. and H.W.: Software and Validation; H.J., H.W. and Y.S.H.: Writing - Review and Editing. All authors read and approved the final manuscript.

---

## Data availability

Source accessions, original publications, persistent identifiers and official repository URLs for all analysed dataset objects are listed in Supplementary Table S1 and the accompanying source register. The benchmark registry contains 37 objects; GSE274389, GSE216606, GSE143437 and GSE291091 are used separately for external validation or biological analyses. GSE171524 is retained in the source register but is not used for a reported endpoint. Access to source data follows the conditions of the original repositories.

## Code availability

The LineaLinkR v0.6.4 source package, benchmark scripts, dataset provenance records and tabular results are available at https://github.com/aeqby770/manuscript-materials. Protocol details are provided in the Supplementary Methods and Table S1. LineaLinkR is released under the MIT license.

## Supplementary material

The accompanying supplement contains the Supplementary Methods, Supplementary Tables S1-S30, Supplementary Notes 1-6 and Supplementary Figures S1-S18.

---

## Table legends and tables

**Table 1.** Benchmark summary for the prespecified protocol (dataset-level means).

| Algorithm | Spearman all (n) | Spearman real (n) | Spearman sim (n) | Lineage ARI (n) | Branch F1 (n) |
|-----------|------------------|-------------------|------------------|-----------------|---------------|
| TSCAN | 0.743 (19) | 0.343 (3) | 0.818 (16) | — | — |
| MultistageOT | 0.700 (19) | 0.517 (3) | 0.734 (16) | — | — |
| **LineaLinkR** | **0.615 (19)** | **0.217 (3)** | **0.690 (16)** | **0.800 (4)** | **0.869 (4)** |
| scTour | 0.524 (19) | −0.006 (3) | 0.623 (16) | — | — |
| scFates | 0.504 (19) | 0.068 (3) | 0.586 (16) | — | — |
| monocle3 | 0.504 (19) | −0.136 (3) | 0.623 (16) | — | — |
| SCORPIUS | 0.456 (19) | 0.498 (3) | 0.448 (16) | — | — |
| slingshot | 0.438 (19) | 0.165 (3) | 0.489 (16) | 0.105 (4) | 0.564 (4) |
| Palantir | 0.388 (19) | 0.251 (3) | 0.414 (16) | — | — |
| pCreode | 0.332 (19) | 0.099 (3) | 0.376 (16) | — | — |
| DPT | 0.142 (19) | 0.224 (3) | 0.127 (16) | — | — |
| StaVia | 0.098 (19) | 0.157 (3) | 0.087 (16) | — | — |

Spearman uses the complete 19-dataset core (3 real + 16 simulated) and is conditional on returned cells after the prespecified alignment rule. ARI and branch F1 use terminal cells from four two-terminal simulations. Root-free comparator orientation details and signed values are in Supplementary Table S3. Dataset-bootstrap 95% CIs are reported in Supplementary Table S8; many-to-one F1 is reported as a sensitivity analysis.

**Table 2.** Static support among top-decile broad-lineage temporal candidates.

| Dataset | Temporal candidates | CellChat | LIANA | Both | Temporal only |
|---------|---------------------|----------|-------|------|---------------|
| pancreas | 31 | 27 | 14 | 14 | 4 |
| dentate gyrus | 26 | 24 | 20 | 20 | 2 |

**Table 3.** Native scope of LineaLinkR and temporal CCC alternatives.

| Method | Required input | Trajectory context | Native output and statistics | Branch-role support |
|--------|----------------|--------------------|------------------------------|---------------------|
| LineaLinkR | Snapshot scRNA-seq; optional cell labels | Blind root or independent user root; soft lineages and branch bins | Source×target×lineage×bin LR activity; descriptive evidence flags, no causal-lag test | Five operational labels; exploratory pathway-prior comparison |
| TICCI | Snapshot scRNA-seq | Internal Louvain/PAGA graph and entropy-directed lineage tree | CCI-weighted trajectory graph; entropy direction, no bin-resolved LR ranking or LR-lag test | No |
| TraSig | Expression, external pseudotime and trajectory edges | Sliding windows along supplied trajectory edges | LR-profile products and permutation support; synchronous profile matching | No |
| TimeTalk | Repeated time-course scRNA-seq | Pseudotime before LR/TF analysis; no bifurcation-relative context | LR scores linked to temporal TFs; Granger direction and correlation filters | No; forward, backward and feedback describe LR–TF direction |

Feature presence does not imply equivalent accuracy. This table compares native method scope and required inputs, not pooled performance; feature definitions follow the primary papers [35, 36, 38].

---

## Figure legends

**Figure 1. LineaLinkR graphical overview.** The graphical abstract follows the analysis from snapshot expression and ligand–receptor priors through expression-only trajectory reconstruction, lineage-weighted equal-mass temporal LR scoring, entropy-based fork interpretation and ranked programs. The lower ribbon reports 37 registered datasets, 24 ordering datasets (19 complete-core plus five scale extensions), four terminal branch-recovery sets, 11 temporal-inference comparators and 12 temporal bins; the remaining branch-labelled simulations are ordering-stress tests. Soft-lineage uncertainty is shown in Supplementary Figures S1 and S2, program-level uncertainty and lag examples in Supplementary Figures S5 and S9, and repeated full-workflow runtime and memory in Supplementary Figure S12.

Alt text: Wide five-stage graphical overview for LineaLinkR. The panels are labelled input, trajectory, scoring, interpretation and evidence, with compact diagrams for snapshot expression, trajectory reconstruction, temporal ligand–receptor scoring, fork roles and evidence calibration. A benchmark ribbon below the workflow reports the 37 registered datasets, 24 ordering datasets, four terminal branch-recovery sets, 11 temporal-inference comparators and 12 bins.

**Figure 2. Ordering accuracy, returned-cell completeness and native order on shared embeddings.** (A) Twelve-method scorecard for the 19-dataset common core. Horizontal bars show mean Spearman; coloured circles show Spearman, Pearson, Kendall and inverted RMSE; coverage segments span the minimum to median returned-cell fraction; and the final column gives the number of complete outputs across the 24 ordering datasets. (B, C) Pancreas and dentate gyrus cell states followed by native TSCAN, Slingshot and LineaLinkR order on identical UMAP coordinates. Dataset-level Spearman is shown for the dentate-gyrus methods.

Alt text: A 12-method scorecard combines mean-Spearman bars, four metric circles, returned-cell coverage ranges and complete-output counts. Two four-panel rows compare cell states, TSCAN, Slingshot and LineaLinkR on shared pancreas and dentate-gyrus embeddings using one common pseudotime colour scale.

**Figure 3. Branch recovery, ordering stress, output completeness and fork-evidence recovery.** (A) The four terminal-cell branch datasets and their aligned terminal-cell counts. (B, C) Dataset-paired LineaLinkR and Slingshot one-to-one branch F1 and lineage ARI, with dataset-bootstrap confidence intervals and exact sign-flip tests summarized in the subtitles. (D) Spearman ordering scores for three- and four-branch, batch and dropout simulations; these datasets are not scored as branch-label benchmarks. (E) Mean Spearman with 95% dataset-bootstrap intervals across the 19-dataset core. (F) Complete outputs and available artifacts across the 24 ordering datasets. (G) Detection rates across 600 module-level LR-score simulations with expected fork signals or null structures; the integrated AUROC is reported in the subtitle. Repeated workflow runtime and memory remain in Supplementary Figure S12.

Alt text: Seven aligned panels show terminal-cell dataset sizes, paired branch F1, paired lineage ARI, an ordering-stress heatmap, a 12-method Spearman forest plot, output completeness and fork-evidence detection rates. Only the four two-terminal datasets contribute branch F1 and ARI; panel G uses 600 separate module-level LR-score simulations.

**Figure 4. Temporal communication programs and orthogonal support.** (A) Lineage-resolved temporal-bin profiles for four leading pancreas and four leading dentate-gyrus programs. (B) Temporal score, CellChat support and LIANA support for the same programs. (C) Candidate counts by communication family. (D) Exact static-support combinations among temporal candidates. (E) Dominant operational roles within prespecified FGF, Wnt and Notch timing classes. (F) Median operational role weights for the same pathway classes. (G) Exploratory expected-role agreement, with point area showing the number of trajectories and colour showing the median expected-role margin.

Alt text: Seven panels connect temporal dynamics to external communication support and operational roles. Eight small-multiple line plots show lineage-specific temporal programs; compact matrices and bars show CellChat and LIANA support, candidate families and exact support combinations. Three pancreas pathway-role panels show role shares, median weights and expected-role agreement.

**Figure 5. Branch-resolved temporal communication during pancreas development.** (A) Cell-state geometry with the inferred lineage skeleton. (B) LineaLinkR pseudotime on the same embedding. (C) Cell-level probabilities for five soft lineages. (D) Mean lineage probability across annotated states. (E) An evidence rail for leading programs combines temporal score, evaluable fork contrast, CellChat support and LIANA support. (F) Eight expression-supported ligand–receptor programs across 12 temporal bins and five lineages. (G) Lineage-specific trends for four exemplar programs. Shared coordinates and lineage colours link cell state, temporal progression and communication dynamics across panels.

Alt text: Seven-panel pancreas-development figure. The top rows show annotated cell states, pseudotime and five lineage-probability maps on a shared embedding. A state-by-lineage matrix and four-column evidence rail precede an eight-program temporal heatmap and four lineage-specific trend panels. Fork evidence is aggregated from 1,318 evaluable rows.

**Figure 6. Branch-resolved temporal communication during dentate-gyrus neurogenesis.** (A) Neurogenic cell states with the inferred lineage skeleton. (B) LineaLinkR pseudotime on the same embedding. (C) Cell-level probabilities for five soft lineages. (D) Mean lineage probability across neurogenic states. (E) An evidence rail for leading programs combines mean and peak temporal score, bin and lineage coverage, CellChat support and LIANA support. (F) Ten expression-supported ligand–receptor programs across 12 temporal bins and five lineages. (G) Lineage-specific trends for four exemplar programs. Figure 6 uses the same panel geometry and colour semantics as Fig. 5.

Alt text: Seven-panel dentate-gyrus neurogenesis figure matching the pancreas layout. The top rows show neurogenic states, pseudotime and five lineage-probability maps. A state-by-lineage matrix and six-column evidence rail precede a ten-program temporal heatmap and four lineage-specific trend panels.

**Figure 7. Blind trajectory inference evaluated against withheld time, signal history and clonal fate.** (A) Two validation designs keep experimental time, Wnt-recorder history and CellTag clone/fate annotations outside trajectory inference; the GSE274389 lane also summarizes 12 withheld Wnt endpoints. (B) Blind LineaLinkR pseudotime on the GSE274389 gastruloid embedding. (C) Recorded Wnt history and endpoint Wnt-target state among 120-h cells on the same coordinates, paired with a four-lineage by 12-bin validation matrix. (D) Blind pseudotime on the annotation-free GSE216606 CellTag-multi embedding. (E) Published day-5 fates revealed after fitting. (F) Eligible clone centroids, with colour and point area encoding day-5 clone size. (G) Published fates distributed across inferred soft lineages. (H) Expansion of the largest eligible clones from day 2.5 to day 5. Both blind orderings have the opposite orientation to recorded time; the CellTag analysis also shows non-random clonal concentration.

Alt text: Eight-panel external-validation figure linking two blinded designs. A paired timeline shows which experimental variables were withheld. Shared gastruloid coordinates compare blind pseudotime with recorded Wnt history and endpoint target state. Shared CellTag coordinates compare blind pseudotime, published fate and clone expansion. A fate-by-lineage heatmap and clone-growth trajectories show non-random clonal structure despite reversed temporal direction.

**Figure 8. Independent biological support for muscle-regeneration and colorectal-cancer communication programs.** (A) Two validation lanes connect an injury time course to FAP depletion and a matched normal–tumour cohort to patient-level composition and expression analyses. (B) Sample-level C3–C3AR1 expression scores across recorded injury days 0, 2, 5 and 7 in GSE143437. (C) Mesenchymal/FAP-like pseudobulk ligand effects during injury and after FAP depletion in GSE291091. (D) Hallmark pathway responses in broad mononuclear phagocytes after depletion. (E) Patient-by-axis heatmap of paired symmetric expression changes for six prespecified colorectal-cancer LR programs in ten matched patients. (F) Mean paired change and bootstrap interval for the same axes. (G) Myofibroblast and SPP1-positive-myeloid fractions across 23 tumour patients. The muscle panels report phase specificity and a system-level FAP-depletion response; the colorectal panels report patient-level composition and expression recurrence.

Alt text: Seven-panel biological validation figure. Panel A shows parallel muscle and colorectal-cancer validation designs. Panels B–D show ten sample-level muscle scores, ligand effects after injury and FAP depletion, and pathway responses in phagocytes. Panel E is a dominant heatmap of paired normal-to-tumour changes for six ligand–receptor axes across ten matched patients. Panel F summarizes those paired changes with bootstrap intervals. Panel G compares myofibroblast and SPP1-positive-myeloid fractions across tumour patients.

---

## References

1. Macosko EZ, Basu A, Satija R, et al. Highly parallel genome-wide expression profiling of individual cells using nanoliter droplets. *Cell*. 2015;161(5):1202–1214. doi:10.1016/j.cell.2015.05.002
2. Zheng GXY, Terry JM, Belgrader P, et al. Massively parallel digital transcriptional profiling of single cells. *Nat Commun*. 2017;8:14049. doi:10.1038/ncomms14049
3. Street K, Risso D, Fletcher RB, et al. Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics. *BMC Genomics*. 2018;19:477. doi:10.1186/s12864-018-4772-0
4. Cao J, Spielmann M, Qiu X, et al. The single-cell transcriptional landscape of mammalian organogenesis. *Nature*. 2019;566:496–502. doi:10.1038/s41586-019-0969-x
5. Wolf FA, Hamey FK, Plass M, et al. PAGA: graph abstraction reconciles clustering with trajectory inference through a topology preserving map of single cells. *Genome Biol*. 2019;20:59. doi:10.1186/s13059-019-1663-x
6. Setty M, Kiseliovas V, Levine J, et al. Characterization of cell fate probabilities in single-cell data with Palantir. *Nat Biotechnol*. 2019;37:451–460. doi:10.1038/s41587-019-0068-4
7. Haghverdi L, Büttner M, Wolf FA, et al. Diffusion pseudotime robustly reconstructs lineage branching. *Nat Methods*. 2016;13:845–848. doi:10.1038/nmeth.3971
8. Herman JS, Sagar, Grün D. FateID infers cell fate bias in multipotent progenitors from single-cell RNA-seq data. *Nat Methods*. 2018;15:379–386. doi:10.1038/nmeth.4662
9. Cannoodt R, Saelens W, Sichien D, et al. SCORPIUS improves trajectory inference and identifies novel modules in dendritic cell development. *bioRxiv*. 2016. doi:10.1101/079509
10. Ji Z, Ji H. TSCAN: pseudo-time reconstruction and evaluation in single-cell RNA-seq analysis. *Nucleic Acids Res*. 2016;44:e117. doi:10.1093/nar/gkw430
11. Faure L, Soldatov R, Kharchenko PV, Adameyko I. scFates: a scalable python package for advanced pseudotime and bifurcation analysis from single-cell data. *Bioinformatics*. 2023;39(1):btac746. doi:10.1093/bioinformatics/btac746
12. Stassen SV, Kobashi M, Lam EY, et al. StaVia: spatially and temporally aware cartography with higher-order random walks for cell atlases. *Genome Biol*. 2024;25:224. doi:10.1186/s13059-024-03347-y
13. Li Q. scTour: a deep learning architecture for robust inference and accurate prediction of cellular dynamics. *Genome Biol*. 2023;24:149. doi:10.1186/s13059-023-02988-9
14. Tronstad M, Karlsson J, Dahlin JS. MultistageOT: multistage optimal transport infers trajectories from a snapshot of single-cell data. *Proc Natl Acad Sci USA*. 2025;122(50):e2516046122. doi:10.1073/pnas.2516046122
15. Trapnell C, Cacchiarelli D, Grimsby J, et al. The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells. *Nat Biotechnol*. 2014;32(4):381–386. doi:10.1038/nbt.2859
16. Luecken MD, Theis FJ. Current best practices in single-cell RNA-seq analysis: a tutorial. *Mol Syst Biol*. 2019;15(6):e8746. doi:10.15252/msb.20188746
17. Wang J, Sun L, Wei N, et al. A semi-supervised Bayesian approach for marker gene trajectory inference from single-cell RNA-seq data. *Bioinformatics*. 2025;41(9):btaf454. doi:10.1093/bioinformatics/btaf454
18. Zhang H, Zhang W, Zheng X, et al. scRDEN: single-cell dynamic gene rank differential expression network and robust trajectory inference. *Sci Rep*. 2025;15:16963. doi:10.1038/s41598-025-01969-1
19. Wang Y, Lin Z, Wang T. InterVelo: a mutually enhancing model for estimating pseudotime and RNA velocity in multi-omic single-cell data. *Bioinformatics*. 2025;41(10):btaf500. doi:10.1093/bioinformatics/btaf500
20. Ge M, Miao J, Qi J, et al. TIVelo: RNA velocity estimation leveraging cluster-level trajectory inference. *Nat Commun*. 2025;16:6258. doi:10.1038/s41467-025-61628-x
21. Kouadri Boudjelthia I, Milite S, El Kazwini N, et al. Interpretable learning of temporal cellular dynamics from single-cell data. *Cell Rep Methods*. 2026;6(3):101342. doi:10.1016/j.crmeth.2026.101342
22. Fang M, Gorin G, Pachter L. Trajectory inference from single-cell genomics data with a process time model. *PLoS Comput Biol*. 2025;21(1):e1012752. doi:10.1371/journal.pcbi.1012752
23. Cheng M, Jansen J, Reimer KC, et al. PHLOWER leverages single-cell multimodal data to infer complex, multi-branching cell differentiation trajectories. *Nat Methods*. 2025;22(11):2328–2336. doi:10.1038/s41592-025-02870-5
24. Huizing GJ, Samaran J, Capocefalo D, et al. STORIES: learning cell fate landscapes from spatial transcriptomics using optimal transport. *Nat Methods*. 2026;23(3):522–531. doi:10.1038/s41592-025-02855-4
25. Leary JR, Dong X, Bacher R. Interpretable trajectory inference with single-cell linear adaptive negative-binomial expression (scLANE) testing. *Nucleic Acids Res*. 2026;54(2):gkaf1494. doi:10.1093/nar/gkaf1494
26. Jin S, Guerrero-Juarez CF, Zhang L, et al. Inference and analysis of cell-cell communication using CellChat. *Nat Commun*. 2021;12:1088. doi:10.1038/s41467-021-21246-9
27. Vento-Tormo R, Efremova M, Botting RA, et al. Single-cell reconstruction of the early maternal–fetal interface in humans. *Nature*. 2018;563:347–353. doi:10.1038/s41586-018-0698-6
28. Troulé K, Petryszak R, Cakir B, et al. CellPhoneDB v5: inferring cell–cell communication from single-cell multiomics data. *Nat Protoc*. 2025;20:3412–3440. doi:10.1038/s41596-024-01137-1
29. Dimitrov D, Türei D, Garrido-Rodriguez M, et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. *Nat Commun*. 2022;13:3224. doi:10.1038/s41467-022-30755-0
30. Hou R, Denisenko E, Ong HT, et al. Predicting cell-to-cell communication networks using NATMI. *Nat Commun*. 2020;11:5011. doi:10.1038/s41467-020-18873-z
31. Browaeys R, Saelens W, Saeys Y. NicheNet: modeling intercellular communication by linking ligands to target genes. *Nat Methods*. 2020;17:159–162. doi:10.1038/s41592-019-0667-5
32. Cabello-Aguilar R, Laugraud A, Colinge J, et al. SingleCellSignalR: inference of intercellular networks from single-cell transcriptomics. *Nucleic Acids Res*. 2020;48:e55. doi:10.1093/nar/gkaa183
33. Ramilowski JA, Goldberg T, Harshbarger J, et al. A draft network of ligand–receptor-mediated multicellular signalling in human. *Nat Commun*. 2015;6:7866. doi:10.1038/ncomms8866
34. Armingol E, Officer A, Harismendy O, Lewis NE. Deciphering cell–cell interactions and communication from gene expression. *Nat Rev Genet*. 2021;22(2):71–88. doi:10.1038/s41576-020-00292-x
35. Fu Y, Qu H, Qu D, Zhao M. Trajectory Inference with Cell–Cell Interactions (TICCI): intercellular communication improves the accuracy of trajectory inference methods. *Bioinformatics*. 2025;41(2):btaf027. doi:10.1093/bioinformatics/btaf027
36. Wang L, Zheng Y, Sun Y, et al. TimeTalk uses single-cell RNA-seq datasets to decipher cell-cell communication during early embryo development. *Commun Biol*. 2023;6:901. doi:10.1038/s42003-023-05283-2
37. Granger CWJ. Investigating causal relations by econometric models and cross-spectral methods. *Econometrica*. 1969;37(3):424–438. doi:10.2307/1912791
38. Li D, Velazquez JJ, Ding J, et al. TraSig: inferring cell-cell interactions from pseudotime ordering of scRNA-Seq data. *Genome Biol*. 2022;23:73. doi:10.1186/s13059-022-02629-7
39. Bastidas-Ponce A, Tritschler S, Dony L, et al. Comprehensive single cell mRNA profiling reveals a detailed roadmap for pancreatic endocrinogenesis. *Development*. 2019;146(12):dev173849. doi:10.1242/dev.173849
40. Hochgerner H, Zeisel A, Lönnerberg P, Linnarsson S. Conserved properties of dentate gyrus neurogenesis across postnatal development revealed by single-cell RNA sequencing. *Nat Neurosci*. 2018;21(2):290–299. doi:10.1038/s41593-017-0056-2
41. Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol*. 2018;19(1):15. doi:10.1186/s13059-017-1382-0
42. Shannon CE. A mathematical theory of communication. *Bell Syst Tech J*. 1948;27(3):379–423. doi:10.1002/j.1538-7305.1948.tb01338.x
43. Demšar J. Statistical comparisons of classifiers over multiple data sets. *J Mach Learn Res*. 2006;7:1–30. https://jmlr.org/papers/v7/demsar06a.html
44. McNamara HM, Solley SC, Adamson B, et al. Recording morphogen signals reveals mechanisms underlying gastruloid symmetry breaking. *Nat Cell Biol*. 2024;26:1832–1844. doi:10.1038/s41556-024-01521-9
45. Jindal K, Adil MT, Yamaguchi N, et al. Single-cell lineage capture across genomic modalities with CellTag-multi reveals fate-specific gene regulatory changes. *Nat Biotechnol*. 2024;42:946–959. doi:10.1038/s41587-023-01931-4
46. De Micheli AJ, Laurilliard EJ, Heinke CL, et al. Single-cell analysis of the muscle stem cell hierarchy identifies heterotypic communication signals involved in skeletal muscle regeneration. *Cell Rep*. 2020;30(10):3583–3595.e5. doi:10.1016/j.celrep.2020.02.067
47. Liu X, Perez Carbajal EE, Hwang YC, et al. Mesenchymal stromal cell-mediated intercellular communication: mapping the interactome for skeletal muscle homeostasis and regeneration. *Adv Sci (Weinh)*. 2026:e07444. doi:10.1002/advs.202507444
48. Lee HO, Hong Y, Etlioglu HE, et al. Lineage-dependent gene expression programs influence the immune landscape of colorectal cancer. *Nat Genet*. 2020;52:594–603. doi:10.1038/s41588-020-0636-z
49. Qi J, Sun H, Zhang Y, et al. Single-cell and spatial analysis reveal interaction of FAP-positive fibroblasts and SPP1-positive macrophages in colorectal cancer. *Nat Commun*. 2022;13:1742. doi:10.1038/s41467-022-29366-6
50. Saelens W, Cannoodt R, Todorov H, Saeys Y. A comparison of single-cell trajectory inference methods. *Nat Biotechnol*. 2019;37:547–554. doi:10.1038/s41587-019-0071-9
51. Bergen V, Lange M, Peidli S, Wolf FA, Theis FJ. Generalizing RNA velocity to transient cell states through dynamical modeling. *Nat Biotechnol*. 2020;38(12):1408–1414. doi:10.1038/s41587-020-0591-3
52. Pijuan-Sala B, Griffiths JA, Guibentif C, et al. A single-cell molecular map of mouse gastrulation and early organogenesis. *Nature*. 2019;566:490–495. doi:10.1038/s41586-019-0933-9
53. Traeuble K, Munz M, Pauli J, et al. Integrated single-cell atlas of human atherosclerotic plaques. *Nat Commun*. 2025;16:8255. doi:10.1038/s41467-025-63202-x
54. Cang Z, Zhao Y, Almet AA, et al. Screening cell–cell communication in spatial transcriptomics via collective optimal transport. *Nat Methods*. 2023;20(2):218–228. doi:10.1038/s41592-022-01728-4
55. Shao X, Li C, Yang H, et al. Knowledge-graph-based cell-cell communication inference for spatially resolved transcriptomic data with SpaTalk. *Nat Commun*. 2022;13:4429. doi:10.1038/s41467-022-32111-8
56. Liu B, Mosienko V, Vaccari Cardoso B, et al. Glio- and neuro-protection by prosaposin is mediated by orphan G-protein coupled receptors GPR37L1 and GPR37. *Glia*. 2018;66(11):2414–2426. doi:10.1002/glia.23480
