# Supplementary Material: LineaLinkR

**Manuscript:** LineaLinkR: lineage-resolved temporal ligand–receptor mapping along branched single-cell trajectories  
**Journal:** Briefings in Bioinformatics (Problem solving protocol)  

---

## Supplementary Methods

### Benchmark input and direction protocol

All algorithms receive the same library-normalized, `log1p` expression matrix and cell identifiers unless their documented interface requires another representation. Ground-truth order, stage and branch labels are withheld from model fitting. Adapter outputs retain the native cell subset; missing predictions are not imputed. Root choice, direction policy and returned-cell fraction are reported per dataset. Slingshot, SCORPIUS and root-free TSCAN longest-path outputs may be reversed only under the prespecified signed-Spearman rule; the signed value and flip flag remain available. LineaLinkR retains its blind direction. Branch metrics use only `simulated_branch`, `sim_2br_500c`, `sim_2br_2000c` and `sim_2br_5000c`, restricted to cells with `true_pseudotime ≥ 0.4`; the other 12 simulated datasets are ordering-only stress tests.

### Expression-only attractor endpoint selection

The expression graph uses 2000 highly variable genes, 30 PCs and 50 neighbours. Candidate terminal attractors are peripheral, low-potency states. A double-precision Rust kernel estimates attractor absorption, after which five min–max-scaled cell features are combined: detected-gene potency (0.25), fate entropy (0.20), nearest-attractor distance (0.30), inverse-distance-weighted source flow (0.15) and one minus maximum commitment (0.10). The reported `detected_genes_ordering` protocol evaluates high- and low-detected-gene potency with both state-aggregate and cellwise-component anchors. It orients each candidate by high detected-gene potency and selects among them using blind alignment and ordering gains. If the maximum blind ordering score is below 0.45, a confidence-regularized PC/k ensemble supplies the selected endpoint; otherwise the perturbation grid reports endpoint stability.

### Graph, lineage and temporal LR details

Locally scaled fuzzy memberships convert neighbour distances into edge costs. The graph is symmetrized; disconnected components are joined by recorded geometric bridges, while cluster-MST edge augmentation is disabled in the reported analysis. Dijkstra distance is scaled to [0,1]. Root-to-leaf cluster-centroid paths define candidates, while cell distances remain on the weighted graph. Gaussian path weights are row-normalized by Rust log-sum-exp and describe path proximity rather than calibrated fate probability.

Counts are scaled to 10,000 per cell and log-transformed. Each lineage is divided into 12 equal-lineage-weight pseudotime quantiles. Ligand and receptor complex means use the minimum subunit value and require detection ≥0.05; the smaller source/target Kish effective size is retained. The Hill scale is the dataset-wide mean positive gene-by-state expression. Static scoring uses the same expression and detection rules without bins. The full temporal grid retains zero and non-evaluable rows, detection, effective size, complex membership, source, target, lineage and bin.

### Fork, role and resampling details

Normalized lineage entropy is smoothed over three bins. A fork requires two consecutive bins below 0.5; insufficient or flat trajectories receive no fork rather than a midpoint. Focal post-minus-pre activity is contrasted with the effective-size-weighted sibling change using two bins per side. Sample-aware analyses recompute activity in fixed pooled bins and bootstrap biological samples; otherwise the contrast is descriptive. Role components are evaluated only after expression, effective-size and fork-boundary gates. Six component scores are exported. The five-class calibration assigns instructive, reinforcing, feedback, lateral or bystander; the auxiliary transition score is retained separately.

The module-level fork-evidence simulation supplies synthetic lineage×bin LR scores directly to the fork-contrast calculation. Six scenarios contain 100 replicates each: shared trend, focal branch, donor-confounded null, low-expression spike, 30% dropout and nested fork. Focal branch, dropout and nested fork are positive scenarios; the other three are null structures. Evidence combines the non-negative pooled sibling contrast, median sample contrast and sample support, with detection at 0.25. Reported endpoints are integrated AUROC, sensitivity, specificity and scenario-level detection rate. Because expression generation and LR-pair ranking are outside this simulation, it evaluates the fork-evidence layer rather than end-to-end LR recovery.

### External analyses and runtime

GSE216606 inference uses an annotation-free H5AD expression matrix and neutral cell identifiers. Day, clone, fate, author state and technical batch are joined after fitting. Clone concentration compares within-clone lineage-pair agreement with the global lineage-pair probability and uses clone bootstrap intervals. GSE274389 time and Wnt-recorder labels are likewise withheld until trajectory fitting is complete.

For GSE143437, recorded injury day is used as the biological time axis after blind LineaLinkR pseudotime is compared with day. C3-C3AR1 expression scores are summarized by biological sample at days 0, 2, 5 and 7. GSE291091 uses two biological replicates per group; differential expression is estimated from pseudobulk counts within mesenchymal/FAP-like cells and broad mononuclear phagocytes, followed by Hallmark rank enrichment and sample-level composition summaries.

For GSE132465, six fibroblast-macrophage expression axes are fixed before analysis. Ten matched normal-tumour patients form the paired inference set. Paired symmetric changes, Wilcoxon tests, Benjamini-Hochberg adjustment and patient-bootstrap intervals use patient as the unit. Tumour-only subtype fractions and expression proxies are summarized across 23 patients, with fine-subtype proxies restricted to patients meeting the minimum cell-count requirement.

GSE210152 is analysed as a pooled Leiden-resolved dataset. Exact barcode joining to the processed plaque atlas supplies partial annotation without expression-based imputation. Temporal LR rankings are recomputed after prespecified cell restrictions and by patient or leave-one-patient-out subsets while retaining the pooled pseudotime and lineage weights. Runtime measurements use fresh processes and include input loading, preprocessing, graph construction, trajectory inference, static LR and temporal/fork scoring; optional taxonomy calibration is excluded. Process-tree peak RSS includes descendants.

---

## Supplementary Tables

### Supplementary Table S1. Benchmark and external-analysis dataset manifest

#### Benchmark registry (37 datasets)

| Dataset | Registry group | Species | Ordering role | Notes |
|---------|----------------|---------|---------------|-------|
| lung | curated | mouse | scale extension | Bleomycin-injured lung time course |
| paul15 | curated | mouse | non-ordering | Myeloid progenitors |
| pancreas | curated | mouse | non-ordering | Embryonic pancreas case study |
| dentategyrus | curated | mouse | complete core | Postnatal neurogenesis |
| cellrank_stemnet | curated | mouse | non-ordering | Derived pancreas STEMNET object |
| cellagentchat_scc | curated | human | non-ordering | Cutaneous squamous-cell carcinoma sample |
| monocle3_worm_l2 | curated | *C. elegans* | non-ordering | L2 atlas |
| gastrulation_erythroid | literature | mouse | complete core | Gastrulation erythroid lineage |
| zebrafish_urd | literature | zebrafish | non-ordering | Embryogenesis URD dataset |
| zebrafish_neural_dev | literature | zebrafish | complete core | pMN neural progenitors |
| celegans_embryo | literature | *C. elegans* | non-ordering | Embryogenesis atlas |
| human_cd34 | literature | human | non-ordering | CD34-positive bone-marrow hematopoiesis |
| larry_invitro | literature | mouse | scale extension | In-vitro lineage tracing |
| simulated_branch | project generated | simulated | complete core | Two-terminal branch benchmark; 1,357 terminal cells scored |
| sim_2br_500c | project generated | simulated | complete core | Two-terminal branch benchmark; 366 terminal cells scored |
| sim_2br_2000c | project generated | simulated | complete core | Two-terminal branch benchmark; 1,470 terminal cells scored |
| sim_2br_5000c | project generated | simulated | complete core | Two-terminal branch benchmark; 3,684 terminal cells scored |
| sim_3br_500c | project generated | simulated | complete core | Ordering-only segment stress, 500 cells |
| sim_3br_2000c | project generated | simulated | complete core | Ordering-only segment stress, 2,000 cells |
| sim_3br_5000c | project generated | simulated | complete core | Ordering-only segment stress, 5,000 cells |
| sim_4br_500c | project generated | simulated | complete core | Ordering-only segment stress, 500 cells |
| sim_4br_2000c | project generated | simulated | complete core | Ordering-only segment stress, 2,000 cells |
| sim_4br_5000c | project generated | simulated | complete core | Ordering-only segment stress, 5,000 cells |
| sim_3br_2000c_batch | project generated | simulated | complete core | Ordering-only segment stress with batch effect |
| sim_3br_2000c_drop30 | project generated | simulated | complete core | Ordering-only segment stress with 30% dropout |
| sim_3br_2000c_drop50 | project generated | simulated | complete core | Ordering-only segment stress with 50% dropout |
| sim_3br_2000c_drop70 | project generated | simulated | complete core | Ordering-only segment stress with 70% dropout |
| sim_3br_2000c_drop90 | project generated | simulated | complete core | Ordering-only segment stress with 90% dropout |
| simulated_complex_3branch | project generated | simulated | complete core | Ordering-only segment stress with disconnected components |
| GSE243917_RA | GEO disease | human | non-ordering | Rheumatoid-arthritis synovium |
| GSE210152_athero | GEO disease | human | non-ordering | Atherosclerotic-plaque case study |
| GSE125527_IBD | GEO disease | human | non-ordering | Ulcerative-colitis immune atlas |
| GSE132465_CRC | GEO disease | human | non-ordering | Colorectal-cancer validation |
| GSE135893_IPF | GEO disease | human | non-ordering | Pulmonary-fibrosis atlas |
| human_fetal_liver_organoid | extension | human | scale extension | Fetal-liver organoid differentiation |
| human_fetal_pancreas | extension | human | scale extension | Developing fetal pancreas |
| zebrafish_tracerseq | extension | zebrafish | scale extension | Embryo TracerSeq atlas |

#### External validation and biological-analysis datasets (5 datasets)

| Dataset | Repository identifiers | Original study | Analysis role |
|---------|------------------------|----------------|---------------|
| GSE274389_signal_recorder | GEO GSE274389; PRJNA1146327 | McNamara et al., *Nature Cell Biology* (2024) | Withheld-time and Wnt-recorder validation |
| GSE216606_celltag_lsk | GEO GSE216606; SuperSeries GSE216521; PRJNA894497 | Jindal et al., *Nature Biotechnology* (2024) | Withheld-day, clone and fate validation |
| GSE143437_muscle_regeneration | GEO GSE143437; PRJNA600499; SRP241205 | De Micheli et al., *Cell Reports* (2020) | Regeneration time-course analysis |
| GSE291091_fap_depletion | GEO GSE291091; PRJNA1231536 | Liu et al., *Advanced Science* (2026) | FAP-depletion perturbation analysis |
| GSE171524_COVID | GEO GSE171524; PRJNA719842; SCP1219 | Melms et al., *Nature* (2021) | Source-register only; no reported endpoint |

The source register supplies the DOI or PMID, official repository URL and source class for each of the 42 objects. The five datasets in the second section are outside the 37-dataset benchmark registry.

---

### Supplementary Table S2. Ordering-method coverage in the 19-dataset complete core

| Algorithm | Complete-core outputs | Fair pairs with LineaLinkR |
|-----------|----------------------:|---------------------------:|
| LineaLinkR | 19 | — |
| slingshot | 19 | 19 |
| DPT | 19 | 19 |
| SCORPIUS | 19 | 19 |
| MultistageOT | 19 | 19 |
| scFates | 19 | 19 |
| pCreode | 19 | 19 |
| monocle3 | 19 | 19 |
| TSCAN | 19 | 19 |
| Palantir | 19 | 19 |
| StaVia | 19 | 19 |
| scTour | 19 | 19 |

The 37-dataset registry contains 24 ordering datasets and 13 native-endpoint or application datasets. Table 1 uses the complete 19-dataset core; the five scale extensions are reported separately.

---

### Supplementary Table S3. Pseudotime signed Spearman and orientation policy

Per-dataset, per-algorithm columns are `Spearman` (orientation-aligned for root-free slingshot, SCORPIUS and pCreode, and for TSCAN's longest path without a specified start endpoint, when signed ρ < 0), `Spearman_signed` (raw), `orientation_flipped`, `root_input_mode`, and `direction_policy`. TSCAN flips on 18/19 complete-core datasets and changes from signed mean −0.717 to aligned mean 0.743; LineaLinkR retains its inferred direction.

---

### Supplementary Table S4. Returned-cell completeness in the 19-dataset core

Scores are conditional on finite cells returned by each method.

| Algorithm | Datasets | Minimum completeness | Median completeness | Outputs ≥90% |
|-----------|---------:|---------------------:|--------------------:|-------------:|
| LineaLinkR | 19 | 1.000 | 1.000 | 19 |
| slingshot | 19 | 1.000 | 1.000 | 19 |
| DPT | 19 | 1.000 | 1.000 | 19 |
| SCORPIUS | 19 | 1.000 | 1.000 | 19 |
| MultistageOT | 19 | 1.000 | 1.000 | 19 |
| scFates | 19 | 1.000 | 1.000 | 19 |
| pCreode | 19 | 1.000 | 1.000 | 19 |
| monocle3 | 19 | 1.000 | 1.000 | 19 |
| TSCAN | 19 | 1.000 | 1.000 | 19 |
| Palantir | 19 | 1.000 | 1.000 | 19 |
| StaVia | 19 | 1.000 | 1.000 | 19 |
| scTour | 19 | 1.000 | 1.000 | 19 |

All 12 methods return complete, cell-aligned outputs on the primary 19-dataset core.

### Supplementary Table S5. Default LineaLinkR hyperparameters for the reported analyses

| Parameter | Value | Description |
|-----------|-------|-------------|
| Package version | 0.6.4 | LineaLinkR release used for the analyses |
| PCA components | 30 | Neighbour graph, root perturbations and lineage geometry |
| kNN *k* | 50 | PCA-space neighbours converted to locally scaled fuzzy edge costs |
| Temporal bins per lineage | 12 | Equal total lineage-weight pseudotime bins |
| Temporal minimum expression fraction | 0.05 | Required for both ligand and receptor complexes |
| MM-LR Hill scale *S* | data-driven | Dataset-level mean of positive gene-by-cell-state expression means |
| LR scale factor | 10 000 | Expression scaling before LR product |
| Root strategy | `attractor_energy` | Expression-only attractor evidence without time or branch labels |
| Root potency mode | `detected_genes_ordering` | Four candidates: high/low detected-gene potency × state/cellwise anchor |
| Root score weights | 0.25, 0.20, 0.30, 0.15, 0.10 | Potency, fate entropy, nearest-attractor distance, source flow and one minus commitment |
| Root endpoint selector | Alignment gain 0.20; ordering gain 0.10; representation gain 0.00 | Selects among four expression-only terminal-endpoint candidates |
| Root ensemble trigger | Maximum blind ordering score <0.45 | Uses the confidence-regularized PC/k root ensemble |
| Root start cell | Minimum selected oriented pseudotime | Ties resolved by high detected-gene potency and then cell index |
| Root-stability k grid | 30, 50, 70 | Neighbourhood perturbations |
| Root-stability PC grid | 10, 20, 30 | Representation perturbations |
| Root-stability resamples | 20 at 80% per cluster | Deterministic cluster-stratified resampling |
| Root-stability margin | 0.10 | Values at or below this margin indicate directional ambiguity |
| Lineage embedding | PCA | Full 30-PC geometry; 8–25 k-means refinement when lineage clusters <8 |
| Lineage normalization | log-sum-exp | Prevents Gaussian-weight underflow in high-dimensional PCA space |
| Entropy smoothing | MA(3) | Per-lineage entropy curve |
| Fork threshold | 0.5 | First two-consecutive-bin crossing after MA(3) smoothing |
| Fork window | 2 bins per side | Focal post-minus-pre change contrasted with sibling change |
| Fork activity floor | 10⁻⁴ | Lower-activity trajectories are non-evaluable |
| Sample bootstrap | 999 | Used only when biological sample labels support rescoring |
| Static LR permutation count | 100 | Used by `ComputeLRScores`; not a case-study candidate filter |
| Role reliability gate | Kish *n*eff ≥ 5 | Lower-support rows collapse to bystander = 1 |
| Reported case-study candidate rule | Top 10% mean temporal score; detected in ≥2 bins and ≥1 lineage | Then annotated for CellChat/LIANA overlap; no formal q-value claim |
| Random seed | 12345 | `set.seed` |
| Threads | `detectCores()` | `OMP_NUM_THREADS` |
| Graph-k ablation values | 10, 30, 50 (default) | Sensitivity when *n* ≤ 30 000 |

---

### Supplementary Table S6. Temporal bin-number sensitivity (pancreas)

| Number of bins | Temporal rows | Pairs | Bifurcation lineages | Top-1 pair | Jaccard top-30 vs 12 |
|---------|---------------|-------|----------------------|------------|----------------------|
| 8 | 174 680 | 313 | 2 | Mdk–Ncl | 0.875 |
| 12 (default) | 272 340 | 345 | 4 | Mdk–Ncl | 1.000 |
| 16 | 380 880 | 354 | 4 | Mdk–Ncl | 0.818 |
| 20 | 482 100 | 364 | 4 | Mdk–Ncl | 0.818 |

### Supplementary Table S7. Full-panel Hill-S multiplier sensitivity

Transform: if score = geo/(S+geo), then score under *S* → *kS* is score / (*k* − (*k*−1)·score). Sixteen of 37 datasets have non-zero scores and are evaluable.

| *k* | Evaluable datasets | Median Jaccard | IQR | Minimum |
|-----|---------------------|----------------|-----|---------|
| 0.5 | 16 | 0.935 | 0.935–1.000 | 0.875 |
| 1.0 | 16 | 1.000 | 1.000–1.000 | 1.000 |
| 2.0 | 16 | 0.935 | 0.935–1.000 | 0.875 |
| 5.0 | 16 | 0.905 | 0.875–1.000 | 0.765 |

### Supplementary Table S8. Dataset-bootstrap 95% CIs for Table 1 means

Intervals use 2000 resamples of datasets within each algorithm×metric×scope.

| Algorithm | Metric | Scope | n | Mean | 95% CI |
|-----------|--------|-------|---|------|--------|
| LineaLinkR | Spearman | all | 19 | 0.615 | [0.369, 0.782] |
| LineaLinkR | Spearman | real | 3 | 0.217 | [0.063, 0.350] |
| LineaLinkR | Spearman | sim | 16 | 0.690 | [0.447, 0.840] |
| LineaLinkR | lineage ARI | terminal branch | 4 | 0.800 | [0.400, 1.000] |
| slingshot | lineage ARI | terminal branch | 4 | 0.105 | [0.000, 0.210] |
| LineaLinkR | Branch F1, one-to-one | terminal branch | 4 | 0.869 | [0.606, 1.000] |
| slingshot | Branch F1, one-to-one | terminal branch | 4 | 0.564 | [0.405, 0.724] |
| LineaLinkR | Branch F1, many-to-one sensitivity | terminal branch | 4 | 0.928 | [0.783, 1.000] |
| slingshot | Branch F1, many-to-one sensitivity | terminal branch | 4 | 0.686 | [0.647, 0.733] |

Intervals use lineage-weighted bootstrap resampling for the leading pancreas temporal LR pairs.

### Supplementary Table S9. Dataset-paired branch metrics versus slingshot

Differences are LineaLinkR minus slingshot over four two-terminal simulations. Confidence intervals use 10 000 dataset bootstrap resamples; exact two-sided *P* values enumerate all 2^4 sign assignments.

| Metric | LLR mean | slingshot mean | Paired difference | 95% bootstrap CI | Exact sign-flip *P* | Positive / zero / negative |
|--------|---------:|----------------:|------------------:|:----------------:|----------------------:|:--------------------------:|
| One-to-one terminal-branch macro F1 | 0.869 | 0.564 | 0.304 | −0.010–0.595 | 0.250 | 3 / 0 / 1 |
| Lineage ARI | 0.800 | 0.105 | 0.695 | 0.294–1.000 | 0.125 | 4 / 0 / 0 |

### Supplementary Table S10. Module-level fork-evidence simulation

| Scenario | Expected class | Replicates | Evaluable fraction | Median evidence score | Detection rate |
|----------|----------------|-----------:|-------------------:|----------------------:|---------------:|
| focal branch | positive | 100 | 1.00 | 0.797 | 1.00 |
| dropout (30%) | positive | 100 | 1.00 | 0.561 | 0.74 |
| nested fork | positive | 100 | 1.00 | 0.799 | 1.00 |
| shared trend | null | 100 | 1.00 | 0.000 | 0.00 |
| donor-confounded null | null | 100 | 1.00 | 0.000 | 0.00 |
| low-expression spike | null | 100 | 0.00 | 0.000 | 0.00 |

Across all 600 replicates, integrated AUROC was 0.975, sensitivity was 0.913 and specificity was 1.000. Inputs are synthetic lineage×bin LR scores rather than simulated expression matrices.

---

### Supplementary Table S11. Module-output dependence check (pancreas exemplar)

| Ablation | LR pairs | Temporal rows | Fork points |
|----------|------------|-----------------|----------------------|
| full | 3785 | 272340 | 5 |
| no LR | 0 | 0 | 0 |
| no temporal binning | 3785 | 0 | 0 |
| no bifurcation | 3785 | 272340 | 0 |

Pancreas contains 3696 cells, 8 types and 5 inferred lineages. The three module removals use the same upstream trajectory object.

---

### Supplementary Table S12. Trajectory-level coarse pathway-prior calibration (pancreas)

| Pathway | Coarse prior | instructive | reinforcing | feedback | lateral | bystander | Agreement with prior |
|---------|--------------|------------:|------------:|---------:|--------:|----------:|---------------------:|
| FGF | instructive | 19 | 17 | 8 | 0 | 397 | 0.043 |
| Wnt | instructive | 31 | 13 | 0 | 0 | 213 | 0.121 |
| Notch | feedback | 79 | 21 | 7 | 0 | 442 | 0.013 |

The unit is one ligand–receptor pair × source × target × lineage trajectory after averaging role weights across structurally valid, expression-supported bins. The final column is agreement with the prespecified pathway timing class.

---

### Supplementary Table S13. Expression-supported trajectory subset of the pathway-prior calibration

| Pathway | Coarse prior | Expression-supported trajectories | Matching trajectories | Prior agreement | Median expected-role weight | Median expected-role margin |
|---------|--------------|----------------------------------:|----------------------:|----------------:|----------------------------:|----------------------------:|
| FGF | instructive | 441 | 19 | 0.043 | 0.141 | −0.644 |
| Wnt | instructive | 257 | 31 | 0.121 | 0.113 | −0.660 |
| Notch | feedback | 549 | 7 | 0.013 | 0.000 | −0.745 |

### Supplementary Table S14. Temporal LR overlap with static CCC (pancreas and dentategyrus)

| Dataset | LLR temporal pairs | CellChat | LIANA | Both | Temporal only |
|---------|-------------------|--------------------|-----------------|------|---------------------------|
| pancreas | 31 | 27 | 14 | 14 | 4 |
| dentate gyrus | 26 | 24 | 20 | 20 | 2 |

Temporal-only (operational): top-decile temporal candidates without CellChat support or top-decile LIANA support; a pair may still appear below the top decile in LIANA.

---

### Supplementary Table S15. Developmental-case output summary

| Dataset | Cells | Inferred lineages | Temporal rows | Fork rows | Fork-evaluable rows |
|---------|------:|------------------:|--------------:|----------:|--------------------:|
| pancreas | 3696 | 5 | 272340 | 22695 | 1318 |
| dentategyrus | 2930 | 5 | 247860 | 20655 | 0 |

Fork rows are output rows rather than independent biological events; zero evaluable dentategyrus rows means the fork-contrast gate was not met.

---

### Supplementary Table S16. Context-matched descriptive lag analysis for temporal-only pairs

Expression is log-normalized and restricted to the highest-scoring LineaLinkR source, target and lineage context. Bins use lineage-weighted equal-mass allocation. Positive peak lag and positive cross-correlation shift mean receptor later than ligand. These summaries measure relative timing without a causal test.

| Dataset | Pair | Source→target; lineage | Ligand peak | Receptor peak | Peak lag | Receptor-later xcorr shift | Status |
|---------|------|------------------------|------------:|--------------:|---------:|---------------------------:|--------|
| dentategyrus | Mif–Ackr3 | Granule immature→Endothelial; L2 | 11 | 12 | +1 | −1 | both profiles non-zero |
| dentategyrus | Ccl7–Ccr1 | Microglia→Microglia; L2 | 12 | 12 | 0 | 0 | both profiles non-zero |
| pancreas | Iapp–Calcrl | Beta→Pre-endocrine; L5 | 1 | 11 | +10 | +5 | both profiles non-zero |
| pancreas | Pyy–Npy1r | Alpha→Pre-endocrine; L4 | 1 | 11 | +10 | +3 | both profiles non-zero |
| pancreas | Mdk–Lrp1 | Ductal→Ngn3-low EP; L3 | 12 | 12 | 0 | −1 | both profiles non-zero |
| pancreas | Pdgfa–Pdgfrb | Ngn3-low EP→Ngn3-low EP; L3 | 9 | 12 | +3 | +4 | both profiles non-zero |

### Supplementary Table S17. Global-expression lag analysis and literature context

**A. All-cell raw-count peak-bin lag** (equal-frequency bins on LineaLinkR pseudotime; 500 receptor-profile permutations). This all-cell representation differs from the source-, target- and lineage-matched analysis reported in Supplementary Table S16 and Figure S9.

| Dataset | Pair | Ligand peak bin | Receptor peak bin | Lag (bins) | Lag > 0 | Best xcorr lag | Perm *p* (lag ≥ obs) |
|---------|------|-----------------|-------------------|------------|---------|----------------|----------------------|
| pancreas | Mdk–Sdc2 | 6 | 1 | −5 | no | −3 | 1.000 |
| pancreas | Efna4–Epha1 | 5 | 3 | −2 | no | −1 | 0.812 |
| pancreas | Try5–Pard3 | 5 | 5 | 0 | no | +1 | 0.676 |
| pancreas | Cck–Cckar | 6 | 4 | −2 | no | −3 | 0.742 |
| dentategyrus | Angpt2–Tek | 1 | 1 | 0 | no | 0 | 1.000 |
| dentategyrus | Fgf1–Fgfr2 | 12 | 12 | 0 | no | 0 | 0.084 |

None of the six pairs has a positive receptor-minus-ligand peak difference in this global raw-count representation. Because the representation differs from the source-, target- and lineage-matched analysis, these rows are excluded from the selected temporal-only set.

**B. Literature context.**

| Dataset | Pair | Evidence strength | Evidence scope |
|---------|------|-------------------|----------------|
| pancreas | Mdk–Sdc2 | moderate_family | Family-level developmental evidence; pair-specific endocrine timing was not tested |
| pancreas | Efna4–Epha1 | moderate_family | Family-level pathway evidence; no pancreas timing study was identified |
| pancreas | Try5–Pard3 | weak_or_artifact_risk | No curated developmental evidence was identified |
| pancreas | Cck–Cckar | moderate_pair | Pair documented in gastrointestinal signalling; no positive lag in this analysis |
| dentategyrus | Angpt2–Tek | moderate_pair | Pair documented in the neurogenic niche; temporal order was not tested |
| dentategyrus | Fgf1–Fgfr2 | moderate_family | FGFR-family evidence is richer for Fgf2 than for Fgf1 |

---

### Supplementary Table S18. Patient-level colorectal-cancer validation of prespecified fibroblast–macrophage axes

The six axes were selected from an independent colorectal-cancer study before analysis of GSE132465. Expression proxies use patient as the inferential unit; receptor complexes use the minimum component expression. Ten patients have matched normal and tumour samples.

| Axis | Direction | Mean paired symmetric change | 95% bootstrap CI | Wilcoxon *P* | BH *q* | Tumour greater |
|------|-----------|-----------------------------:|:----------------:|---------------:|-------:|---------------:|
| WNT5A–FZD2 | fibroblast→macrophage | 0.688 | 0.427–0.937 | 0.00391 | 0.00469 | 9/10 |
| CCL3–CCR5 | fibroblast→macrophage | 0.456 | 0.305–0.614 | 0.00391 | 0.00469 | 9/10 |
| RARRES2–CMKLR1 | fibroblast→macrophage | 0.275 | 0.062–0.529 | 0.0371 | 0.0371 | 8/10 |
| IL1A–IL1R1/IL1RAP | macrophage→fibroblast | 0.469 | 0.266–0.660 | 0.00391 | 0.00469 | 9/10 |
| IL1B–IL1R1/IL1RAP | macrophage→fibroblast | 0.396 | 0.265–0.560 | 0.00195 | 0.00469 | 10/10 |
| SPP1–CD44 | macrophage→fibroblast | 0.912 | 0.821–1.000 | 0.00195 | 0.00469 | 10/10 |

Myofibroblasts as a fraction of stromal cells increased in all ten matched patients (median 0.002 to 0.379; mean difference 0.404, 95% CI 0.287–0.534; *P* = 0.00195). The SPP1-positive fraction of myeloid cells was not significantly increased (median 0.436 to 0.496; *P* = 0.275). Across 23 tumour patients, the two fractions were correlated (Spearman ρ = 0.491, *P* = 0.0174).

### Supplementary Table S19. Blind-root versus optional user-root Spearman (6 real GT datasets)

User-root = median cell among earliest ground-truth time labels on the same kNN graph.

| Dataset | Spearman blind | Spearman user-root | Δ (user − blind) |
|---------|----------------|--------------------|------------------|
| lung | −0.068 | −0.083 | −0.015 |
| dentategyrus | 0.350 | −0.029 | −0.378 |
| monocle3_worm_l2 | −0.007 | 0.006 | +0.013 |
| gastrulation_erythroid | 0.238 | 0.849 | +0.611 |
| zebrafish_neural_dev | 0.063 | −0.194 | −0.257 |
| larry_invitro | −0.033 | 0.347 | +0.379 |
| **mean** | **0.091** | **0.149** | **+0.059** |

### Supplementary Table S20. Orthogonal dynamic axes withheld from LineaLinkR

| Dataset | Axis | Mode | n | ρ | τ |
|---------|------|------|---|---|---|
| pancreas | stochastic scVelo pseudotime | blind root | 3696 | 0.936 | 0.779 |
| gastrulation erythroid | embryo stage E7.0–E8.5 | blind root | 9815 | 0.238 | 0.200 |
| gastrulation erythroid | embryo stage E7.0–E8.5 | user-root sensitivity | 9815 | 0.849 | 0.708 |

The scVelo axis is model-derived. Gastrulation shows weak blind-root and strong user-root concordance with withheld stage.

### Supplementary Table S21. Fresh full-workflow wallclock and peak RSS (n=3)

Each run starts from MatrixMarket input and includes Python preprocessing and kNN, followed by LineaLinkR trajectory, static LR and temporal/fork scoring. Peak RSS is sampled over the process tree on Windows.

| Dataset | Cells | Total wallclock, mean±SD (s) | kNN, mean±SD (s) | LLR analysis, mean±SD (s) | Static LR, mean±SD (s) | Rust pseudotime, mean±SD (s) | Peak RSS, mean±SD (MB) |
|---------|-------|------------------------------|-------------------------|---------------------------|-------------------------|-------------------------------|------------------------|
| paul15 | 2730 | 10.636 ± 0.813 | 5.391 ± 0.591 | 2.273 ± 0.051 | 0.313 ± 0.006 | 0.217 ± 0.015 | 291.189 ± 18.830 |
| dentategyrus | 2930 | 15.363 ± 0.382 | 4.657 ± 0.157 | 6.767 ± 0.205 | 1.137 ± 0.021 | 0.210 ± 0.010 | 334.672 ± 1.865 |
| pancreas | 3696 | 18.802 ± 1.058 | 5.868 ± 1.069 | 6.657 ± 0.165 | 1.830 ± 0.092 | 0.330 ± 0.006 | 519.308 ± 1.698 |

### Supplementary Table S22. GSE274389 external time-course and Wnt-recorder analysis

Sampling time and Wnt-recorder state were withheld from LineaLinkR inference. AUC and Spearman summarize separation and ordering. Cluster-adjusted odds ratios use the 5002 cells in Leiden clusters containing at least 20 cells of each recorder state.

| Endpoint | Estimate | 95% CI | n | Interpretation |
|----------|---------:|:------:|--:|----------------|
| Selected-root signed pseudotime versus 96/120 h, Spearman ρ | −0.788 | — | 7942 | Reversed sampling-time direction |
| Selected-root pseudotime, 96 h versus 120 h AUC | 0.027 | — | 7942 | Opposite orientation to sampling time |
| Alternative-root signed pseudotime versus 96/120 h, Spearman ρ | — | — | 7942 | Alternative root not produced under default retention |
| Alternative-root pseudotime, 96 h versus 120 h AUC | — | — | 7942 | Alternative root not produced under default retention |
| Pseudotime, Wnt-recorder AUC at 120 h | 0.453 | — | 5049 | Weak discrimination |
| Wnt-target expression score, recorder AUC | 0.604 | — | 5049 | Modest positive control |
| LineaLinkR Wnt temporal score, recorder AUC | 0.507 | — | 5049 | Weak aggregate discrimination |
| Cluster-adjusted recorder OR per pseudotime IQR | 0.633 | 0.577–0.696 | 5002 | Cell-level association; direction follows the reversed ordering |
| Cluster-adjusted recorder OR per Wnt-target IQR | 1.467 | 1.345–1.601 | 5002 | Cell-level positive control only |
| Cluster-adjusted recorder OR per LineaLinkR Wnt-score IQR | 1.023 | 0.948–1.104 | 5002 | Model-based cell-level association; not replicated |
| Lineage-bin recorder fraction versus Wnt-target score, ρ | 0.419 | — | 5049 | Descriptive positive-control concordance |
| Lineage-bin recorder fraction versus LineaLinkR Wnt score, ρ | −0.010 | — | 5049 | Weak aggregate/bin concordance |

### Supplementary Table S23. CellTag-multi holdout endpoints

Inference used an annotation-free GSE216606 H5AD. Author state, technical batch, day, clone and fate fields were withheld and joined only after fitting. The four-candidate selector retained `detected_genes_high_state`.

| Endpoint | Single returned root | Interpretation |
|----------|---------------------:|----------------|
| Signed Spearman versus day | −0.383 | Root direction is reversed |
| Day-5 AUC, signed | 0.050 | Opposite orientation to day |
| Day-5 AUC, direction-agnostic | 0.950 | Strong separability after orientation is ignored |
| State ARI / NMI | 0.257 / 0.318 | Modest agreement with withheld state labels |
| Fate ARI / NMI | 0.232 / 0.270 | Limited agreement with withheld fate labels |
| Within-clone pair overlap | 0.532 | Greater than global pair overlap 0.255 |
| Normalized clone excess | 0.371 | 95% CI 0.354–0.390; 5066 clones |

---

### Supplementary Table S24. Cross-dataset skeletal-muscle regeneration and FAP depletion

GSE143437 recorded sampling day anchors the C3–C3AR1 expression analysis because blind LineaLinkR ordering does not track day. GSE291091 estimates are sample-level pseudobulk contrasts with two biological replicates per group. Broad phagocyte enrichment is composition-sensitive.

| Layer | Endpoint | Estimate | P or FDR | Evidence note |
|-------|----------|---------:|---------:|---------------|
| GSE143437 blind ordering | Pseudotime versus recorded day, Spearman ρ | −0.060 | — | Recorded day is used as the biological time axis |
| GSE143437 samples | Mean C3–C3AR1 score at days 0/2/5/7 | 0.844/0.661/0.709/0.625 | — | 3/2/3/2 samples per day |
| GSE291091 mesenchymal/FAP-like compartment | C3, depleted minus control at 3 dpi, log2FC | −0.640 | 0.484 | No FDR-supported C3 reduction; *n* = 2+2 |
| GSE291091 mononuclear phagocytes | C3, depleted minus control at 3 dpi, log2FC | −0.531 | 0.0738 | Nominal reduction, not FDR-supported; *n* = 2+2 |
| GSE291091 mononuclear phagocytes | Il6, depleted minus control at 3 dpi, log2FC | −2.360 | 0.00952 | FDR-supported reduction; *n* = 2+2 |
| GSE291091 mononuclear phagocytes | IL6/JAK/STAT3 Hallmark NES | −1.880 | 3.82 × 10−4 | Negative rank enrichment after depletion |
| GSE291091 mononuclear phagocytes | Inflammatory-response Hallmark NES | −2.415 | 1.09 × 10−11 | Negative rank enrichment after depletion |
| GSE291091 mononuclear phagocytes | IFN-γ-response Hallmark NES | −3.155 | 1.96 × 10−42 | Negative rank enrichment after depletion |
| GSE291091 mononuclear phagocytes | TNF/NF-κB Hallmark NES | −2.042 | 3.03 × 10−7 | Negative rank enrichment after depletion |
| GSE291091 differential-expression summary | FAP / phagocyte genes at FDR < 0.05 | 317 / 1360 | — | Whole-population depletion response; not C3-specific |
| GSE291091 control mesenchymal/FAP-like compartment | C3, 3 dpi minus uninjured, log2FC | −0.435 | 0.442 | No mouse 3-dpi induction; differs from the elderly-human spatial result |
| GSE291091 phagocyte composition | Mixed mono/macrophage fraction, control versus depleted | 0.49–0.50 vs 0.61–0.74 | — | Captured subtype composition shifts after depletion |
| GSE291091 phagocyte composition | Pro-inflammatory fraction, control versus depleted | 0.11–0.16 vs 0.06–0.07 | — | Broad-compartment pathway NES is composition-sensitive |
| GSE291091 phagocyte composition | Anti-inflammatory fraction, control versus depleted | 0.24–0.27 vs 0.10–0.20 | — | Broad-compartment pathway NES is composition-sensitive |

### Supplementary Table S25. GSE210152 patient×Leiden purity

| Leiden | n cells | Max patient fraction | Dominant patient | Norm. patient entropy |
|--------|---------|----------------------|------------------|-----------------------|
| 0 | 1300 | 0.248 | P7 | 0.954 |
| 1 | 8005 | 0.248 | P6 | 0.974 |
| 2 | 2771 | 0.257 | P7 | 0.976 |
| 3 | 518 | 0.498 | P3 | 0.689 |
| 4 | 21340 | 0.248 | P7 | 0.977 |
| 5 | 260 | 0.415 | P3 | 0.848 |
| 6 | 262 | 0.981 | P1 | 0.063 |

### Supplementary Table S26. GSE210152 exact-barcode plaque-atlas mapping

Fractions are calculated among mapped cells; clusters with low mapping coverage are not assigned from the mapped-cell majority alone.

| ID | n local | n mapped | Map frac. | Major | Major frac. |
|----|---------|----------|-----------|-------|-------------|
| 0 | 1300 | 1104 | 0.849 | B cell | 0.592 |
| 1 | 8005 | 6324 | 0.790 | Macrophage | 0.488 |
| 2 | 2771 | 18 | 0.006 | T cell | 0.778 |
| 3 | 518 | 160 | 0.309 | Plasma cell | 0.963 |
| 4 | 21340 | 17382 | 0.815 | T cell | 0.813 |
| 5 | 260 | 205 | 0.788 | Mast cell | 0.985 |
| 6 | 262 | 0 | 0.000 | unmapped | — |

### Supplementary Table S27. GSE210152 processed-atlas membership and basic QC

The basic QC rule is detected genes >200, counts >500 and <40 000, and mitochondrial fraction <10%, matching the atlas paper's stated cell-level thresholds. Atlas processing additionally included ambient-RNA correction and doublet detection.

| Group | n | QC pass n | QC pass | Counts med. | Genes med. | mt% med. |
|-------|---|-----------|---------|-------------|------------|----------|
| Atlas-mapped | 25193 | 25181 | 0.9995 | 2841 | 978 | 2.624 |
| Not in atlas | 9263 | 3461 | 0.3736 | 551 | 267 | 6.849 |

Association summary: Leiden×membership Cramér's V = 0.521; patient×membership V = 0.062; basic-QC pass versus membership odds ratio = 3517.8. P values are reported in the source CSV and are not interpreted as biological effect sizes.

### Supplementary Table S28. GSE210152 atlas-coverage temporal-LR sensitivity

Temporal scores are recomputed on retained cells with the pooled pseudotime and lineage weights fixed, isolating LR-score sensitivity from trajectory reconstruction.

| Setting | n cells | n traj. | n common | Rank ρ | Top-30 Jaccard |
|---------|--------:|--------:|---------:|-------:|---------------:|
| Full | 34456 | 3698 | 3698 | 1.000 | 1.000 |
| Exclude Leiden 2/3 | 31167 | 3158 | 3046 | 0.901 | 0.500 |
| Exclude Leiden 3 | 33938 | 3517 | 3449 | 0.958 | 0.714 |
| Exact-atlas cells only | 25193 | 3465 | 2904 | 0.612 | 0.071 |

### Supplementary Table S29. GSE210152 patient and leave-one-patient-out rank stability

Patient effects are temporal-LR rescoring results for the pooled top-100 trajectory pool; the table reports the top 30. Rank checks use temporal amplitude, without a temporal-permutation or causal-significance test.

| Patient | n eval. | Patient ρ | Patient J15 | n LOO | LOO ρ | LOO J15 |
|---------|--------:|----------:|------------:|------:|------:|--------:|
| P1 | 26 | −0.049 | 0.364 | 26 | 0.989 | 1.000 |
| P2 | 27 | −0.153 | 0.250 | 27 | 0.998 | 1.000 |
| P3 | 28 | −0.121 | 0.364 | 28 | 0.878 | 1.000 |
| P5 | 27 | 0.070 | 0.429 | 27 | 0.814 | 0.765 |
| P6 | 29 | 0.078 | 0.364 | 29 | 0.936 | 0.875 |
| P7 | 24 | −0.348 | 0.200 | 24 | 0.999 | 0.765 |

Across patients, median patient-versus-pooled ρ is −0.085 and median top-15 Jaccard is 0.364. Across LOO analyses, median ρ is 0.963 and median Jaccard is 0.938. Median fixed-effect heterogeneity *I*² is 0.863.

### Supplementary Table S30. Primary-source feature comparison of temporal CCC methods

Feature definitions follow the primary method papers [35, 36, 38].

| Method | Trajectory handling | Branch context | Repeated times | Native LR output | Lag/direction model | Bifurcation-relative role |
|--------|---------------------|----------------|----------------|------------------|---------------------|---------------------------|
| LineaLinkR | internal blind/user-root TI; lineage weights | source×target×lineage×bin | no | saturating bin-resolved LR activity | descriptive lag summary | five operational labels plus an auxiliary transition score; pancreas prior agreement in Tables S12 and S13 |
| TICCI | internal Louvain/PAGA, scEntropy direction and fitted time | trajectory topology, not LR profiles | no | CCI-weighted neighbourhood graph for TI | no LR lag test | no |
| TraSig | external trajectory; sliding windows on edges | yes | no; uses supplied pseudotime | LR profile dot products with permutation support | synchronous profile matching | no |
| TimeTalk | reconstructed pseudotime plus temporal TF network | cell-type/pseudotime | yes for intended workflow | co-varying LR scores linked to TF activity | LR↔TF Granger direction and correlation filters | no; its forward/backward/feedback labels are LR–TF directions |

## Supplementary Notes

### Supplementary Note 1. Benchmark eligibility

Algorithm-dataset pairs are eligible when the required native inputs are available and the method returns a finite, cell-aligned endpoint. Ordering comparisons use common evaluated datasets; endpoint-only outputs are reported separately. Terminal-branch comparisons use the four two-terminal simulations listed in Table S1 and cells with reference pseudotime ≥0.4.

### Supplementary Note 2. Ablation design

The pancreas module ablation removes LR scoring, temporal binning and fork detection independently from the same upstream trajectory. Graph-*k* variants on `sim_2br_500c` recompute the trajectory.

### Supplementary Note 3. Main-figure scope

Figure 1 presents the workflow; Figures 2–4 cover trajectory geometry and benchmark comparisons; Figures 5–6 show developmental applications; Figures 7–8 present external temporal, clonal, perturbational and patient-level evidence.

### Supplementary Note 4. Software versions and environment

| Component | Version |
|-----------|---------|
| R | 4.5.3 |
| LineaLinkR | 0.6.4 |
| CellPhoneDB | 5.0 |
| CellChat | 2.0.0 |
| LIANA | 1.0.0 |
| slingshot | 2.18.0 |
| monocle3 | 1.4.0 |
| Random seed | 12345 |

The analysis environment and dependency versions are recorded with the software release.

---

### Supplementary Note 5. MM-LR Hill scale and sensitivity

MM-LR uses a data-driven Hill scale *S*. Sixteen datasets with non-zero LR scores are evaluable in the multiplier sensitivity sweep. At *k*=5, median top-30 Jaccard versus *k*=1 is 0.905 (minimum 0.765); the remaining datasets have zero LR scores and are not evaluable for this sweep.

---

### Supplementary Note 6. GSE210152 partial plaque-atlas annotation

The deposited GEO H5AD lacks a cell-type field. Exact barcode joining to the processed 2025 plaque atlas maps 25 193 of 34 456 cells (73.1%) without imputation. Basic QC passes for 25 181/25 193 mapped cells (99.95%) versus 3461/9263 non-members (37.36%). Leiden 0, 4 and 5 have B-cell, T-cell and mast-cell majorities; Leiden 1 is mixed macrophage enriched, Leiden 3 has partial plasma-cell support, Leiden 2 has only 18 mapped cells and Leiden 6 is unmapped. Patient dominance is strongest in Leiden 6 (P1 fraction 0.981).

Coverage perturbations recompute temporal LR scores with pooled pseudotime and lineage weights fixed. Excluding Leiden 3 gives rank ρ = 0.958 and top-30 Jaccard 0.714; excluding Leiden 2/3 gives 0.901 and 0.500. Exact-atlas-only restriction gives 0.612 and 0.071. Individual-patient versus pooled rank Spearman has median −0.085 and top-15 Jaccard 0.364, whereas leave-one-patient-out rankings have median ρ = 0.963 and top-15 Jaccard 0.938; median *I*² = 0.863. These are pooled Leiden-level results.

---

## Supplementary Figure legends

**Supplementary Figure S1.** Pancreas soft-lineage assignment: UMAP coloured by categorical dominant lineage and by maximum lineage weight. Maximum lineage weight quantifies lineage concentration.

Alt text: Two pancreas UMAP panels. The first uses a discrete five-lineage legend; the second uses a continuous colour scale for the maximum soft-lineage weight.

**Supplementary Figure S2.** Soft-lineage uncertainty in pancreas. (A) Lineage-weight glyphs on a fixed cell subset show local mixtures on the UMAP geometry. (B) Hexagonal cell density and binned medians relate maximum lineage weight to LineaLinkR pseudotime.

Alt text: UMAP overlay of pie glyphs encoding soft lineage probabilities per cell.

**Supplementary Figure S3.** Temporal bin-number sensitivity on pancreas: Jaccard of top-30 pairs versus the default 12-bin run, and number of lineages with a detected fork bin, for bin counts {8, 12, 16, 20}.

Alt text: Two-panel sensitivity plot showing top-pair Jaccard stability and bifurcation-lineage counts across four temporal bin settings.

**Supplementary Figure S4.** Full-panel Hill-S multiplier sensitivity (*k* ∈ {0.5, 1, 2, 5}). Panel A shows top-30 Jaccard distributions for 16 datasets with non-zero LR scores; the dashed 0.70 line is descriptive. Panel B shows the dataset-resolved matrix.

Alt text: Box-and-point distributions and a 16-dataset heatmap show strong but non-uniform top-pair stability as the Hill scale changes; the lowest Jaccard is 0.765 at five-fold scale.

**Supplementary Figure S5.** Pancreas top temporal LR lineage-bootstrap CIs.

Alt text: Forest plot of per-pair mean scores with 95% confidence intervals for pancreas top temporal programs.

**Supplementary Figure S6.** Spearman correlation matrix for the 19 complete-core datasets and 12 TI algorithms. The upper summary shows each method's mean across the same core; grey cells mark missing outputs.

Alt text: Heatmap of pseudotime Spearman correlations across 19 datasets and 12 algorithms, with grey missing-output cells and a method-mean summary.

**Supplementary Figure S7.** Rank distributions for the 19-dataset, 12-method complete core. Dataset-level ranks are shown with method summaries, top-three frequency and win fraction; lower rank is better. The Friedman test gives χ² = 74.13 with 11 df and *p* = 1.99 × 10⁻¹¹.

Alt text: Dataset-level rank distributions, top-three frequency and win fraction for 12 trajectory algorithms on the same 19 datasets.

**Supplementary Figure S8.** Module-output dependence across 17 evaluable datasets. The matrix reports log10(1 + observed count) for static LR pairs, temporal rows and fork points in the full analysis and after removing LR scoring, temporal binning or fork detection. This is an implementation-level output check rather than an accuracy comparison.

Alt text: A 17-dataset matrix shows that removing LR scoring eliminates all communication outputs, removing temporal binning preserves static LR pairs only, and removing fork detection preserves static and temporal outputs while eliminating fork points.

**Supplementary Figure S9.** Context-matched descriptive lag analysis for six selected temporal-only pairs. Panel A shows source/target-cell-type expression profiles across approximate lineage-weighted pseudotime bins. Panel B compares peak-bin differences with receptor-later cross-correlation shifts; positive values indicate a later receptor profile.

Alt text: Paired ligand and receptor heatmap rows show context-matched profiles for six selected temporal-only pairs; a dot plot compares peak-bin and cross-correlation lag summaries.

**Supplementary Figure S10.** Fine-subtype CRC expression proxies across 23 tumour patients. Values are within-axis z scores for SPP1-positive myeloid-to-myofibroblast or reverse-direction expression proxies; grey cells have fewer than five cells in either subtype.

Alt text: Heatmap of six prespecified LR expression proxies across tumour patients, with missing cells explicitly shown when either author-defined subtype is sparse.

**Supplementary Figure S11.** Orthogonal dynamic-axis validation with axes withheld from LineaLinkR inference. Pancreas stochastic scVelo pseudotime is shown on UMAP and against blind LineaLinkR pseudotime; gastrulation boxplots compare blind-root and user-root pseudotime against embryo sampling stage.

Alt text: Four panels show strong pancreas concordance with RNA velocity, weak positive blind-root gastrulation-stage concordance (Spearman 0.24), and stronger user-root concordance (0.85).

**Supplementary Figure S12.** Three independent workflow runs per dataset on paul15, dentategyrus and pancreas. Panels show total wallclock mean±SD, mean stage composition and process-tree peak RSS mean±SD. Each run includes preprocessing, kNN, R trajectory inference and static/temporal/fork LR scoring.

Alt text: Three horizontal panels compare repeated total wallclock, stage composition and peak memory for three datasets. Total wallclock is approximately 11–32 s and peak RSS approximately 281–563 MB.

**Supplementary Figure S13.** GSE291091 captured macrophage-subtype composition at 3 days post injury. Circles are sample-level fractions within the broad mononuclear-phagocyte compartment for two control and two FAP-depleted samples; diamonds are group means and horizontal segments connect the means. Subtypes are auxiliary GSE143437 label-transfer classes at confidence ≥0.50; lower-confidence cells are retained as a separate category.

Alt text: Dot plot of five transferred macrophage categories in four samples. Depleted samples have a larger mixed mono/macrophage fraction and generally smaller pro- and anti-inflammatory fractions than controls, showing that broad-compartment pathway contrasts are composition-sensitive.

**Supplementary Figure S14.** GSE210152 atherosclerosis temporal programs ranked by mean temporal score. Point area gives the number of source–target contexts and colour denotes communication family. Outputs remain Leiden-resolved; exact-barcode plaque-atlas context is reported separately in Supplementary Figure S16.

Alt text: Lollipop plot ranking 15 temporal ligand–receptor pairs, with source–target context count encoded by point area and communication family by colour.

**Supplementary Figure S15.** GSE210152 temporal profiles for the same 15 LR pairs. Values are scaled within each pair across 12 pseudotime bins to show temporal shape independently of absolute score.

Alt text: Heatmap showing each of 15 temporal ligand–receptor pairs across 12 pseudotime bins after within-pair scaling.

**Supplementary Figure S16.** GSE210152 reference annotation by exact barcode join to the 2025 plaque atlas. Panels show Leiden-by-level-1 composition, per-cluster mapping coverage and major-label fraction among mapped cells. Leiden 2 and 6 have insufficient mapping; Leiden 3 is partially mapped and plasma-cell enriched among mapped cells.

Alt text: A heatmap and two horizontal-bar panels show partial plaque-atlas annotation of seven local Leiden clusters, including high-confidence B-cell, T-cell and mast-cell majorities and two insufficiently mapped clusters.

**Supplementary Figure S17.** GSE210152 processed-atlas membership and QC. ECDF panels compare detected genes and mitochondrial fraction for exact-atlas members versus non-members; the heatmap resolves mapping coverage by patient and Leiden; the final panel compares mapping and basic-QC pass fractions by Leiden.

Alt text: Four panels show that non-members of the processed plaque-atlas subset have markedly lower detected genes and higher mitochondrial fractions, while mapping gaps concentrate in Leiden clusters 2 and 6 with partial coverage of cluster 3 rather than in one patient.

**Supplementary Figure S18.** GSE210152 atlas-coverage and patient-rescoring sensitivity. Panel A compares temporal-LR rank and top-set stability under three cell restrictions; panel B shows heterogeneous patient-specific rankings; panel C shows fixed-effect leave-one-patient-out meta rankings; panel D displays random-effects temporal amplitudes for leading programs.

Alt text: Four panels show that exact-atlas-only restriction changes the pooled LR ranking, individual patients have heterogeneous rankings, exclusion of any one patient leaves the fixed-effect meta ranking largely stable, and random-effects estimates for leading pooled amplitudes have wide uncertainty.
