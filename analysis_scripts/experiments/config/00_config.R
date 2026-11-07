# =============================================================================
# LineaLinkR 0.6.5 analysis configuration
# =============================================================================

CONFIG_FILE <- local({
  frame_paths <- vapply(sys.frames(), function(fr) {
    if (exists("ofile", envir = fr, inherits = FALSE)) {
      get("ofile", envir = fr)
    } else {
      NA_character_
    }
  }, character(1))
  frame_paths <- frame_paths[!is.na(frame_paths)]
  normalizePath(frame_paths[length(frame_paths)], mustWork = TRUE)
})
WORK_ROOT <- normalizePath(Sys.getenv(
  "LLR_WORK_ROOT", dirname(dirname(dirname(CONFIG_FILE)))
), winslash = "/", mustWork = TRUE)
INPUT_DATA_ROOT <- normalizePath(
  file.path(WORK_ROOT, "inputs", "benchmark"),
  winslash = "/", mustWork = FALSE
)
READONLY_REVISION_ROOT <- normalizePath(
  file.path(WORK_ROOT, "inputs", "revision"), winslash = "/", mustWork = FALSE
)
CANONICAL_RESULT_ROOT <- normalizePath(
  file.path(READONLY_REVISION_ROOT, "results"),
  winslash = "/", mustWork = FALSE
)
SOURCE_ASSETS_ROOT <- WORK_ROOT
ROOT <- WORK_ROOT
ASSETS <- WORK_ROOT
CONFIG <- list(

  dir = list(
    root         = WORK_ROOT,
    assets       = WORK_ROOT,
    source_assets = SOURCE_ASSETS_ROOT,
    pipeline     = file.path(WORK_ROOT, "experiments"),
    manuscript   = file.path(WORK_ROOT, "manuscript"),
    exploratory  = file.path(WORK_ROOT, "experiments", "exploratory"),
    qa           = file.path(WORK_ROOT, "experiments", "qa"),
    data         = INPUT_DATA_ROOT,
    curated      = file.path(INPUT_DATA_ROOT, "curated"),
    literature   = file.path(INPUT_DATA_ROOT, "literature_benchmarks"),
    canonical_input = file.path(INPUT_DATA_ROOT, "canonical"),
    canonical    = file.path(WORK_ROOT, "inputs", "project_canonical"),
    algorithms   = file.path(WORK_ROOT, "tools", "competitors", "optimized"),
    benchmarks   = file.path(WORK_ROOT, "benchmarks"),
    docs         = file.path(WORK_ROOT, "docs"),
    output       = file.path(WORK_ROOT, "results"),
    competitor_baseline = CANONICAL_RESULT_ROOT,
    paper_tables = file.path(WORK_ROOT, "results", "paper_tables"),
    table1_main  = file.path(WORK_ROOT, "results", "paper_tables", "table1_main"),
    reference    = file.path(WORK_ROOT, "results", "paper_tables", "reference"),
    evaluation   = file.path(WORK_ROOT, "results", "paper_tables", "table1_main"),
    evaluation_audit = file.path(WORK_ROOT, "results", "paper_tables", "audit"),
    supplement   = file.path(WORK_ROOT, "results", "paper_tables", "supplement"),
    supplement_tables = file.path(WORK_ROOT, "results", "paper_tables", "supplement", "tables"),
    algorithm_out= file.path(SOURCE_ASSETS_ROOT, "results", "by_algorithm"),
    algorithm_misc = file.path(WORK_ROOT, "results", "misc"),
    cache        = file.path(WORK_ROOT, "cache"),
    linealinkr   = file.path(WORK_ROOT, "cache", "feasible_clock"),
    r_export     = file.path(SOURCE_ASSETS_ROOT, "cache", "r_export"),
    log          = Sys.getenv("LLR_LOG_DIR", file.path(WORK_ROOT, "logs", "analysis")),
    archive      = file.path(WORK_ROOT, "archive"),
    benchmark_exclusions = file.path(WORK_ROOT, "experiments", "config", "benchmark_exclusions.tsv")
  ),

  lr_resources = list(
    celegans = file.path(WORK_ROOT, "experiments", "resources", "celegans", "cell2cell_2022_celegans_lr.csv")
  ),

  datasets_curated = list(
    lung                     = "lung_benchmark.h5ad",
    paul15                   = "paul15_benchmark.h5ad",
    pancreas                 = "pancreas_benchmark.h5ad",
    dentategyrus             = "dentategyrus_benchmark.h5ad",
    cellagentchat_scc        = "cellagentchat_scc_processed_benchmark.h5ad",
    monocle3_worm_l2         = "monocle3_worm_l2_benchmark.h5ad"
  ),

  datasets_literature = list(
    gastrulation_erythroid = "gastrulation_erythroid/gastrulation_erythroid.h5ad",
    zebrafish_urd          = "zebrafish_urd/zebrafish_urd.h5ad",
    zebrafish_neural_dev   = "zebrafish_neural_dev/zebrafish_neural_dev.h5ad",
    human_cd34             = "human_cd34_bone_marrow/human_cd34_bone_marrow.h5ad",
    larry_invitro          = "larry_invitro/larry_invitro.h5ad",
    simulated_branch       = "simulated_branch/simulated_branch.h5ad"
  ),

  datasets_project_canonical = list(
    human_fetal_liver_organoid = "human_fetal_liver_organoid.h5ad",
    human_fetal_pancreas = "human_fetal_pancreas.h5ad",
    zebrafish_tracerseq = "zebrafish_tracerseq.h5ad"
  ),

  # Raw GEO bundles (pre-conversion); h5ad lives under geo_candidates / geo_narrative
  datasets_geo_raw = list(
    GSE125527_IBD   = "geo_curated/GSE125527",
    GSE132465_CRC   = "geo_curated/GSE132465",
    GSE135893_IPF   = "geo_curated/GSE135893",
    GSE171524_COVID = "geo_curated/GSE171524"
  ),

  datasets_simulated = list(
    sim_2br_500c             = "sim_2br_500c_canonical.h5ad",
    sim_2br_2000c            = "sim_2br_2000c_canonical.h5ad",
    sim_2br_5000c            = "sim_2br_5000c_canonical.h5ad",
    sim_3br_500c             = "sim_3br_500c_canonical.h5ad",
    sim_3br_2000c            = "sim_3br_2000c_canonical.h5ad",
    sim_3br_5000c            = "sim_3br_5000c_canonical.h5ad",
    sim_4br_2000c            = "sim_4br_2000c_canonical.h5ad",
    sim_4br_500c             = "sim_4br_500c_canonical.h5ad",
    sim_4br_5000c            = "sim_4br_5000c_canonical.h5ad",
    sim_3br_2000c_batch      = "sim_3br_2000c_batch_canonical.h5ad",
    sim_3br_2000c_drop30     = "sim_3br_2000c_drop30_canonical.h5ad",
    sim_3br_2000c_drop50     = "sim_3br_2000c_drop50_canonical.h5ad",
    sim_3br_2000c_drop70     = "sim_3br_2000c_drop70_canonical.h5ad",
    sim_3br_2000c_drop90     = "sim_3br_2000c_drop90_canonical.h5ad",
    simulated_complex_3branch = "simulated_complex_3branch_canonical.h5ad"
  ),

  # Revision 0.6.5 same-source counts->LR->fork simulations.
  # Inputs remain in READONLY_REVISION_ROOT/sim/.
  datasets_revision_sim = list(
    lrsim_shared_trend     = "lrsim_shared_trend_canonical.h5ad",
    lrsim_focal_branch     = "lrsim_focal_branch_canonical.h5ad",
    lrsim_donor_confounded = "lrsim_donor_confounded_canonical.h5ad",
    lrsim_dropout          = "lrsim_dropout_canonical.h5ad",
    lrsim_sampling_thinning = "lrsim_sampling_thinning_canonical.h5ad",
    lrsim_nested_fork      = "lrsim_nested_fork_canonical.h5ad",
    lrsim_transient_focal  = "lrsim_transient_focal_canonical.h5ad",
    lrsim_donor_ligand_only = "lrsim_donor_ligand_only_canonical.h5ad",
    lrsim_donor_receptor_only = "lrsim_donor_receptor_only_canonical.h5ad",
    lrsim_donor_nested_fork = "lrsim_donor_nested_fork_canonical.h5ad",
    lrsim_donor_transient_focal = "lrsim_donor_transient_focal_canonical.h5ad"
  ),

  datasets_geo_candidate = list(
    GSE243917_RA     = "geo_candidates/GSE243917/GSE243917_genesXcells_ALL_raw.h5ad",
    GSE210152_athero = "geo_candidates/GSE210152/GSE210152_raw.h5ad",
    GSE125527_IBD    = "geo_candidates/GSE125527/GSE125527_IBD_raw.h5ad",
    GSE132465_CRC    = "geo_candidates/GSE132465/GSE132465_CRC_raw.h5ad"
  ),

  datasets_revision_author_cohort = list(
    GSE135893_IPF = "../author_cohort/GSE135893_IPF_author_cohort.h5ad",
    # Historical key retained; the unchanged matrix is Cao 2017 L2, not embryo.
    celegans_embryo = "../author_cohort/cao2017_l2_author_annotations.h5ad"
  ),

  # Biological cases and processed versions of existing benchmark objects.
  datasets_geo_narrative = list(
    cellrank_stemnet = "curated/cellrank_stemnet_benchmark.h5ad",
    GSE171524_COVID = "geo_narrative/GSE171524/GSE171524_COVID_raw.h5ad",
    GSE274389_signal_recorder = "canonical/GSE274389_signal_recorder_canonical.h5ad",
    GSE143437_muscle_regeneration = "canonical/GSE143437_muscle_regeneration_canonical.h5ad"
  ),

  # Independent longitudinal source and receiver geometry views.
  # The receiver views below are not additional independent datasets.
  datasets_revision_narrative = list(
    GSE216451_FGFR2_P7 = "sim/GSE216451_fgfr2_perturbation/GSE216451_P7_canonical.h5ad",
    GSE216451_FGFR2_P21 = "sim/GSE216451_fgfr2_perturbation/GSE216451_P21_canonical.h5ad",
    GSE216451_FGFR2_P7_secretory = "sim/GSE216451_fgfr2_perturbation/GSE216451_P7_secretory.h5ad",
    GSE216451_FGFR2_P21_secretory = "sim/GSE216451_fgfr2_perturbation/GSE216451_P21_secretory.h5ad",
    GSE216451_FGFR2_P21_epithelial = "sim/GSE216451_fgfr2_perturbation/GSE216451_P21_epithelial.h5ad",
    GSE270807_HLBO = "cache/biology_inputs/GSE270807_HLBO.h5ad",
    GSE207889_ULA = "cache/biology_inputs/GSE207889_ULA.h5ad",
    GSE302266_kidney = "cache/biology_inputs/GSE302266_kidney.h5ad",
    GSE302266_kidney_author_cells = "cache/biology_inputs/GSE302266_kidney_author_cells.h5ad",
    GSE300309_DMSO = "cache/biology_inputs/GSE300309_DMSO.h5ad",
    GSE300309_DAPT = "cache/biology_inputs/GSE300309_DAPT.h5ad",
    GSE265803_mdx5cv_quadriceps = "cache/biology_inputs/GSE265803_mdx5cv_quadriceps.h5ad",
    GSE265803_mdx5cv_diaphragm = "cache/biology_inputs/GSE265803_mdx5cv_diaphragm.h5ad",
    GSE265803_Ccr2KO_quadriceps = "cache/biology_inputs/GSE265803_Ccr2KO_quadriceps.h5ad",
    GSE265803_Ccr2KO_diaphragm = "cache/biology_inputs/GSE265803_Ccr2KO_diaphragm.h5ad",
    GSE118184_Takasato_timecourse = "cache/biology_inputs/GSE118184_Takasato_timecourse.h5ad",
    GSE181595_kidney = "cache/biology_inputs/GSE181595_kidney.h5ad",
    GSE124472_kidney = "cache/biology_inputs/GSE124472_kidney.h5ad",
    GSE184928_kidney = "cache/biology_inputs/GSE184928_kidney.h5ad",
    GSE263556_kidney = "cache/biology_inputs/GSE263556_kidney.h5ad",
    GSE164564_kidney = "cache/biology_inputs/GSE164564_kidney.h5ad",
    GSE143437_phagocyte_receiver = "cache/biology_inputs/GSE143437_phagocyte_receiver.h5ad",
    GSE143437_myogenic_receiver = "cache/biology_inputs/GSE143437_myogenic_receiver.h5ad",
    GSE135893_epithelial_receiver = "cache/biology_inputs/GSE135893_epithelial_receiver.h5ad",
    GSE135893_alveolar_receiver = "cache/biology_inputs/GSE135893_alveolar_receiver.h5ad",
    GSE271604_gastruloid_42h = "cache/biology_inputs/GSE271604/GSE271604_42h_rep2_canonical.h5ad",
    GSE136689_foregut_DE_SM = "cache/biology_inputs/GSE136689/GSE136689_foregut_DE_SM.h5ad",
    GSE115956_intestinal_organoid = "cache/biology_inputs/GSE115956/GSE115956_intestinal_organoid.h5ad",
    GSE239582_FT_WNT7A_FZD5 = "cache/biology_inputs/GSE239582/GSE239582_FT_WNT7A_FZD5.h5ad"
  ),

  n_threads = as.integer(Sys.getenv("LLR_THREADS", as.character(parallel::detectCores()))),
  seed      = 12345L,

  exe = list(
    rscript = file.path(R.home("bin"),ifelse(.Platform$OS.type=="windows","Rscript.exe","Rscript")),
    python  = Sys.getenv("LLR_PYTHON",unset="python"),
    mingw_bin = "E:/software/QTsoftware/Tools/mingw1310_64/bin"
  ),

  llr = list(
    benchmark_protocol  = "linealinkr065_common_grid",
    # Bump when LineaLinkR R/C++ source changes; invalidates cached linealinkr_complete.rds
    result_revision     = "linealinkr065_cellchat_complex_unique_pairs",
    package_version     = "0.6.5",
    root_strategy       = "attractor_energy",
    root_potency_mode   = "detected_genes_ordering",
    clock_depth_strata_columns = c("donor", "sample_id", "orig.ident",
      "sample_name", "sample", "replicate", "batch"),
    root_endpoint_alignment_min_gain = 0.20,
    root_endpoint_require_terminal_leaf = TRUE,
    root_ordering_min_gain = 0.10,
    root_representation_alignment_min_gain = 0.00,
    root_ordering_fallback_threshold = 0.45,
    root_fallback_k_max = 25L,
    root_fallback_k_min = 8L,
    root_fallback_cells_per_cluster = 50L,
    root_fallback_nstart = 10L,
    root_fallback_kmeans_algorithm = "Hartigan-Wong",
    root_fallback_selection_score_mode = "confidence_regularized",
    root_fallback_weights = c(0.30, 0.25, 0.20, 0.15, 0.10),
    lineage_embedding   = "pca",
    lineage_topology_k  = 15L,
    lineage_topology_nstart = 10L,
    root_top_fraction   = 0.05,
    root_n_attractors   = 0L,
    root_max_iter       = 300L,
    root_tol            = 1e-5,
    root_weights        = c(0.25, 0.20, 0.30, 0.15, 0.10),
    root_direction_n_pc = 3L,
    root_direction_min_potency_cor = 0.10,
    root_direction_endpoint_fraction = 0.05,
    root_direction_min_potency_delta = 0.10,
    root_direction_min_blind_delta = 0.05,
    root_ensemble_k_grid = c(30L, 50L, 70L),
    root_ensemble_pc_grid = c(10L, 20L, 30L),
    root_ensemble_resamples = 20L,
    root_ensemble_min_resamples = 3L,
    root_ensemble_adaptive_resampling = TRUE,
    root_ensemble_subsample_fraction = 0.80,
    root_ensemble_direction_margin = 0.10,
    n_bins = 12L,
    min_n_eff = 5.0,
    fork_window_bins = 2L,
    fork_activity_floor = 1e-4,
    fork_sample_boot = 999L,
    start_cell_method   = "auto",
    prefer_root_label   = FALSE,
    rank_normalize      = FALSE,
    force_kmeans_refine = TRUE,
    root_partition_max_clusters = 8L,
    crossing_min_retained_variance = 0.75,
    knn_neighbors       = 50L,
    lr_scale_factor     = 10000,
    # ConnectKnnComponents supplies the minimum bridges needed for coverage.
    # Additional cluster-MST edges are opt-in because they can shortcut time.
    augment_mst_bridges = FALSE,
    connect_components  = TRUE
  ),

  lr_db = list(
    human = list(
      cellchat_db = system.file("extdata", "CellChatDB.human.rds",
                        package = "CellChat", mustWork = FALSE),
      omni_path   = file.path(SOURCE_ASSETS_ROOT, "02_vendor_algorithms", "LIANA",
                        "src", "liana", "resource", "omni_resource.csv")
    ),
    mouse = list(
      cellchat_db = system.file("extdata", "CellChatDB.mouse.rds",
                        package = "CellChat", mustWork = FALSE)
    )
  ),

  algorithm_groups = list(
    ti_trajectory = c("slingshot", "DPT", "SCORPIUS", "MultistageOT",
                      "scFates", "pCreode", "monocle3", "TSCAN", "Palantir",
                      "StaVia", "scTour"),
    ti_nonbenchmark = c("FateID", "Tempora", "cellrank2", "SLICER", "STREAM",
                        "Wishbone", "PHLOWER"),
    ccc_static    = c("SingleCellSignalR", "CellChat", "NATMI", "CellPhoneDB",
                      "LIANA", "nichenetr", "CytoTalk"),
    ccc_temporal  = c("TICCI", "TraSig", "TimeTalk"),
    ti_de_other   = c("PAGA", "tradeSeq")
  ),

  # Full-panel execution order: light (O(n)–O(n log n) graph) → heavy (O(n²) OT / neural).
  # Used by 01_run_other_algos.R; --algos= keeps this order among selected algos.
  algorithm_run_order = c(
    "PAGA", "slingshot", "DPT", "monocle3", "SCORPIUS", "TSCAN",
    "Palantir", "scFates", "NATMI", "LIANA", "nichenetr",
    "SingleCellSignalR", "CellChat", "pCreode",
    "StaVia", "scTour", "TICCI", "TimeTalk", "TraSig",
    "CytoTalk", "CellPhoneDB", "MultistageOT"
  ),

  ticci = list(
    root    = file.path(ROOT, "tools", "competitors", "optimized", "TICCI"),
    py      = file.path(ROOT, "tools", "competitors", "optimized", "TICCI", "Python_tools"),
    r_tools = file.path(ROOT, "tools", "competitors", "optimized", "TICCI", "R_tools")
  ),

  # Only registered experimental sampling times enter calendar inference.
  dataset_calendar_time_columns = c(
    GSE143437_muscle_regeneration = "time_days",
    GSE253000_ulipstic_mLN_unsorted = "time_hours",
    GSE48968_BMDC_WT = "time_hours"
  ),
  dataset_calendar_time_units = c(
    GSE143437_muscle_regeneration = "day",
    GSE253000_ulipstic_mLN_unsorted = "hour",
    GSE48968_BMDC_WT = "hour"
  ),

  dataset_calendar_observation_columns = c(GSE48968_BMDC_WT = "experiment_id"),
  dataset_calendar_condition_columns = c(GSE48968_BMDC_WT = "stimulation"),
  dataset_calendar_reference_conditions = c(GSE48968_BMDC_WT = "none"),
  dataset_celltype_columns = c(GSE48968_BMDC_WT = "cell_type"),

  # Confirmed sample roles: the liver and WNT sample labels encode time/sorting.
  dataset_sample_columns = c(
    GSE48968_BMDC_WT = NA_character_,
    GSE141259_lung_whole = "sample_id",
    GSE141259_epithelial = "sample_id",
    GSE253000_ulipstic_mLN_unsorted = "sample_id",
    larry_cytokine_perturbation = NA_character_,
    human_fetal_liver_organoid = NA_character_,
    GSE216451_FGFR2_P7 = NA_character_,
    GSE216451_FGFR2_P21 = NA_character_,
    GSE216451_FGFR2_P7_secretory = NA_character_,
    GSE216451_FGFR2_P21_secretory = NA_character_,
    GSE216451_FGFR2_P21_epithelial = NA_character_,
    GSE270807_HLBO = NA_character_,
    GSE207889_ULA = "sample_id",
    GSE302266_kidney = NA_character_,
    GSE302266_kidney_author_cells = NA_character_,
    GSE300309_DMSO = NA_character_,
    GSE300309_DAPT = NA_character_,
    GSE118184_Takasato_timecourse = NA_character_,
    GSE181595_kidney = "sample_id",
    GSE239582_FT_WNT7A_FZD5 = "patient_id",
    GSE124472_kidney = NA_character_,
    GSE184928_kidney = NA_character_,
    GSE263556_kidney = NA_character_,
    GSE164564_kidney = NA_character_,
    GSE274389_signal_recorder = NA_character_,
    GSE135893_IPF = "sample_id",
    GSE135893_epithelial_receiver = "sample_id",
    GSE135893_alveolar_receiver = "sample_id",
    monocle3_worm_l2 = NA_character_,
    celegans_embryo = NA_character_,
    GSE143437_muscle_regeneration = "sample_id",
    GSE143437_phagocyte_receiver = "sample_id",
    GSE143437_myogenic_receiver = "sample_id",
    lrsim_shared_trend = "donor", lrsim_focal_branch = "donor",
    lrsim_donor_confounded = "donor", lrsim_dropout = "donor",
    lrsim_nested_fork = "donor", lrsim_sampling_thinning = "donor",
    lrsim_transient_focal = "donor", lrsim_donor_ligand_only = "donor",
    lrsim_donor_receptor_only = "donor", lrsim_donor_nested_fork = "donor",
    lrsim_donor_transient_focal = "donor"
  ),

  # The CellAgentChat SCC download supplies processed expression and no raw layer.
  dataset_lr_input_scale = c(cellagentchat_scc = "provided_expression",
    GSE48968_BMDC_WT = "tpm"),

  dataset_species = list(
    lung = "mouse", paul15 = "mouse", pancreas = "mouse", dentategyrus = "mouse",
    cellrank_stemnet = "mouse",
    cellagentchat_scc = "human", monocle3_worm_l2 = "celegans",
    gastrulation_erythroid = "mouse",
    zebrafish_urd = "zebrafish", zebrafish_neural_dev = "zebrafish",
    celegans_embryo = "celegans", human_cd34 = "human",
    human_fetal_liver_organoid = "human", human_fetal_pancreas = "human",
    GSE270807_HLBO = "human",
    GSE207889_ULA = "human",
    GSE302266_kidney = "human",
    GSE302266_kidney_author_cells = "human",
    GSE300309_DMSO = "human",
    GSE300309_DAPT = "human",
    GSE118184_Takasato_timecourse = "human",
    GSE181595_kidney = "human",
    GSE124472_kidney = "human",
    GSE184928_kidney = "human",
    GSE263556_kidney = "human",
    GSE164564_kidney = "human",
    zebrafish_tracerseq = "zebrafish",
    larry_invitro = "mouse", larry_cytokine_perturbation = "mouse",
    GSE253000_ulipstic_mLN_unsorted = "mouse", GSE48968_BMDC_WT = "mouse", simulated_branch = "simulated",
    sim_2br_500c = "simulated", sim_2br_2000c = "simulated", sim_2br_5000c = "simulated",
    sim_3br_500c = "simulated", sim_3br_2000c = "simulated", sim_3br_5000c = "simulated",
    sim_4br_500c = "simulated", sim_4br_2000c = "simulated", sim_4br_5000c = "simulated",
    sim_3br_2000c_batch = "simulated", sim_3br_2000c_drop30 = "simulated",
    sim_3br_2000c_drop50 = "simulated", sim_3br_2000c_drop70 = "simulated",
    sim_3br_2000c_drop90 = "simulated", simulated_complex_3branch = "simulated",
    GSE125527_IBD = "human", GSE132465_CRC = "human", GSE135893_IPF = "human",
    GSE171524_COVID = "human", GSE243917_RA = "human", GSE210152_athero = "human",
    GSE136689_foregut_DE_SM = "mouse",
    GSE115956_intestinal_organoid = "mouse",
    GSE239582_FT_WNT7A_FZD5 = "human",
    GSE274389_signal_recorder = "mouse",
    GSE216451_FGFR2_P7 = "mouse",
    GSE216451_FGFR2_P21 = "mouse",
    GSE216451_FGFR2_P7_secretory = "mouse",
    GSE216451_FGFR2_P21_secretory = "mouse",
    GSE216451_FGFR2_P21_epithelial = "mouse",
    GSE143437_muscle_regeneration = "mouse",
    GSE143437_phagocyte_receiver = "mouse",
    GSE143437_myogenic_receiver = "mouse",
    GSE265803_mdx5cv_quadriceps = "mouse",
    GSE265803_mdx5cv_diaphragm = "mouse",
    GSE265803_Ccr2KO_quadriceps = "mouse",
    GSE265803_Ccr2KO_diaphragm = "mouse",
    GSE135893_epithelial_receiver = "human",
    GSE135893_alveolar_receiver = "human",
    lrsim_shared_trend = "human", lrsim_focal_branch = "human",
    lrsim_donor_confounded = "human", lrsim_dropout = "human",
    lrsim_sampling_thinning = "human", lrsim_nested_fork = "human",
    lrsim_transient_focal = "human", lrsim_donor_ligand_only = "human",
    lrsim_donor_receptor_only = "human", lrsim_donor_nested_fork = "human",
    lrsim_donor_transient_focal = "human"
  )
)

BENCHMARK_PANEL_NAME <- "independent_benchmark"
BENCHMARK_EXCLUSIONS_PATH <- CONFIG$dir$benchmark_exclusions
BENCHMARK_EXCLUSIONS <- c(
  cellrank_stemnet = "processed pancreas subset: all 2531 cells and 3000 genes overlap; identical nonzero expression mask",
  smoke_human_cd34_500c = "500-cell smoke subset of the bone-marrow panel",
  cellrank_bone_marrow = "same 5,780-cell object as human_cd34 and bonemarrow_scvelo",
  bonemarrow_scvelo = "same 5,780-cell object as human_cd34 and cellrank_bone_marrow"
)
BENCHMARK_PSEUDOTIME_ALGORITHMS <- c(
  "LineaLinkR",
  CONFIG$algorithm_groups$ti_trajectory
)
BENCHMARK_ORDERING_DATASETS <- c(
  "lung", "dentategyrus", "gastrulation_erythroid", "zebrafish_neural_dev",
  "larry_invitro", "simulated_branch",
  "sim_2br_500c", "sim_2br_2000c", "sim_2br_5000c",
  "sim_3br_500c", "sim_3br_2000c", "sim_3br_5000c",
  "sim_4br_500c", "sim_4br_2000c", "sim_4br_5000c",
  "sim_3br_2000c_batch", "sim_3br_2000c_drop30",
  "sim_3br_2000c_drop50", "sim_3br_2000c_drop70",
  "sim_3br_2000c_drop90", "simulated_complex_3branch",
  "human_fetal_liver_organoid", "human_fetal_pancreas",
  "zebrafish_tracerseq"
)
COMPLETE_CORE_MAX_CELLS <- 10000L
BENCHMARK_ORDERING_N_OBS <- c(
  lung = 32559L,
  dentategyrus = 2930L,
  gastrulation_erythroid = 9815L,
  zebrafish_neural_dev = 8331L,
  larry_invitro = 49302L,
  simulated_branch = 2000L,
  sim_2br_500c = 500L,
  sim_2br_2000c = 2000L,
  sim_2br_5000c = 5000L,
  sim_3br_500c = 500L,
  sim_3br_2000c = 2000L,
  sim_3br_5000c = 5000L,
  sim_4br_500c = 500L,
  sim_4br_2000c = 2000L,
  sim_4br_5000c = 5000L,
  sim_3br_2000c_batch = 2000L,
  sim_3br_2000c_drop30 = 2000L,
  sim_3br_2000c_drop50 = 2000L,
  sim_3br_2000c_drop70 = 2000L,
  sim_3br_2000c_drop90 = 2000L,
  simulated_complex_3branch = 2000L,
  human_fetal_liver_organoid = 22874L,
  human_fetal_pancreas = 28430L,
  zebrafish_tracerseq = 36749L
)
BENCHMARK_COMPLETE_CORE_DATASETS <- c(
  "dentategyrus", "gastrulation_erythroid", "zebrafish_neural_dev",
  "simulated_branch",
  "sim_2br_500c", "sim_2br_2000c", "sim_2br_5000c",
  "sim_3br_500c", "sim_3br_2000c", "sim_3br_5000c",
  "sim_4br_500c", "sim_4br_2000c", "sim_4br_5000c",
  "sim_3br_2000c_batch", "sim_3br_2000c_drop30",
  "sim_3br_2000c_drop50", "sim_3br_2000c_drop70",
  "sim_3br_2000c_drop90", "simulated_complex_3branch"
)
BENCHMARK_SCALE_EXTENSION_DATASETS <- c(
  "lung", "larry_invitro", "human_fetal_liver_organoid",
  "human_fetal_pancreas", "zebrafish_tracerseq"
)
BENCHMARK_BRANCH_TERMINAL_COUNTS <- c(
  simulated_branch = 2L,
  sim_2br_500c = 2L,
  sim_2br_2000c = 2L,
  sim_2br_5000c = 2L
)
BENCHMARK_BRANCH_TERMINAL_MIN_PSEUDOTIME <- c(
  simulated_branch = 0.4,
  sim_2br_500c = 0.4,
  sim_2br_2000c = 0.4,
  sim_2br_5000c = 0.4
)
BENCHMARK_BRANCH_TERMINAL_EXPECTED_N <- c(
  simulated_branch = 1357L,
  sim_2br_500c = 366L,
  sim_2br_2000c = 1470L,
  sim_2br_5000c = 3684L
)
BENCHMARK_BRANCH_DATASETS <- names(BENCHMARK_BRANCH_TERMINAL_COUNTS)
BENCHMARK_ORDERING_ONLY_SEGMENT_DATASETS <- c(
  "sim_3br_500c", "sim_3br_2000c", "sim_3br_5000c",
  "sim_4br_500c", "sim_4br_2000c", "sim_4br_5000c",
  "sim_3br_2000c_batch", "sim_3br_2000c_drop30",
  "sim_3br_2000c_drop50", "sim_3br_2000c_drop70",
  "sim_3br_2000c_drop90", "simulated_complex_3branch"
)
BENCHMARK_FORCE_AUTO_CLUSTER_DATASETS <- BENCHMARK_ORDERING_ONLY_SEGMENT_DATASETS
stopifnot(
  identical(names(BENCHMARK_BRANCH_TERMINAL_MIN_PSEUDOTIME), BENCHMARK_BRANCH_DATASETS),
  identical(names(BENCHMARK_BRANCH_TERMINAL_EXPECTED_N), BENCHMARK_BRANCH_DATASETS),
  all(BENCHMARK_BRANCH_TERMINAL_EXPECTED_N > 0L),
  all(BENCHMARK_BRANCH_DATASETS %in% BENCHMARK_ORDERING_DATASETS),
  all(BENCHMARK_ORDERING_ONLY_SEGMENT_DATASETS %in% BENCHMARK_ORDERING_DATASETS),
  identical(BENCHMARK_FORCE_AUTO_CLUSTER_DATASETS, BENCHMARK_ORDERING_ONLY_SEGMENT_DATASETS),
  !any(BENCHMARK_BRANCH_DATASETS %in% BENCHMARK_ORDERING_ONLY_SEGMENT_DATASETS)
)

InitDirs <- function() {
  for (d in list(
    CONFIG$dir$table1_main, CONFIG$dir$reference, CONFIG$dir$evaluation_audit,
    CONFIG$dir$supplement_tables,
    CONFIG$dir$linealinkr, CONFIG$dir$canonical,
    CONFIG$dir$log, CONFIG$dir$algorithm_out, CONFIG$dir$algorithm_misc,
    CONFIG$dir$r_export
  )) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
}

BenchmarkTablePath <- function(...) {
  rel <- file.path(...)
  if (length(list(...)) >= 1L) {
    first <- list(...)[[1L]]
    if (is.character(first) && grepl("^table1", first)) {
      return(file.path(CONFIG$dir$table1_main, ...))
    }
  }
  file.path(CONFIG$dir$supplement_tables, ...)
}

StageScript <- function(basename, tier = c("manuscript", "exploratory")) {
  tier <- match.arg(tier)
  roots <- c(
    file.path(CONFIG$dir$pipeline, c("run", "evaluation", "biology", "figures", "tests", "ablation", "data")),
    file.path(CONFIG$dir$root, "research", c("biology", "geometry", "communication", "time", "root", "lineage"))
  )
  candidates <- unlist(lapply(roots, list.files, recursive = TRUE, full.names = TRUE))
  hits <- candidates[base::basename(candidates) == basename]
  if (length(hits) != 1L) stop("Expected one exact script name: ", basename)
  hits[[1L]]
}

BuildDatasetRegistry <- function() {
  source_relative_paths <- c(
    setNames(file.path("curated", unlist(CONFIG$datasets_curated, use.names = FALSE)),
             names(CONFIG$datasets_curated)),
    setNames(file.path("literature_benchmarks", unlist(CONFIG$datasets_literature, use.names = FALSE)),
             names(CONFIG$datasets_literature)),
    setNames(file.path("canonical", unlist(CONFIG$datasets_simulated, use.names = FALSE)),
             names(CONFIG$datasets_simulated)),
    setNames(unlist(CONFIG$datasets_geo_candidate, use.names = FALSE),
             names(CONFIG$datasets_geo_candidate))
  )
  project_paths <- setNames(
    file.path(CONFIG$dir$canonical, unlist(CONFIG$datasets_project_canonical, use.names = FALSE)),
    names(CONFIG$datasets_project_canonical)
  )
  revision_paths <- c(
    setNames(file.path(READONLY_REVISION_ROOT, "sim", unlist(CONFIG$datasets_revision_sim, use.names = FALSE)),
             names(CONFIG$datasets_revision_sim)),
    setNames(file.path(READONLY_REVISION_ROOT, unlist(CONFIG$datasets_revision_author_cohort, use.names = FALSE)),
             names(CONFIG$datasets_revision_author_cohort))
  )
  registry_names <- c(names(source_relative_paths), names(project_paths), names(revision_paths))
  if (!length(registry_names) || anyDuplicated(registry_names)) {
    stop("Active dataset registry must contain unique independent dataset names", call. = FALSE)
  }
  source_absolute_paths <- vapply(
    file.path(CONFIG$dir$data, unname(source_relative_paths)),
    normalizePath,
    character(1L),
    winslash = "/",
    mustWork = FALSE
  )
  names(source_absolute_paths) <- names(source_relative_paths)
  project_absolute_paths <- vapply(
    unname(project_paths),
    normalizePath,
    character(1L),
    winslash = "/",
    mustWork = FALSE
  )
  names(project_absolute_paths) <- names(project_paths)
  revision_absolute_paths <- vapply(
    unname(revision_paths),
    normalizePath,
    character(1L),
    winslash = "/",
    mustWork = FALSE
  )
  names(revision_absolute_paths) <- names(revision_paths)
  absolute_paths <- c(source_absolute_paths, project_absolute_paths, revision_absolute_paths)
  if (anyDuplicated(tolower(absolute_paths))) {
    stop("Active dataset registry contains duplicate canonical paths", call. = FALSE)
  }
  absolute_paths
}

BuildNarrativeDatasetRegistry <- function() {
  source_paths <- setNames(
    file.path(CONFIG$dir$data, unlist(CONFIG$datasets_geo_narrative, use.names = FALSE)),
    names(CONFIG$datasets_geo_narrative)
  )
  revision_paths <- setNames(
    file.path(READONLY_REVISION_ROOT, unlist(CONFIG$datasets_revision_narrative, use.names = FALSE)),
    names(CONFIG$datasets_revision_narrative)
  )
  external_paths <- c(
    larry_cytokine_perturbation = file.path(WORK_ROOT,
      "inputs/external_validation/LARRY/larry_cytokine_perturbation.h5ad"),
    GSE253000_ulipstic_mLN_unsorted = file.path(WORK_ROOT,
      "inputs/external_validation/GSE253000/GSE253000_ulipstic_mLN_unsorted.h5ad"),
    GSE48968_BMDC_WT = file.path(WORK_ROOT,
      "inputs/external_validation/GSE48968/GSE48968_BMDC_WT.h5ad"))
  absolute_paths <- vapply(
    c(source_paths, revision_paths, external_paths),
    normalizePath,
    character(1L),
    winslash = "/",
    mustWork = FALSE
  )
  absolute_paths
}

CONFIG$dataset_registry <- BuildDatasetRegistry()
CONFIG$dataset_narrative_registry <- BuildNarrativeDatasetRegistry()
BENCHMARK_NON_ORDERING_DATASETS <- setdiff(
  names(CONFIG$dataset_registry),
  BENCHMARK_ORDERING_DATASETS
)
stopifnot(
  length(CONFIG$dataset_registry) == 36L + length(CONFIG$datasets_revision_sim),
  length(BENCHMARK_PSEUDOTIME_ALGORITHMS) == 12L,
  !anyDuplicated(BENCHMARK_PSEUDOTIME_ALGORITHMS),
  identical(
    setdiff(BENCHMARK_PSEUDOTIME_ALGORITHMS, "LineaLinkR"),
    CONFIG$algorithm_groups$ti_trajectory
  ),
  length(BENCHMARK_ORDERING_DATASETS) == 24L,
  length(BENCHMARK_NON_ORDERING_DATASETS) == 12L + length(CONFIG$datasets_revision_sim),
  !any(BENCHMARK_NON_ORDERING_DATASETS %in% BENCHMARK_ORDERING_DATASETS),
  setequal(
    c(BENCHMARK_ORDERING_DATASETS, BENCHMARK_NON_ORDERING_DATASETS),
    names(CONFIG$dataset_registry)
  ),
  identical(names(BENCHMARK_ORDERING_N_OBS), BENCHMARK_ORDERING_DATASETS),
  all(is.finite(BENCHMARK_ORDERING_N_OBS)),
  all(BENCHMARK_ORDERING_N_OBS > 0L),
  length(BENCHMARK_COMPLETE_CORE_DATASETS) == 19L,
  length(BENCHMARK_SCALE_EXTENSION_DATASETS) == 5L,
  identical(
    BENCHMARK_COMPLETE_CORE_DATASETS,
    BENCHMARK_ORDERING_DATASETS[
      BENCHMARK_ORDERING_N_OBS[BENCHMARK_ORDERING_DATASETS] <= COMPLETE_CORE_MAX_CELLS
    ]
  ),
  identical(
    BENCHMARK_SCALE_EXTENSION_DATASETS,
    BENCHMARK_ORDERING_DATASETS[
      BENCHMARK_ORDERING_N_OBS[BENCHMARK_ORDERING_DATASETS] > COMPLETE_CORE_MAX_CELLS
    ]
  ),
  identical(
    BENCHMARK_ORDERING_DATASETS,
    BENCHMARK_ORDERING_DATASETS[
      BENCHMARK_ORDERING_DATASETS %in% c(
        BENCHMARK_COMPLETE_CORE_DATASETS,
        BENCHMARK_SCALE_EXTENSION_DATASETS
      )
    ]
  ),
  length(BENCHMARK_BRANCH_DATASETS) > 0L,
  all(BENCHMARK_ORDERING_DATASETS %in% names(CONFIG$dataset_registry)),
  all(BENCHMARK_BRANCH_DATASETS %in% BENCHMARK_ORDERING_DATASETS)
)

BuildDatasetList <- function() {
  as.list(CONFIG$dataset_registry)
}

DatasetPath <- function(ds_name) {
  registered_paths <- c(CONFIG$dataset_registry, CONFIG$dataset_narrative_registry)
  p <- registered_paths[[ds_name]]
  if (is.null(p)) {
    stop("Dataset is not explicitly registered: ", ds_name, call. = FALSE)
  }
  normalizePath(p,winslash="/",mustWork=TRUE)
}

inherited_ablation <- Sys.getenv("LLR_INHERITED_ABLATION", "full")
stopifnot(inherited_ablation %in% c("full", "graph_k_10", "graph_k_30",
  "no_kmeans_refinement"))
CONFIG$llr$knn_neighbors <- switch(inherited_ablation,
  graph_k_10 = 10L, graph_k_30 = 30L, CONFIG$llr$knn_neighbors)
CONFIG$llr$force_kmeans_refine <- CONFIG$llr$force_kmeans_refine &&
  inherited_ablation != "no_kmeans_refinement"
CONFIG$llr$inherited_ablation <- inherited_ablation

InitDirs()
datasets <- BuildDatasetList()
