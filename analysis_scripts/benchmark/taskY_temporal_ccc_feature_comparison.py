# -*- coding: utf-8 -*-
"""Primary-source feature comparison for temporal CCC methods."""
from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_LLR_CFG = Path(__file__).resolve().parents[1] / "config"
sys.path.insert(0, str(_LLR_CFG))
from llr_paths import assets_dir as _llr_assets_dir  # noqa: E402

_LLR_ASSETS = _llr_assets_dir()

CONFIG = {
    "paths": {
        "output": _LLR_ASSETS
        / "results/paper_tables/supplement/tables/temporal_ccc_feature_comparison"
    },
    "threads": os.cpu_count() or 1,
    "seed": 12345,
    "plot": {},
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
plt.rcParams.update(CONFIG["plot"])
CONFIG["paths"]["output"].mkdir(parents=True, exist_ok=True)

comparison = pd.DataFrame(
    {
        "method": ["LineaLinkR", "TICCI", "TraSig", "TimeTalk"],
        "primary_input": [
            "snapshot scRNA-seq; cell-type labels optional but used for CCC",
            "snapshot scRNA-seq",
            "expression plus a pseudotime trajectory with edges or groups",
            "time-course scRNA-seq integrated along pseudotime",
        ],
        "trajectory_handling": [
            "internal blind-root TI or user root; lineage weights retained",
            "internal Louvain/PAGA, scEntropy direction, directed tree and fitted time",
            "external TI output; paper uses CSHMM and permits other TI methods",
            "trajectory reconstructed before LR/TF analysis; paper uses Monocle workflows",
        ],
        "branch_context": [
            "yes; source-target-lineage-bin records",
            "yes for trajectory topology; not an LR branch profile",
            "yes; sliding-window profiles on trajectory edges",
            "cell-type/pseudotime context; no bifurcation-relative LR role",
        ],
        "repeated_timepoints_required": [
            "no",
            "no",
            "no; validity depends on trusted pseudotime",
            "yes for the intended time-course workflow",
        ],
        "native_lr_output": [
            "saturating LR activity by source, target, lineage and bin",
            "CCI-weighted cell-neighborhood graph used to improve TI; no bin-resolved LR ranking",
            "sliding-window ligand-receptor profile dot products with permutation significance",
            "co-varying LR interaction scores linked to temporal TF activity",
        ],
        "lag_or_direction_claim": [
            "none; descriptive expression-lag audit is external to the score",
            "trajectory direction from scEntropy; no ligand-receptor lag test",
            "synchronous profile matching; no ligand-leading causal test",
            "Granger direction between LR activity and TF clusters; authors require further validation",
        ],
        "bifurcation_relative_role": [
            "five operational labels; pilot calibration only",
            "no",
            "no",
            "no; forward/backward/feedback labels refer to LR-TF Granger direction",
        ],
        "spatial_input": ["no", "no", "no", "no"],
        "native_statistical_support": [
            "score and structural/expression flags; no event-level causal p-value",
            "trajectory benchmarks and marker/stage concordance",
            "permutation tests for temporal profile matches",
            "Granger tests and LR-TF correlation filters",
        ],
        "primary_source": [
            "current LineaLinkR protocol and source code",
            "Fu et al., Bioinformatics 2025, doi:10.1093/bioinformatics/btaf027",
            "Li et al., Genome Biology 2022, doi:10.1186/s13059-022-02629-7",
            "Wang et al., Communications Biology 2023, doi:10.1038/s42003-023-05283-2",
        ],
        "audit_interpretation": [
            "integrated snapshot temporal-CCC protocol; not causal TI or spatial validation",
            "CCI-informed TI comparator, not a temporal LR-scoring substitute",
            "closest comparator for pseudotime-resolved LR profile matching",
            "stronger LR-TF temporal-direction model but needs repeated temporal sampling and sensitivity checks",
        ],
    }
)

comparison.to_csv(
    CONFIG["paths"]["output"] / "temporal_ccc_method_feature_comparison.csv",
    index=False,
)

print(comparison.to_string(index=False))
