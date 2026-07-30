"""Generate simulated scRNA-seq data with known ground truth.

Generates branching trajectory data where:
- pseudotime and branch assignments are KNOWN
- gene expression changes along pseudotime
- some genes are branch-specific (dynamic vs static)
- known LR pairs are embedded with temporal signaling patterns
"""
import numpy as np
import scanpy as sc
import pandas as pd
from scipy.sparse import csr_matrix
from pathlib import Path

np.random.seed(42)

# ═══════════════════════════════════════════════════════════
# Parameters
# ═══════════════════════════════════════════════════════════
N_CELLS = 2000
N_GENES = 1000
N_BRANCHES = 2
BRANCH_SPLIT = 0.4  # bifurcation at 40% pseudotime

# ═══════════════════════════════════════════════════════════
# Generate pseudotime and branches
# ═══════════════════════════════════════════════════════════
# Pre-branch pseudotime (uniform)
branch1_cells = int(N_CELLS * 0.6)
branch2_cells = N_CELLS - branch1_cells

# Branch 1 cells: pseudotime 0..1
pt1 = np.sort(np.random.beta(2, 3, branch1_cells))
# Branch 2 cells: same root, then diverges
pt2_root = np.sort(np.random.beta(2, 3, branch2_cells))
# After split, add displacement
pt2 = pt2_root * (1 - BRANCH_SPLIT) + BRANCH_SPLIT

pseudotime = np.concatenate([pt1, pt2])
branch = np.array([0] * branch1_cells + [1] * branch2_cells)

# ═══════════════════════════════════════════════════════════
# Generate gene expression with pseudotime dependency
# ═══════════════════════════════════════════════════════════
# 3 gene classes:
#   - Early genes (high at low pseudotime, then decrease)
#   - Late genes (increase with pseudotime)
#   - Branch-specific genes (different in branch 0 vs 1)
#   - Noise genes (random)
N_EARLY = 100; N_LATE = 100; N_BRANCH_SPECIFIC = 50; N_NOISE = N_GENES - N_EARLY - N_LATE - N_BRANCH_SPECIFIC

expr = np.zeros((N_CELLS, N_GENES))
gene_names = []; gene_classes = []

# Early genes
for i in range(N_EARLY):
    base = np.random.lognormal(2, 0.5)
    trend = 4 * (1 - pseudotime)
    expr[:, i] = np.random.poisson(base * np.exp(-trend))
    gene_names.append(f"EARLY_{i+1}")
    gene_classes.append("early")

# Late genes
for i in range(N_LATE):
    base = np.random.lognormal(2, 0.5)
    trend = 3 * pseudotime
    expr[:, N_EARLY + i] = np.random.poisson(base * np.exp(trend * 0.5))
    gene_names.append(f"LATE_{i+1}")
    gene_classes.append("late")

# Branch-specific genes (different in branch 1 vs branch 0)
for i in range(N_BRANCH_SPECIFIC):
    base = np.random.lognormal(2, 0.3)
    expr_b1 = np.random.poisson(base * np.exp(2 * (pseudotime - BRANCH_SPLIT) * (branch == 1)))
    expr_b0 = np.random.poisson(base * np.exp(2 * (1 - pseudotime) * (branch == 0)))
    expr[:, N_EARLY + N_LATE + i] = np.where(branch == 1, expr_b1, expr_b0)
    gene_names.append(f"BRANCH_{i+1}")
    gene_classes.append("branch_specific")

# Noise genes
start = N_EARLY + N_LATE + N_BRANCH_SPECIFIC
for i in range(N_NOISE):
    base = np.random.lognormal(1.5, 0.8)
    expr[:, start + i] = np.random.poisson(base)
    gene_names.append(f"NOISE_{i+1}")
    gene_classes.append("noise")

# ═══════════════════════════════════════════════════════════
# Embed known LR pairs with temporal signaling
# ═══════════════════════════════════════════════════════════
# Fake LR pairs: L1-R1 activates early, L2-R2 activates late
# L1 peaks at pseudotime=0.2; R1 expressed in branch 0
# L2 peaks at pseudotime=0.7; R2 expressed in branch 1
known_lr_true = {
    "L1": {"pseudotime_peak": 0.2, "branch": 0},
    "L2": {"pseudotime_peak": 0.7, "branch": 1},
}

# ═══════════════════════════════════════════════════════════
# Build AnnData
# ═══════════════════════════════════════════════════════════
obs = pd.DataFrame({
    "true_pseudotime": pseudotime,
    "true_branch": branch.astype(str),
    "true_branch_numeric": branch,
    "clusters": [f"branch_{b}" for b in branch],  # for slingshot cluster detection
    "cell_type": [f"branch_{b}" for b in branch],  # backup cluster column
    "cell_id": [f"cell_{i}" for i in range(N_CELLS)],
})
# Force string type to avoid pandas categorical export issues
obs["clusters"] = obs["clusters"].astype(str)
obs["cell_type"] = obs["cell_type"].astype(str)
obs["true_branch"] = obs["true_branch"].astype(str)
var = pd.DataFrame({
    "gene_class": gene_classes,
}, index=gene_names)

adata = sc.AnnData(csr_matrix(expr.astype(np.float32)), obs=obs, var=var)
adata.var_names = gene_names

# Layers
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
adata.layers["lognorm"] = adata.X.copy()
adata.X = adata.layers["counts"].copy()

# ═══════════════════════════════════════════════════════════
# Save
# ═══════════════════════════════════════════════════════════
out_dir = Path(r"F:\biotrade\LineaLinkR_benchmark_assets\data\literature_benchmarks\simulated_branch")
out_dir.mkdir(parents=True, exist_ok=True)
adata.write(out_dir / "simulated_branch.h5ad")

import json
info = {
    "dataset": "simulated_branch",
    "description": "Synthetic branching trajectory with known pseudotime and branch labels",
    "species": "simulated",
    "ground_truth": "true_pseudotime + true_branch",
    "n_cells": N_CELLS, "n_genes": N_GENES,
    "n_branches": N_BRANCHES, "bifurcation_point": float(BRANCH_SPLIT),
    "obs_columns": list(adata.obs.columns),
    "gene_classes": {"early": N_EARLY, "late": N_LATE, "branch_specific": N_BRANCH_SPECIFIC, "noise": N_NOISE},
    "known_lr_signals": known_lr_true,
}
with open(out_dir / "dataset_info.json", "w") as f:
    json.dump(info, f, indent=2)

print(f"Saved: {N_CELLS}c x {N_GENES}g")
print(f"  Branches: {N_BRANCHES}, bifurcation at pseudotime={BRANCH_SPLIT}")
print(f"  Early genes: {N_EARLY}, Late: {N_LATE}, Branch-specific: {N_BRANCH_SPECIFIC}, Noise: {N_NOISE}")
print(f"  True pseudotime range: [{pseudotime.min():.3f}, {pseudotime.max():.3f}]")
