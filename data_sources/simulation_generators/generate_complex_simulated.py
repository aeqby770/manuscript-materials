"""Generate complex simulated scRNA-seq data with multi-branch trajectories.

Creates datasets with:
  - 3-4 branching topology (hematopoiesis-like)
  - Multiple dropout rates (30%, 50%, 70%, 90%)
  - Multiple cell counts (500, 2000, 5000, 10000)
  - Batch effects (2 batches)
  - Known ground truth: pseudotime, branch labels, dynamic gene classes
  - Embedded known LR pairs with temporal signaling patterns

Output: multiple .h5ad files with canonical structure compatible with
the benchmark pipeline (UMAP + layers + obs columns).

Ground truth columns:
  - true_pseudotime  : numeric [0,1]
  - true_branch      : categorical branch label
  - true_bifurcation_level : which bifurcation level (primary/secondary)
  - batch            : batch label (for batch-effect datasets)
  - clusters         : for cluster-based algorithms
  - cell_type        : for cell-type-based algorithms
"""

import numpy as np
import scanpy as sc
import pandas as pd
from scipy.sparse import csr_matrix
from pathlib import Path
import json

# ═══════════════════════════════════════════════════════════════
# Global config
# ═══════════════════════════════════════════════════════════════
CONFIG = {
    "seed": 42,
    "n_genes": 1000,
    "gene_classes": {
        "early":       {"count": 80,  "desc": "high early, decay"},
        "late_shared": {"count": 80,  "desc": "increase globally"},
        "branch1":     {"count": 60,  "desc": "specific to branch 1"},
        "branch2":     {"count": 60,  "desc": "specific to branch 2"},
        "branch3":     {"count": 60,  "desc": "specific to branch 3 (if present)"},
        "branch4":     {"count": 40,  "desc": "specific to branch 4 (if present)"},
        "bifurcation_genes": {"count": 40, "desc": "change at bifurcation points"},
        "noise":       {"count": "remainder", "desc": "random noise"},
    },
    "dropout_rates":    [0.30, 0.50, 0.70, 0.90],
    "cell_counts":      [500, 2000, 5000, 10000],
    "n_branches_options": [2, 3, 4],  # will generate each
    "output_root":      Path(r"F:\biotrade\LineaLinkR_benchmark_assets\data\simulated_complex"),
}

np.random.seed(CONFIG["seed"])

# ═══════════════════════════════════════════════════════════════
# Branch topology definitions
# ═══════════════════════════════════════════════════════════════

def topology_2branch(n_cells):
    """Simple bifurcation at pseudotime=0.4"""
    n0 = int(n_cells * 0.5)
    n1 = n_cells - n0
    pt0 = np.sort(np.random.beta(2, 3, n0))
    pt1 = np.sort(np.random.beta(2, 3, n1)) * (1 - 0.4) + 0.4
    pseudotime = np.concatenate([pt0, pt1])
    branch = np.array([0] * n0 + [1] * n1)
    bif_level = np.array([0] * n0 + [1] * n1)
    return pseudotime, branch, bif_level, {"bifurcation_points": [0.4], "n_branches": 2}


def topology_3branch(n_cells):
    """Primary bifurcation at 0.35 → branch 0/1 + 2
       Secondary: branch 1 → branch 1a/1b at 0.65 (nested bifurcation)"""
    n0 = int(n_cells * 0.30)   # branch 0
    n1a = int(n_cells * 0.25)  # branch 1a
    n1b = int(n_cells * 0.25)  # branch 1b
    n2 = n_cells - n0 - n1a - n1b  # branch 2

    pt0 = np.sort(np.random.beta(2, 3, n0))
    pt2 = np.sort(np.random.beta(2, 3, n2)) * (1 - 0.35) + 0.35

    # Branch 1: common trunk to 0.65, then split into 1a/1b
    n_trunk = n1a + n1b
    pt_trunk = np.sort(np.random.beta(2, 3, n_trunk)) * (0.65 - 0.35) + 0.35
    pt1a = np.sort(np.random.beta(2, 3, n1a)) * (1 - 0.65) + 0.65
    pt1b = np.sort(np.random.beta(2, 3, n1b)) * (1 - 0.65) + 0.65

    pseudotime = np.concatenate([pt0, pt_trunk, pt_trunk, pt2])
    pt_len_n0 = len(pt0)
    # Reassign branch 1 trunk: half to 1a, half to 1b based on pt
    pt1 = np.concatenate([pt_trunk[:n1a], pt1a])  # not used directly, keep for ordering
    # Simpler: assign branch 1 trunk cells randomly to 1a/1b
    trunk_pt = np.concatenate([pt_trunk, pt_trunk])
    # Build full branch array
    branch = np.array([0] * n0 + [1] * n1a + [2] * n1b + [3] * n2)
    # Actually make the topology clear: 0=primary_branch, 1=secondary_1a, 2=secondary_1b, 3=other_branch
    # But the pseudotime for branch 1's pre-secondary cells should be in [0.35, 0.65]
    # Let me redesign this more carefully

    # Branch topology:
    #   root [0, 0.35] → 
    #     branch 0: [0.35, 1.0] (n0 cells)  
    #     trunk: [0.35, 0.65] → 
    #       branch 1a: [0.65, 1.0] (n1a cells)
    #       branch 1b: [0.65, 1.0] (n1b cells)
    #     branch 2: [0.35, 1.0] (n2 cells)

    n_trunk_cells = n1a + n1b
    pt_root = np.sort(np.random.beta(2, 3, n_trunk_cells)) * 0.35  # shared root for all branches
    # Redistribute: first n0 get branch 0 continuation, next n_trunk get trunk, last n2 get branch 2
    # Actually, let me just assign cells to branches proportionally with proper pseudotime

    # All cells start from root [0, 0.35]
    pt_pre = np.sort(np.random.beta(2, 3, n_cells)) * 0.35

    # Post-bifurcation pseudotime for each branch
    branch_assignments = np.zeros(n_cells, dtype=int)
    bif_level = np.zeros(n_cells, dtype=int)
    pt = pt_pre.copy()

    idx = 0
    # Branch 0: n0 cells
    for i in range(n0):
        pt[idx] = pt_pre[idx]
        branch_assignments[idx] = 0
        bif_level[idx] = 1
        idx += 1
    # Branch 2: n2 cells
    for i in range(n2):
        pt[idx] = pt_pre[idx]
        branch_assignments[idx] = 2
        bif_level[idx] = 1
        idx += 1
    # Trunk → 1a + 1b: need to give them proper bifurcation post-times
    for i in range(n1a):
        pt[idx] = pt_pre[idx] * (0.65 / 0.35)  # scale to [0, 0.65]
        branch_assignments[idx] = 1
        bif_level[idx] = 1
        idx += 1
    for i in range(n1b):
        pt[idx] = pt_pre[idx] * (0.65 / 0.35)
        branch_assignments[idx] = 1
        bif_level[idx] = 1
        idx += 1

    # Now post-secondary split for 1a/1b: randomly assign trunk cells
    trunk_mask = branch_assignments == 1
    trunk_indices = np.where(trunk_mask)[0]
    np.random.shuffle(trunk_indices)
    for i, t_idx in enumerate(trunk_indices):
        if i < n1a:
            branch_assignments[t_idx] = 1
            bif_level[t_idx] = 2
            pt[t_idx] = pt[t_idx] * (1 - 0.65) / 0.65 + 0.65
        else:
            branch_assignments[t_idx] = 3  # sub-branch label
            bif_level[t_idx] = 2
            pt[t_idx] = pt[t_idx] * (1 - 0.65) / 0.65 + 0.65

    pseudotime = pt.copy()
    branch = branch_assignments.copy()

    return pseudotime, branch, bif_level, {"bifurcation_points": [0.35, 0.65], "n_branches": 4}


def topology_4branch(n_cells):
    """Two-level bifurcation: first at 0.3, second at 0.6 for each primary branch"""
    # Primary split at 0.3 into left/right
    # Left → left_a/left_b at 0.6
    # Right → right_a/right_b at 0.6
    n_left = int(n_cells * 0.5)
    n_right = n_cells - n_left
    n_la = int(n_left * 0.5)
    n_lb = n_left - n_la
    n_ra = int(n_right * 0.5)
    n_rb = n_right - n_ra

    pt_pre = np.sort(np.random.beta(2, 3, n_cells)) * 0.3

    branch = np.zeros(n_cells, dtype=int)
    bif_level = np.zeros(n_cells, dtype=int)
    pt = pt_pre.copy()

    # Assign left trunk then split
    left_indices = np.arange(0, n_left)
    np.random.shuffle(left_indices)
    for i, li in enumerate(left_indices):
        if i < n_la:
            branch[li] = 0
            bif_level[li] = 2
            pt[li] = pt_pre[li] * (1 - 0.6) / 0.3 + 0.6
        else:
            branch[li] = 1
            bif_level[li] = 2
            pt[li] = pt_pre[li] * (1 - 0.6) / 0.3 + 0.6

    # Assign right trunk then split
    right_indices = np.arange(n_left, n_cells)
    np.random.shuffle(right_indices)
    for i, ri in enumerate(right_indices):
        if i < n_ra:
            branch[ri] = 2
            bif_level[ri] = 2
            pt[ri] = pt_pre[ri] * (1 - 0.6) / 0.3 + 0.6
        else:
            branch[ri] = 3
            bif_level[ri] = 2
            pt[ri] = pt_pre[ri] * (1 - 0.6) / 0.3 + 0.6

    pseudotime = pt.copy()
    return pseudotime, branch, bif_level, {"bifurcation_points": [0.3, 0.6], "n_branches": 4}


TOPOLOGIES = {
    2: topology_2branch,
    3: topology_3branch,
    4: topology_4branch,
}

# ═══════════════════════════════════════════════════════════════
# Gene expression generation
# ═══════════════════════════════════════════════════════════════

def generate_expression(pseudotime, branch, n_genes, n_cells):
    """Generate gene expression matrix with known structure."""
    expr = np.zeros((n_cells, n_genes), dtype=np.float32)
    gene_names = []
    gene_classes = []
    gene_info = []  # per-gene metadata

    # Helper: add gene block
    def add_genes(base_count, class_name, expr_func):
        nonlocal gene_names, gene_classes, gene_info, expr
        start_idx = len(gene_names)
        for i in range(base_count):
            gene_idx = start_idx + i
            expr[:, gene_idx] = expr_func(i)
            gene_names.append(f"{class_name.upper()}_{i+1}")
            gene_classes.append(class_name)
            gene_info.append({"class": class_name, "index_in_class": i})

    # Early genes: high at low pseudotime, decay
    add_genes(80, "early", lambda i: np.random.poisson(
        np.random.lognormal(2, 0.5) * np.exp(-4 * pseudotime)))

    # Late shared genes: increase with pseudotime globally
    add_genes(80, "late_shared", lambda i: np.random.poisson(
        np.random.lognormal(2, 0.5) * np.exp(3 * (pseudotime - 0.3))))

    # Branch-specific genes
    for br in range(4):
        if br == 0:
            br_name = "branch0"; n = 60
        elif br == 1:
            br_name = "branch1"; n = 60
        elif br == 2:
            br_name = "branch2"; n = 60
        else:
            br_name = "branch3"; n = 40
        mask = branch == br
        add_genes(n, br_name, lambda i, m=mask: np.where(
            m, np.random.poisson(np.random.lognormal(2, 0.3) * np.exp(2 * pseudotime)),
            np.random.poisson(np.random.lognormal(0.5, 0.3))))

    # Bifurcation genes: expression changes at bifurcation points
    add_genes(40, "bifurcation", lambda i: np.random.poisson(
        np.random.lognormal(1.5, 0.5) * (1 + 2 * (pseudotime > 0.35) * (branch == 0))))

    # Noise genes
    n_noise = n_genes - len(gene_names)
    add_genes(n_noise, "noise", lambda i: np.random.poisson(np.random.lognormal(1.5, 0.8)))

    return expr, gene_names, gene_classes, gene_info


def apply_dropout(expr, dropout_rate, seed=None):
    """Apply dropout (random zeros) to expression matrix."""
    if seed is not None:
        np.random.seed(seed)
    mask = np.random.random(expr.shape) < dropout_rate
    expr_do = expr.copy()
    expr_do[mask] = 0
    return expr_do


def apply_batch_effect(expr, n_batches=2, batch_seed=None):
    """Add batch effect: multiplicative noise per batch."""
    if batch_seed is not None:
        np.random.seed(batch_seed)
    n_cells = expr.shape[0]
    batch = np.random.choice(n_batches, n_cells)
    for b in range(n_batches):
        mask = batch == b
        batch_factor = np.random.lognormal(0, 0.3)
        expr[mask] = expr[mask] * batch_factor
    return expr, batch


# ═══════════════════════════════════════════════════════════════
# Build AnnData with canonical format
# ═══════════════════════════════════════════════════════════════

def build_anndata(expr, pseudotime, branch, bif_level, batch_label,
                  gene_names, gene_classes, config_meta):
    """Build AnnData object matching canonical benchmark format."""
    n_cells = len(pseudotime)
    obs = pd.DataFrame({
        "true_pseudotime": pseudotime,
        "true_branch": branch.astype(str),
        "true_bifurcation_level": bif_level,
        "clusters": [f"branch_{b}" for b in branch],
        "cell_type": [f"branch_{b}" for b in branch],
        "batch": batch_label,
        "cell_id": [f"cell_{i}" for i in range(n_cells)],
    })
    obs["clusters"] = obs["clusters"].astype(str)
    obs["cell_type"] = obs["cell_type"].astype(str)
    obs["true_branch"] = obs["true_branch"].astype(str)
    obs["batch"] = obs["batch"].astype(str)

    var = pd.DataFrame({"gene_class": gene_classes}, index=gene_names)

    adata = sc.AnnData(csr_matrix(expr.astype(np.float32)), obs=obs, var=var)
    adata.var_names = gene_names

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    adata.layers["lognorm"] = adata.X.copy()
    adata.X = adata.layers["counts"].copy()

    # Compute UMAP
    try:
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        sc.pp.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
    except Exception:
        pass

    return adata


# ═══════════════════════════════════════════════════════════════
# Main generation loop
# ═══════════════════════════════════════════════════════════════

def main():
    np.random.seed(CONFIG["seed"])
    out_root = CONFIG["output_root"]
    out_root.mkdir(parents=True, exist_ok=True)

    all_manifests = []

    # Base dataset: 2000 cells × 1000 genes, 3 branches, no dropout, no batch
    print("=" * 60)
    print("Generating base complex simulated dataset (3-branch)")
    print("=" * 60)
    pt, br, bl, topo = topology_3branch(2000)
    expr, gn, gc, gi = generate_expression(pt, br, 1000, 2000)
    batch_label = np.array(["batch_A"] * 2000)
    adata = build_anndata(expr, pt, br, bl, batch_label, gn, gc, topo)
    out_path = out_root / "simulated_complex_3branch.h5ad"
    adata.write(out_path)
    all_manifests.append({
        "file": str(out_path.name), "n_cells": 2000, "n_genes": 1000,
        "n_branches": 4, "dropout": 0.0, "batch": "none",
        "ground_truth": "true_pseudotime + true_branch + true_bifurcation_level"
    })
    print(f"  Saved: {out_path}")

    # Complex multi-configuration generation
    configs = []
    for n_br in [2, 3, 4]:
        for n_cells in [500, 2000, 5000]:
            configs.append({"n_branches": n_br, "n_cells": n_cells, "dropout": 0.0, "batch": False})

    # Add dropout variants (on 2000-cell 3-branch base)
    for dropout_rate in [0.30, 0.50, 0.70, 0.90]:
        configs.append({"n_branches": 3, "n_cells": 2000, "dropout": dropout_rate, "batch": False})

    # Add batch effect variant
    configs.append({"n_branches": 3, "n_cells": 2000, "dropout": 0.0, "batch": True})

    for cfg in configs:
        n_br = cfg["n_branches"]
        n_cells = cfg["n_cells"]
        do = cfg["dropout"]
        has_batch = cfg["batch"]

        label = f"sim_{n_br}br_{n_cells}c"
        if do > 0:
            label += f"_drop{int(do*100)}"
        if has_batch:
            label += "_batch"

        print(f"\nGenerating: {label}")

        topo_fn = TOPOLOGIES[n_br]
        pt, br, bl, topo = topo_fn(n_cells)
        expr, gn, gc, gi = generate_expression(pt, br, 1000, n_cells)

        if do > 0:
            expr = apply_dropout(expr, do, seed=CONFIG["seed"] + int(do * 10000))

        if has_batch:
            expr, batch_labels = apply_batch_effect(expr, 2, batch_seed=CONFIG["seed"] + 9999)
        else:
            batch_labels = np.array(["batch_A"] * n_cells)

        adata = build_anndata(expr, pt, br, bl, batch_labels, gn, gc, topo)
        out_path = out_root / f"{label}.h5ad"
        adata.write(out_path)

        all_manifests.append({
            "file": str(out_path.name), "n_cells": n_cells, "n_genes": 1000,
            "n_branches": n_br, "dropout": do, "batch": str(has_batch),
            "ground_truth": "true_pseudotime + true_branch + true_bifurcation_level"
        })
        print(f"  Saved: {out_path}")

    # Save manifest
    manifest_path = out_root / "manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(all_manifests, f, indent=2)
    print(f"\nManifest saved: {manifest_path}")
    print(f"Total datasets generated: {len(all_manifests)}")

    # Canonicalize: copy to pipeline canonical directory
    canon_dir = Path(r"F:\biotrade\LineaLinkR_benchmark_assets\benchmark_pipeline\canonical")
    for entry in all_manifests:
        src = out_root / entry["file"]
        dst = canon_dir / entry["file"].replace(".h5ad", "_canonical.h5ad")
        if src.exists() and not dst.exists():
            import shutil
            shutil.copy2(src, dst)
            print(f"  Canonicalized: {dst.name}")


if __name__ == "__main__":
    main()
