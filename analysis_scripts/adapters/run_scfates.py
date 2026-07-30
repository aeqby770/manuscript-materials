"""scFates adapter: principal tree + pseudotime + lineage.

Usage: python run_scfates.py <canonical_h5ad> <output_csv>
Outputs: <output_csv> (pseudotime), <output_dir>/scfates_lineage_*.csv (tree + segments)
"""
import sys
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scFates
import _mem_guard
from _benchmark_contract import benchmark_track, expression_root_position, root_principal_point
_mem_guard.cap_self(55)

CONFIG = {
    "paths": {},
    "threads": os.cpu_count() or 1,
    "seed": 12345,
    "plot": {},
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])

def run_scfates(h5ad_path, output_csv, nodes=50):
    adata = sc.read(h5ad_path)
    canonical_cell_ids = adata.obs_names.astype(str).to_numpy(copy=True)
    if canonical_cell_ids.size != adata.n_obs or pd.Index(canonical_cell_ids).has_duplicates:
        raise RuntimeError("scFates: canonical obs_names must be complete and unique")
    
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    if hasattr(adata.X, "toarray"):
        adata.X = adata.X.toarray()
    adata.X = np.nan_to_num(adata.X, nan=0, posinf=0, neginf=0)
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata_hvg = adata[:, adata.var.highly_variable].copy()
    sc.pp.pca(adata_hvg, n_comps=30)
    sc.pp.neighbors(adata_hvg)
    sc.tl.umap(adata_hvg)
    
    scFates.tl.tree(adata_hvg, Nodes=nodes, use_rep='X_pca',
                     method='ppt', ppt_sigma=0.1, ppt_lambda=1,
                     ppt_nsteps=50, seed=12345)
    scFates.tl.cleanup(adata_hvg)
    scFates.tl.merge_empty_segments(adata_hvg)
    graph = adata_hvg.uns["graph"]
    tips = np.asarray(graph.get("tips", []), dtype=int)
    if benchmark_track() == "common_root":
        root_pp = root_principal_point(adata_hvg)
    elif benchmark_track() == "annotation_free":
        root_pp = root_principal_point(
            adata_hvg,
            root_position=expression_root_position(adata_hvg),
        )
    else:
        root_pp = int(tips[0] if tips.size else 0)
    scFates.tl.root(adata_hvg, root_pp)
    scFates.tl.pseudotime(adata_hvg, n_map=1, seed=12345)
    
    if not np.array_equal(adata_hvg.obs_names.astype(str).to_numpy(), canonical_cell_ids):
        raise RuntimeError("scFates: preprocessing changed canonical cell order")
    pseudotime = pd.to_numeric(
        adata_hvg.obs.get('t', pd.Series(np.nan, index=adata_hvg.obs_names)),
        errors="coerce",
    ).to_numpy(dtype=float)
    if pseudotime.size != canonical_cell_ids.size:
        raise RuntimeError(
            f"scFates: pseudotime length {pseudotime.size} != canonical cells {canonical_cell_ids.size}"
        )
    if not np.isfinite(pseudotime).all() or np.ptp(pseudotime) <= 0:
        raise RuntimeError("scFates: pseudotime must be finite for every cell and have non-zero range")
    pseudotime = (pseudotime - pseudotime.min()) / np.ptp(pseudotime)
    
    out = pd.DataFrame({"cell_id": canonical_cell_ids, "pseudotime": pseudotime})
    out.to_csv(output_csv, index=False)
    written = pd.read_csv(output_csv, dtype={"cell_id": str})
    written_pt = pd.to_numeric(written["pseudotime"], errors="coerce").to_numpy(dtype=float)
    if (
        written.columns.tolist() != ["cell_id", "pseudotime"]
        or not np.array_equal(written["cell_id"].astype(str).to_numpy(), canonical_cell_ids)
        or written["cell_id"].astype(str).duplicated().any()
        or not np.isfinite(written_pt).all()
        or np.ptp(written_pt) <= 0
    ):
        raise RuntimeError("scFates: written CSV failed the full-cell pseudotime contract")
    
    # Export tree topology + lineage info (dataset-specific filenames)
    out_dir = os.path.dirname(output_csv)
    ds_name = os.path.basename(output_csv).replace('scfates_', '').replace('_pseudotime.csv', '')
    graph = adata_hvg.uns.get('graph', {})
    ppt = adata_hvg.uns.get('ppt', {})
    
    # Save adjacency matrix as edge list
    if 'B' in graph:
        B = graph['B']
        n_nodes = B.shape[0]
        edges = []
        for i in range(n_nodes):
            for j in range(i + 1, n_nodes):
                if B[i, j] > 0:
                    edges.append({'from_node': i, 'to_node': j})
        if edges:
            edges_df = pd.DataFrame(edges)
            edges_df.to_csv(os.path.join(out_dir, f'scfates_lineage_edges_{ds_name}.csv'), index=False)
    
    # Save segment info (pp_seg: segments connecting principal points)
    if 'pp_seg' in graph:
        pp_seg = np.array(graph['pp_seg'])
        seg_df = pd.DataFrame(pp_seg, columns=['from_pp', 'to_pp', 'length', 'segment_id'])
        seg_df.to_csv(os.path.join(out_dir, f'scfates_lineage_segments_{ds_name}.csv'), index=False)
    
    # Save per-cell segment assignments
    if 'seg' in adata_hvg.obs.columns:
        cell_seg = pd.DataFrame({
            'cell': adata_hvg.obs_names,
            'segment': adata_hvg.obs['seg'].values,
            'pseudotime': pseudotime
        })
        cell_seg.to_csv(os.path.join(out_dir, f'scfates_lineage_cell_segments_{ds_name}.csv'), index=False)
    
    tips = graph.get('tips', np.array([]))
    forks = graph.get('forks', np.array([]))
    n_tips = len(tips) if hasattr(tips, '__len__') else 0
    n_forks = len(forks) if hasattr(forks, '__len__') else 0
    n_pp = graph['B'].shape[0] if 'B' in graph else 0
    n_seg = graph['pp_seg'].shape[0] if 'pp_seg' in graph else 0
    
    print(f"scFates tree: {n_pp} PP, {n_seg} segments, {n_tips} tips, {n_forks} forks, root={graph.get('root','?')}, track={benchmark_track()}")
    
    n_valid = int(np.isfinite(pseudotime).sum())
    print(f"scFates OK: {adata_hvg.n_obs} cells, {n_valid} with pseudotime, nodes={nodes}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python run_scfates.py <canonical_h5ad> <output_csv>")
        sys.exit(1)
    run_scfates(sys.argv[1], sys.argv[2])
