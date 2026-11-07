#!/usr/bin/env python
# coding: utf-8
"""Run the TICCI CellChat, CCI-weighted trajectory, and fitted-time workflow."""
import os
import shutil
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import inspect
import igraph as ig

# Keep zero-edge graph coordinates two-dimensional for current SciPy.
native_sparse_name=ig.Graph.get_adjacency_sparse.__name__
native_sparse_source=inspect.getsource(ig.Graph.get_adjacency_sparse)
assert native_sparse_source.count("list(zip(*edges))")==1
native_sparse_source=native_sparse_source.replace(
    "list(zip(*edges))","tuple(np.asarray(edges,dtype=np.int64).reshape(-1,2).T)")
native_sparse_scope=dict(ig.Graph.get_adjacency_sparse.__globals__,np=np)
exec(native_sparse_source,native_sparse_scope)
ig.Graph.get_adjacency_sparse=native_sparse_scope[native_sparse_name]

import scanpy._utils as sc_utils
native_igraph_source=inspect.getsource(sc_utils.get_igraph_from_adjacency)
assert native_igraph_source.count("weights = adjacency[sources, targets]")==1
native_igraph_source=native_igraph_source.replace(
    "weights = adjacency[sources, targets]",
    "weights = sparse.csr_matrix(adjacency[sources, targets]).toarray().ravel()")
native_igraph_scope=dict(sc_utils.get_igraph_from_adjacency.__globals__,sparse=sparse)
exec("from __future__ import annotations\n"+native_igraph_source,native_igraph_scope)
sc_utils.get_igraph_from_adjacency=native_igraph_scope["get_igraph_from_adjacency"]


matplotlib.use("Agg")

CONFIG = {
    "paths": {
        "vendor": Path(os.environ["LLR_TICCI_ROOT"]).resolve(),
    },
    "threads": int(os.environ.get("LLR_NATIVE_THREADS", os.cpu_count() or 1)),
    "seed": 12345,
    "plot": {},
    "cci_weight": 1.0,
    "expression_threshold": 1.0,
    "n_top_genes": 2000,
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])

inpath = Path(sys.argv[1]).resolve()
outdir = Path(sys.argv[2]).resolve()
outdir.mkdir(parents=True, exist_ok=True)
dataset = outdir.name
cellchat_csv = outdir / "ticci_cellchat_interactions.csv"

adata = sc.read_h5ad(inpath)
adata.X = sparse.csr_matrix(adata.X, dtype=np.float64)
sc.pp.normalize_total(adata, target_sum=10000)
sc.pp.log1p(adata)
sc.pp.filter_genes(adata, min_cells=3)

saved_entropy=pd.read_csv(outdir/"ticci_functional_entropy.csv",dtype={"cell_id":str}).set_index("cell_id")
assert set(saved_entropy.index)==set(adata.obs_names)
entropy=saved_entropy.reindex(adata.obs_names)["entropy"].to_numpy(dtype=float)
assert np.isfinite(entropy).all()

mean = np.asarray(adata.X.mean(axis=0)).ravel()
sq_mean = np.asarray(adata.X.multiply(adata.X).mean(axis=0)).ravel()
variance = np.nan_to_num(sq_mean - np.square(mean), nan=0.0, posinf=0.0, neginf=0.0)
n_top = min(CONFIG["n_top_genes"], adata.n_vars)
top_index = np.argsort(variance)[::-1][:n_top]
work = adata[:, top_index].copy()
work.obs["entropy"] = entropy

sc.pp.pca(work, n_comps=min(50, work.n_vars - 1, work.n_obs - 1))
sc.pp.neighbors(work, n_neighbors=10, n_pcs=min(50, work.obsm["X_pca"].shape[1]))

saved_groups = pd.read_csv(outdir / "ticci_cellchat_groups.csv", dtype=str).set_index("cell_id")
assert saved_groups.index.is_unique and set(saved_groups.index) == set(work.obs_names)
groups = saved_groups.reindex(work.obs_names)["Classification"].to_numpy()
assert pd.notna(groups).all()

gene_binding = pd.concat(
    [pd.DataFrame({"original_gene_id": adata.var_names, "resource_symbol": adata.var_names}),
     *map(pd.read_csv, outdir.glob("ticci_gene_id_binding.csv"))],
    ignore_index=True,
).drop_duplicates("original_gene_id", keep="last").set_index("original_gene_id")
resource_symbols = gene_binding.reindex(adata.var_names)["resource_symbol"]
assert resource_symbols.notna().all() and resource_symbols.is_unique
adata.var["resource_symbol"] = resource_symbols.to_numpy()

net = pd.read_csv(cellchat_csv, encoding="utf-8-sig", dtype={"source": str, "target": str})
parts = (
    net["interaction_name_2"]
    .astype(str)
    .str.replace(r"[()]", "", regex=True)
    .str.split(r"\s+-\s+", n=1, expand=True)
    .reindex(columns=[0, 1]).astype("string")
)
net = net.assign(
    cci_id=np.arange(len(net), dtype=int),
    ligand_complex=parts[0],
    receptor_complex=parts[1],
)

ligand_requirements = net[["cci_id", "ligand_complex"]].rename(
    columns={"ligand_complex": "gene"}
)
ligand_requirements["side"] = "ligand"
ligand_requirements["gene"] = ligand_requirements["gene"].str.split("+")
ligand_requirements = ligand_requirements.explode("gene", ignore_index=True)

receptor_requirements = net[["cci_id", "receptor_complex"]].rename(
    columns={"receptor_complex": "gene"}
)
receptor_requirements["side"] = "receptor"
receptor_requirements["gene"] = receptor_requirements["gene"].str.split("+")
receptor_requirements = receptor_requirements.explode("gene", ignore_index=True)

requirements = pd.concat(
    [ligand_requirements, receptor_requirements], ignore_index=True
)
requirements["gene"] = requirements["gene"].str.strip()
requirements["present"] = requirements["gene"].isin(resource_symbols)
required_counts = requirements.groupby(["cci_id", "side"]).size().unstack(fill_value=0)
present_counts = (
    requirements.groupby(["cci_id", "side"])["present"].sum().unstack(fill_value=0)
)
valid_cci = required_counts.index[(required_counts == present_counts).all(axis=1)]
net = net[net["cci_id"].isin(valid_cci)].copy()
requirements = requirements[requirements["cci_id"].isin(valid_cci)].copy()

connectivity = work.obsp["connectivities"].tocoo()
base_edges = pd.DataFrame(
    {"row": connectivity.row, "col": connectivity.col, "weight": connectivity.data}
)
base_edges = base_edges[base_edges["row"] < base_edges["col"]].reset_index(drop=True)
base_edges["edge_id"] = np.arange(len(base_edges), dtype=int)
forward = base_edges[["edge_id", "row", "col"]].copy()
reverse = base_edges[["edge_id", "col", "row"]].rename(
    columns={"col": "row", "row": "col"}
)
directed_edges = pd.concat([forward, reverse], ignore_index=True)
directed_edges["source"] = groups[directed_edges["row"].to_numpy()]
directed_edges["target"] = groups[directed_edges["col"].to_numpy()]

candidates = directed_edges.merge(
    net[["cci_id", "source", "target", "prob"]],
    on=["source", "target"],
    how="inner",
)
candidates["candidate_id"] = np.arange(len(candidates), dtype=int)
candidate_requirements = candidates[["candidate_id", "edge_id", "row", "col", "cci_id"]].merge(
    requirements[["cci_id", "side", "gene"]], on="cci_id", how="inner"
)
gene_index = pd.Series(np.arange(adata.n_vars, dtype=int), index=resource_symbols.to_numpy())
cell_index = np.where(
    candidate_requirements["side"].eq("ligand"),
    candidate_requirements["row"],
    candidate_requirements["col"],
).astype(int)
required_gene_index = gene_index.loc[candidate_requirements["gene"]].to_numpy(dtype=int)
expression_values = sparse.csr_matrix(adata.X[cell_index, required_gene_index]).toarray().ravel()
candidate_requirements["expressed"] = expression_values > CONFIG["expression_threshold"]
valid_expression = (
    candidate_requirements.groupby(["candidate_id", "cci_id"], sort=False)["expressed"]
    .all()
    .reset_index()
)
valid_expression = valid_expression[valid_expression["expressed"]]
valid_candidates = candidates.merge(
    valid_expression[["candidate_id", "cci_id"]],
    on=["candidate_id", "cci_id"],
    how="inner",
)

edge_delta = valid_candidates.groupby("edge_id", sort=False)["prob"].sum()
weighted_edges = base_edges.loc[edge_delta.index]
delta_values = edge_delta.to_numpy(dtype=float) * CONFIG["cci_weight"]
delta = sparse.coo_matrix(
    (
        np.concatenate([delta_values, delta_values]),
        (
            np.concatenate([weighted_edges["row"], weighted_edges["col"]]),
            np.concatenate([weighted_edges["col"], weighted_edges["row"]]),
        ),
    ),
    shape=work.obsp["connectivities"].shape,
).tocsr()
work.obsp["connectivities"] = (work.obsp["connectivities"] + delta).tocsr()

sc.tl.louvain(
    work,
    use_weights=True,
    flavor="igraph",
    resolution=1.0,
    random_state=CONFIG["seed"],
    key_added="louvain",
)
work.obs["louvain"] = work.obs["louvain"].astype(int).astype("category")
sc.tl.paga(work, groups="louvain")
sc.settings.figdir = outdir
sc.pl.paga(work, plot=False, show=False)

sys.path.insert(0, str(CONFIG["paths"]["vendor"].parent))
sys.path.insert(0, str(CONFIG["paths"]["vendor"]))
from Python_tools.ticci import ticci
import matplotlib.pyplot as plt

original_colorbar = plt.colorbar
plt.colorbar = lambda mappable, **kwargs: original_colorbar(
    mappable, ax=plt.gca(), **kwargs
)
lineage_path = outdir / "ticci_lineage.svg"
ticci.getLineage(work, "_ticci_stable_entropy", str(lineage_path))
plt.colorbar = original_colorbar

root_cluster = int(np.argmax(work.uns["louvain_entropy"]))
root_mask = work.obs["louvain"].astype(int).to_numpy() == root_cluster
root_index = int(np.argmax(np.where(root_mask, entropy, -np.inf)))
work.uns["iroot"] = root_index
sc.tl.diffmap(work)
sc.tl.dpt(work)

tree = work.uns["paga"]["directed_tree"]
tree_rows, tree_cols = tree.nonzero()
tree_weights = sparse.csr_matrix(tree[tree_rows, tree_cols]).toarray().ravel()
pd.DataFrame(
    {"from_cluster": tree_rows, "to_cluster": tree_cols, "weight": tree_weights}
).to_csv(outdir / "ticci_lineage_edges.csv", index=False)

pd.DataFrame(
    {
        "cell_id": work.obs_names,
        "pseudotime": work.obs["dpt_pseudotime"].to_numpy(dtype=float),
        "cluster": work.obs["louvain"].astype(int).to_numpy(),
        "classification": groups,
        "entropy": entropy,
    }
).to_csv(outdir / "ticci_pseudotime.csv", index=False)

summary = pd.DataFrame(
    [
        {
            "n_cells": int(work.n_obs),
            "n_clusters": int(work.obs["louvain"].nunique()),
            "n_paga_edges": int(work.uns["paga"]["connectivities"].nnz),
            "n_lineage_tree_edges": int(len(tree_rows)),
            "n_cellchat_interactions": int(len(net)),
            "n_cci_weighted_edges": int(len(edge_delta)),
            "cci_weight": CONFIG["cci_weight"],
            "root_cluster": root_cluster,
            "method": "ticci_SLICE_functional_entropy_fixed_k",
            "bootstrap_B":100,"bootstrap_genes":1000,"seed":12345,
            "n_finite_pseudotime":int(np.isfinite(work.obs["dpt_pseudotime"]).sum()),
        }
    ]
)
summary.to_csv(outdir / "ticci_summary.csv", index=False)
work.write_h5ad(outdir/"ticci_author_result.h5ad",compression="gzip")

stable_entropy_source = outdir / "paga_ticci_stable_entropy.pdf"
stable_entropy_target = outdir / "ticci_stable_entropy.pdf"
shutil.move(stable_entropy_source, stable_entropy_target)

print(f"OK (ticci_SLICE_functional_entropy_fixed_k): {outdir / 'ticci_summary.csv'}")
