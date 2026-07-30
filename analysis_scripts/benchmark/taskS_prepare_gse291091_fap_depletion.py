from pathlib import Path
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.neighbors import KNeighborsClassifier


CONFIG = {
    "paths": {
        "input": Path(
            r"F:\biotrade\LineaLinkR_benchmark_assets\01_datasets\external_validation\GSE291091_fap_depletion"
        ),
        "reference": Path(
            r"F:\biotrade\LineaLinkR_benchmark_assets\01_datasets\canonical\GSE143437_muscle_regeneration_canonical.h5ad"
        ),
        "output": Path(
            r"F:\biotrade\0.64LLR\cache\biology\GSE291091_fap_depletion_canonical.h5ad"
        ),
        "qc_table": Path(
            r"F:\biotrade\0.64LLR\results\biology\muscle\gse291091_qc_and_annotation_summary.csv"
        ),
    },
    "threads": os.cpu_count(),
    "seed": 12345,
    "plot": {
        "font.family": "Arial",
        "font.size": 8,
        "axes.linewidth": 0.7,
        "figure.dpi": 120,
        "savefig.dpi": 400,
    },
    "qc": {"min_genes": 200, "max_genes": 2500, "max_pct_mt": 20.0},
    "analysis": {
        "target_sum": 10000,
        "n_top_genes": 3000,
        "n_pcs": 40,
        "n_neighbors": 25,
        "min_label_confidence": 0.50,
    },
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
os.environ["MKL_NUM_THREADS"] = str(CONFIG["threads"])
os.environ["OPENBLAS_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
plt.rcParams.update(CONFIG["plot"])

input_dir = CONFIG["paths"]["input"]
reference_path = CONFIG["paths"]["reference"]
output_path = CONFIG["paths"]["output"]
qc_path = CONFIG["paths"]["qc_table"]
output_path.parent.mkdir(parents=True, exist_ok=True)
qc_path.parent.mkdir(parents=True, exist_ok=True)

ct_1u = sc.read_10x_mtx(input_dir, prefix="GSM8829039_CT_1U_", var_names="gene_symbols")
ct_2u = sc.read_10x_mtx(input_dir, prefix="GSM8829040_CT_2U_", var_names="gene_symbols")
ab_1u = sc.read_10x_mtx(input_dir, prefix="GSM8829041_AB_1U_", var_names="gene_symbols")
ab_2u = sc.read_10x_mtx(input_dir, prefix="GSM8829042_AB_2U_", var_names="gene_symbols")
ct_1inj = sc.read_10x_mtx(input_dir, prefix="GSM8829043_CT_11NJ_", var_names="gene_symbols")
ct_2inj = sc.read_10x_mtx(input_dir, prefix="GSM8829044_CT_2INJ_", var_names="gene_symbols")
ab_1inj = sc.read_10x_mtx(input_dir, prefix="GSM8829045_AB_1INJ_", var_names="gene_symbols")
ab_2inj = sc.read_10x_mtx(input_dir, prefix="GSM8829046_AB_2INJ_", var_names="gene_symbols")

query = sc.concat(
    {
        "CT_1U": ct_1u,
        "CT_2U": ct_2u,
        "AB_1U": ab_1u,
        "AB_2U": ab_2u,
        "CT_1INJ": ct_1inj,
        "CT_2INJ": ct_2inj,
        "AB_1INJ": ab_1inj,
        "AB_2INJ": ab_2inj,
    },
    label="sample_id",
    index_unique="-",
    join="inner",
)
query.obs["fap_status"] = np.where(
    query.obs["sample_id"].str.startswith("AB"), "depleted", "control"
)
query.obs["injury_status"] = np.where(
    query.obs["sample_id"].str.endswith("INJ"), "3dpi", "uninjured"
)
query.obs["replicate"] = query.obs["sample_id"].str.extract(r"_(\d)", expand=False)
query.var_names_make_unique()
query.var["mt"] = query.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(query, qc_vars=["mt"], inplace=True)

qc_before = (
    query.obs.groupby(["sample_id", "fap_status", "injury_status"], observed=True)
    .agg(
        n_cells=("sample_id", "size"),
        median_counts=("total_counts", "median"),
        median_genes=("n_genes_by_counts", "median"),
        median_pct_mt=("pct_counts_mt", "median"),
    )
    .reset_index()
    .assign(stage="before_qc")
)

keep = (
    (query.obs["n_genes_by_counts"] >= CONFIG["qc"]["min_genes"])
    & (query.obs["n_genes_by_counts"] <= CONFIG["qc"]["max_genes"])
    & (query.obs["pct_counts_mt"] <= CONFIG["qc"]["max_pct_mt"])
)
query = query[keep].copy()
sc.pp.filter_genes(query, min_cells=4)
query.layers["counts"] = query.X.copy()
sc.pp.normalize_total(query, target_sum=CONFIG["analysis"]["target_sum"])
sc.pp.log1p(query)

qc_after = (
    query.obs.groupby(["sample_id", "fap_status", "injury_status"], observed=True)
    .agg(
        n_cells=("sample_id", "size"),
        median_counts=("total_counts", "median"),
        median_genes=("n_genes_by_counts", "median"),
        median_pct_mt=("pct_counts_mt", "median"),
    )
    .reset_index()
    .assign(stage="after_qc")
)

reference = sc.read_h5ad(reference_path)
common_genes = reference.var_names.intersection(query.var_names)
reference = reference[:, common_genes].copy()
query_map = query[:, common_genes].copy()
sc.pp.normalize_total(reference, target_sum=CONFIG["analysis"]["target_sum"])
sc.pp.log1p(reference)
sc.pp.highly_variable_genes(
    reference,
    n_top_genes=CONFIG["analysis"]["n_top_genes"],
    flavor="seurat",
)
reference = reference[:, reference.var["highly_variable"]].copy()
query_map = query_map[:, reference.var_names].copy()
sc.tl.pca(
    reference,
    n_comps=CONFIG["analysis"]["n_pcs"],
    svd_solver="arpack",
    random_state=CONFIG["seed"],
)
sc.pp.neighbors(
    reference,
    n_neighbors=CONFIG["analysis"]["n_neighbors"],
    n_pcs=CONFIG["analysis"]["n_pcs"],
    random_state=CONFIG["seed"],
)
sc.tl.ingest(query_map, reference, obs="cell_type", embedding_method="pca")

classifier = KNeighborsClassifier(
    n_neighbors=CONFIG["analysis"]["n_neighbors"],
    weights="distance",
    n_jobs=CONFIG["threads"],
)
classifier.fit(reference.obsm["X_pca"], reference.obs["cell_type"].astype(str))
label_probabilities = classifier.predict_proba(query_map.obsm["X_pca"])
query.obs["predicted_cell_type"] = classifier.classes_[label_probabilities.argmax(axis=1)]
query.obs["label_confidence"] = label_probabilities.max(axis=1)
query.obsm["X_pca_transfer"] = query_map.obsm["X_pca"]
sc.pp.highly_variable_genes(
    query,
    n_top_genes=CONFIG["analysis"]["n_top_genes"],
    flavor="seurat",
)
sc.tl.pca(
    query,
    n_comps=CONFIG["analysis"]["n_pcs"],
    mask_var="highly_variable",
    svd_solver="arpack",
    random_state=CONFIG["seed"],
)
sc.pp.neighbors(
    query,
    n_neighbors=CONFIG["analysis"]["n_neighbors"],
    n_pcs=CONFIG["analysis"]["n_pcs"],
    random_state=CONFIG["seed"],
)
sc.tl.umap(query, random_state=CONFIG["seed"])
sc.tl.leiden(
    query,
    resolution=0.6,
    random_state=CONFIG["seed"],
    key_added="leiden",
    flavor="igraph",
    n_iterations=2,
    directed=False,
)
query.obs["cell_type"] = np.select(
    [
        query.obs["leiden"].eq("0"),
        query.obs["leiden"].isin(["1", "2", "12"]),
        query.obs["leiden"].isin(["3", "4", "5", "6", "14"]),
        query.obs["leiden"].eq("7"),
        query.obs["leiden"].eq("8"),
        query.obs["leiden"].eq("10"),
        query.obs["leiden"].eq("11"),
        query.obs["leiden"].eq("13"),
    ],
    [
        "Smooth muscle cells",
        "Endothelial",
        "Mononuclear phagocytes",
        "Quiescent MuSCs",
        "FAPs",
        "Lymphoid/APCs",
        "FAP-like/Tenocytes",
        "Neutrophils",
    ],
    default="unassigned",
)
query.obs["cluster"] = query.obs["cell_type"].astype(str)
query.obs["clusters"] = query.obs["cell_type"].astype(str)

annotation_summary = (
    query.obs.groupby(
        ["sample_id", "fap_status", "injury_status", "cell_type"], observed=True
    )
    .agg(
        n_cells=("cell_type", "size"),
        median_label_confidence=("label_confidence", "median"),
    )
    .reset_index()
    .assign(stage="annotation")
)
pd.concat([qc_before, qc_after, annotation_summary], ignore_index=True).to_csv(
    qc_path, index=False, encoding="utf-8", lineterminator="\n"
)
query.uns["source_accession"] = "GSE291091"
query.uns["source_doi"] = "10.1002/advs.202507444"
query.uns["annotation_method"] = (
    "Consensus annotation combined de novo HVG/PCA/Leiden clusters, lineage marker "
    "genes reported by Liu et al. (2026), and kNN transfer from the author-annotated "
    "GSE143437 atlas. "
    "The mixed activated-myogenic/neural cluster was labeled unassigned."
)
query.write_h5ad(output_path, compression="gzip")

print(query)
print(annotation_summary.to_string(index=False))
