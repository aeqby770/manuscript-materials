from pathlib import Path
import os

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse


CONFIG = {
    "paths": {
        "input": Path(
            r"F:\biotrade\0.64LLR\cache\biology\GSE291091_fap_depletion_canonical.h5ad"
        ),
        "output": Path(
            r"F:\biotrade\0.64LLR\results\biology\muscle"
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
    "contexts": {
        "Mononuclear phagocytes": ["Mononuclear phagocytes"],
        "FAP compartment": ["FAPs", "FAP-like/Tenocytes"],
    },
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
os.environ["MKL_NUM_THREADS"] = str(CONFIG["threads"])
os.environ["OPENBLAS_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
plt.rcParams.update(CONFIG["plot"])

input_path = CONFIG["paths"]["input"]
output_dir = CONFIG["paths"]["output"]
output_dir.mkdir(parents=True, exist_ok=True)

adata = ad.read_h5ad(input_path)
obs = adata.obs.copy()
selected = obs["cell_type"].isin(
    CONFIG["contexts"]["Mononuclear phagocytes"]
    + CONFIG["contexts"]["FAP compartment"]
)
selected_obs = obs.loc[selected].copy()
selected_obs["context"] = np.where(
    selected_obs["cell_type"].eq("Mononuclear phagocytes"),
    "Mononuclear phagocytes",
    "FAP compartment",
)
selected_obs["group_id"] = (
    selected_obs["context"].astype(str) + "|" + selected_obs["sample_id"].astype(str)
)
group_factor = pd.Categorical(selected_obs["group_id"], ordered=True)
counts = adata.layers["counts"][selected.to_numpy(), :].tocsr()
group_matrix = sparse.csr_matrix(
    (
        np.ones(selected_obs.shape[0], dtype=np.float64),
        (group_factor.codes, np.arange(selected_obs.shape[0])),
    ),
    shape=(len(group_factor.categories), selected_obs.shape[0]),
)
pseudobulk = group_matrix @ counts
pseudobulk_table = pd.DataFrame.sparse.from_spmatrix(
    pseudobulk.T,
    index=adata.var_names.astype(str),
    columns=group_factor.categories.astype(str),
)
pseudobulk_table.index.name = "gene"

group_metadata = pd.DataFrame({"group_id": group_factor.categories.astype(str)})
group_metadata[["context", "sample_id"]] = group_metadata["group_id"].str.split(
    "|", regex=False, expand=True
)
sample_metadata = obs[
    ["sample_id", "fap_status", "injury_status", "replicate"]
].drop_duplicates()
group_metadata = group_metadata.merge(sample_metadata, on="sample_id", how="left")
group_metadata["n_cells"] = selected_obs.groupby("group_id", observed=True).size().reindex(
    group_metadata["group_id"]
).to_numpy()

cell_composition = (
    obs.groupby(
        ["sample_id", "fap_status", "injury_status", "replicate", "cell_type"],
        observed=True,
    )
    .size()
    .rename("n_cells")
    .reset_index()
)
cell_composition["sample_total"] = cell_composition.groupby(
    "sample_id", observed=True
)["n_cells"].transform("sum")
cell_composition["captured_fraction"] = (
    cell_composition["n_cells"] / cell_composition["sample_total"]
)

macrophage_obs = obs.loc[
    obs["cell_type"].eq("Mononuclear phagocytes")
    & obs["injury_status"].eq("3dpi"),
    [
        "sample_id",
        "fap_status",
        "replicate",
        "predicted_cell_type",
        "label_confidence",
    ],
].copy()
macrophage_obs["subtype_audit_label"] = np.where(
    macrophage_obs["label_confidence"].ge(0.50),
    macrophage_obs["predicted_cell_type"].astype(str),
    "Low-confidence transfer",
)
macrophage_subtype_composition = (
    macrophage_obs.groupby(
        ["sample_id", "fap_status", "replicate", "subtype_audit_label"],
        observed=True,
    )
    .agg(
        n_cells=("subtype_audit_label", "size"),
        median_label_confidence=("label_confidence", "median"),
    )
    .reset_index()
)
macrophage_subtype_composition["sample_total"] = (
    macrophage_subtype_composition.groupby("sample_id", observed=True)[
        "n_cells"
    ].transform("sum")
)
macrophage_subtype_composition["within_phagocyte_fraction"] = (
    macrophage_subtype_composition["n_cells"]
    / macrophage_subtype_composition["sample_total"]
)

pseudobulk_table.to_csv(
    output_dir / "gse291091_pseudobulk_counts.csv",
    encoding="utf-8",
    lineterminator="\n",
)
group_metadata.to_csv(
    output_dir / "gse291091_pseudobulk_metadata.csv",
    index=False,
    encoding="utf-8",
    lineterminator="\n",
)
cell_composition.to_csv(
    output_dir / "gse291091_cell_composition.csv",
    index=False,
    encoding="utf-8",
    lineterminator="\n",
)
macrophage_subtype_composition.to_csv(
    output_dir / "gse291091_macrophage_subtype_composition.csv",
    index=False,
    encoding="utf-8",
    lineterminator="\n",
)

print(group_metadata.to_string(index=False))
