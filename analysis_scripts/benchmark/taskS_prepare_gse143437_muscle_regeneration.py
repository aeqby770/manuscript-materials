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
            r"F:\biotrade\LineaLinkR_benchmark_assets\01_datasets\external_validation\GSE232106_muscle_aging"
        ),
        "output": Path(
            r"F:\biotrade\LineaLinkR_benchmark_assets\01_datasets\canonical\GSE143437_muscle_regeneration_canonical.h5ad"
        ),
        "qc_table": Path(
            r"F:\biotrade\0.64LLR\results\biology\muscle\gse143437_input_summary.csv"
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
    "time_days": {"Day 0": 0.0, "Day 2": 2.0, "Day 5": 5.0, "Day 7": 7.0},
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
os.environ["MKL_NUM_THREADS"] = str(CONFIG["threads"])
os.environ["OPENBLAS_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
plt.rcParams.update(CONFIG["plot"])

input_dir = CONFIG["paths"]["input"]
output_path = CONFIG["paths"]["output"]
qc_path = CONFIG["paths"]["qc_table"]
output_path.parent.mkdir(parents=True, exist_ok=True)
qc_path.parent.mkdir(parents=True, exist_ok=True)

metadata = pd.read_csv(
    input_dir / "GSE143437_DeMicheli_MuSCatlas_metadata.txt.gz",
    sep="\t",
    index_col=0,
    compression="gzip",
    encoding="utf-8-sig",
)
raw = pd.read_csv(
    input_dir / "GSE143437_DeMicheli_MuSCatlas_rawdata.txt.gz",
    sep="\t",
    index_col=0,
    compression="gzip",
    encoding="utf-8-sig",
).astype(np.uint32, copy=False)

metadata = metadata.loc[raw.columns].copy()
metadata["time_days"] = metadata["injury"].map(CONFIG["time_days"]).astype(float)
metadata["sample_id"] = metadata["sampleID"].astype(str)
metadata["cell_type"] = metadata["cell_annotation"].astype(str)
metadata["cluster"] = metadata["cell_annotation"].astype(str)
metadata["clusters"] = metadata["cell_annotation"].astype(str)

counts = sparse.csr_matrix(raw.to_numpy(dtype=np.uint32, copy=False).T).astype(
    np.int32, copy=False
)
genes = pd.DataFrame(index=raw.index.astype(str))
adata = ad.AnnData(X=counts, obs=metadata, var=genes)
adata.uns["source_accession"] = "GSE143437"
adata.uns["source_doi"] = "10.1016/j.celrep.2020.02.067"
adata.uns["validation_design"] = (
    "Blind LineaLinkR inference: time_days and injury are retained only for post hoc "
    "validation and are not used for root, graph, lineage, or ligand-receptor scoring."
)

input_summary = (
    metadata.groupby(["injury", "time_days", "sample_id", "cell_type"], observed=True)
    .size()
    .rename("n_cells")
    .reset_index()
    .sort_values(["time_days", "sample_id", "cell_type"])
)
input_summary.to_csv(qc_path, index=False, encoding="utf-8", lineterminator="\n")
adata.write_h5ad(output_path, compression="gzip")

print(adata)
print(input_summary.groupby(["injury", "cell_type"], observed=True)["n_cells"].sum())
