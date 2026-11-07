import os
from pathlib import Path
CONFIG = {
    "paths": {"root": os.environ.get("LLR_WORK_ROOT", str(Path(__file__).resolve().parents[2]))},
    "threads": os.cpu_count(),
    "plot": {},
    "n_variable_genes": 2000,
    "n_pcs": 20,
    "n_clusters": 20,
    "min_features": 200,
    "max_features": 2500,
    "max_mito_fraction": 0.20,
    "scale_factor": 10000
}
os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
os.environ["OPENBLAS_NUM_THREADS"] = str(CONFIG["threads"])
os.environ["MKL_NUM_THREADS"] = str(CONFIG["threads"])
import gzip
import tarfile
from functools import partial
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread, mmwrite
from scipy.linalg import eigh
from scipy.cluster.hierarchy import linkage, fcluster

root = Path(CONFIG["paths"]["root"])
data = root / "inputs/comparison/biology/muscle"
out = root / "results/validation/gse291091"
out.mkdir(parents=True, exist_ok=True)
raw = data / "raw"
raw.mkdir(exist_ok=True)
archive = tarfile.open(data / "GSE291091_RAW.tar")
members = pd.Series(archive.getnames(), name="filename")
assert members.str.fullmatch(r"GSM\d+_.+_(barcodes.tsv.gz|features.tsv.gz|matrix.mtx.gz)").all()
archive.extractall(raw, filter="data")
files = members.to_frame()
files[["accession", "sample_file", "kind"]] = files.filename.str.extract(r"^(GSM\d+)_(.+)_(barcodes.tsv.gz|features.tsv.gz|matrix.mtx.gz)$")
files["sample_id"] = files.sample_file.str.replace("11NJ", "1INJ", regex=False)
samples = files.pivot(index=["accession", "sample_id"], columns="kind", values="filename").reset_index()
samples["condition"] = np.where(samples.sample_id.str.startswith("CT"), "control", "depleted")
samples["injury"] = np.where(samples.sample_id.str.endswith("INJ"), "3dpi", "uninjured")
samples.to_csv(out / "samples.csv", index=False, encoding="utf-8")
inj = samples.loc[samples.injury.eq("3dpi")].sort_values(["condition", "sample_id"])
matrix_paths = raw / inj["matrix.mtx.gz"]
matrices = list(map(mmread, map(gzip.open, matrix_paths)))
counts = sparse.hstack(matrices, format="csr")
features = pd.read_csv(raw / inj["features.tsv.gz"].iloc[0], sep="\t", header=None)
gene_symbols = pd.Index(features[1].astype(str))
feature_tables = list(map(partial(pd.read_csv, sep="\t", header=None), raw / inj["features.tsv.gz"]))
assert pd.concat(feature_tables, keys=inj.sample_id.tolist()).groupby(level=1).nunique().to_numpy().max() == 1
barcodes = pd.concat(map(partial(pd.read_csv, sep="\t", header=None), raw / inj["barcodes.tsv.gz"]), ignore_index=True)[0]
sizes = np.array(list(map(lambda x: x.shape[1], matrices)))
meta = pd.DataFrame({"barcode": barcodes, "sample_id": np.repeat(inj.sample_id.to_numpy(), sizes), "condition": np.repeat(inj.condition.to_numpy(), sizes)})
meta["cell_id"] = meta.sample_id + ":" + meta.barcode
lib = np.asarray(counts.sum(axis=0)).ravel()
nfeat = np.asarray((counts > 0).sum(axis=0)).ravel()
mito = np.asarray(counts[gene_symbols.str.startswith("mt-")].sum(axis=0)).ravel() / np.maximum(lib, 1)
meta["library_size"] = lib
meta["n_features"] = nfeat
meta["mito_fraction"] = mito
keep = (nfeat >= CONFIG["min_features"]) & (nfeat <= CONFIG["max_features"]) & (mito <= CONFIG["max_mito_fraction"])
meta["retained"] = keep
meta.to_csv(out / "all_cell_qc.csv.gz", index=False, encoding="utf-8")
counts = counts[:, keep]
meta = meta.loc[keep].reset_index(drop=True)
lib = lib[keep]
expr = counts @ sparse.diags(CONFIG["scale_factor"] / lib)
expr.data = np.log1p(expr.data)
mu = np.asarray(expr.mean(axis=1)).ravel()
variance = np.asarray(expr.power(2).mean(axis=1)).ravel() - mu * mu
detected = np.asarray((counts > 0).sum(axis=1)).ravel()
eligible = np.flatnonzero(detected > 3)
variable = eligible[np.argsort(variance[eligible], kind="stable")[-CONFIG["n_variable_genes"]:]]
dense = expr[variable].T.toarray()
dense -= dense.mean(axis=0)
dense /= np.maximum(dense.std(axis=0, ddof=1), np.finfo(float).eps)
dense = np.clip(dense, -10, 10)
_, vectors = eigh(dense.T @ dense, subset_by_index=[len(variable) - CONFIG["n_pcs"], len(variable) - 1])
pcs = dense @ vectors[:, ::-1]
tree = linkage(pcs, method="ward", optimal_ordering=False)
meta["cluster"] = fcluster(tree, t=CONFIG["n_clusters"], criterion="maxclust")
levels = np.sort(meta.cluster.unique())
membership = sparse.csc_matrix((np.ones(len(meta)), (np.arange(len(meta)), np.searchsorted(levels, meta.cluster))), shape=(len(meta), len(levels)))
cluster_sizes = np.asarray(membership.sum(axis=0)).ravel()
means = np.asarray((expr @ membership).todense()) / cluster_sizes
fractions = np.asarray(((counts > 0) @ membership).todense()) / cluster_sizes
rest_means = (np.asarray(expr.sum(axis=1)).ravel()[:, None] - means * cluster_sizes) / (len(meta) - cluster_sizes)
markers = pd.DataFrame({"gene_id": np.tile(features[0].to_numpy(), len(levels)), "gene": np.tile(gene_symbols, len(levels)), "cluster": np.repeat(levels, len(features)), "mean_log1p": means.T.ravel(), "pct": fractions.T.ravel(), "mean_rest": rest_means.T.ravel()})
markers["mean_difference"] = markers.mean_log1p - markers.mean_rest
markers["rank"] = markers.groupby("cluster").mean_difference.rank(ascending=False, method="min")
markers.loc[markers["rank"].le(25)].to_csv(out / "cluster_top_markers.csv", index=False, encoding="utf-8")
panel = ["Pdgfra", "Col1a1", "Dcn", "Lum", "Scx", "Tnmd", "Pecam1", "Cdh5", "Kdr", "Rgs5", "Acta2", "Myh11", "Ptprc", "Lyz2", "Csf1r", "Adgre1", "C1qa", "C1qb", "Apoe", "Gpnmb", "Spp1", "Cd68", "Mrc1", "Cd163", "Folr2", "Lyve1", "Ccr2", "Ly6c2", "S100a8", "S100a9", "Csf3r", "Retnlg", "Cd74", "H2-Ab1", "Flt3", "Clec9a", "Cd3d", "Cd3e", "Cd79a", "Ncr1", "Pax7", "Myod1", "Myog", "Des", "Sox10", "Plp1", "Il6ra", "Il6st"]
markers.loc[markers.gene.isin(panel)].to_csv(out / "cluster_marker_panel.csv", index=False, encoding="utf-8")
meta.to_csv(out / "cell_clusters.csv.gz", index=False, encoding="utf-8")
sparse.save_npz(out / "retained_counts.npz", counts)
features.to_csv(out / "features.tsv", sep="\t", header=False, index=False, encoding="utf-8")
pd.crosstab(meta.cluster, meta.sample_id).to_csv(out / "cluster_sample_counts.csv", encoding="utf-8")
print(pd.crosstab(meta.cluster, meta.sample_id).to_string(), flush=True)
print(markers.loc[markers["rank"].le(8)].groupby("cluster").gene.agg(",".join).to_string(), flush=True)
