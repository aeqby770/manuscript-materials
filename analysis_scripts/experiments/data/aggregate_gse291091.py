import os
from pathlib import Path
CONFIG = {"paths": {"root": os.environ.get("LLR_WORK_ROOT", str(Path(__file__).resolve().parents[2]))}, "threads": os.cpu_count(), "plot": {}}
os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
os.environ["OPENBLAS_NUM_THREADS"] = str(CONFIG["threads"])
import numpy as np
import pandas as pd
from scipy import sparse

root = Path(CONFIG["paths"]["root"])
out = root / "results/validation/gse291091"
meta = pd.read_csv(out / "cell_clusters.csv.gz")
features = pd.read_csv(out / "features.tsv", sep="\t", header=None)
counts = sparse.load_npz(out / "retained_counts.npz")
annotation = pd.DataFrame({
    "cluster": np.arange(1, 21),
    "cell_type": ["Endothelial", "Myogenic", "FAPs", "FAPs", "C1qa macrophages", "C1qa macrophages", "Gpnmb macrophages", "Gpnmb macrophages", "Gpnmb macrophages", "Neutrophils", "MNP unresolved", "APC", "MNP activated", "C1qa macrophages", "Gpnmb macrophages", "Inflammatory MNP", "MNP IFN state", "Inflammatory MNP", "Inflammatory MNP", "Inflammatory MNP"],
    "positive_evidence": ["Pecam1,Cdh5,Egfl7,Fabp4", "Des,Myod1,Pax7", "Col1a1,Col1a2,Col3a1,Dcn", "Col1a1,Col1a2,Pdgfra,Dcn", "C1qa,Csf1r,Adgre1,Lyz2", "C1qa,Csf1r,Adgre1,Mrc1,Folr2", "Gpnmb,C1qa,Adgre1,Csf1r", "Gpnmb,C1qa,Adgre1,Csf1r", "Gpnmb,Spp1,Lyz2,Adgre1", "S100a8,S100a9,Retnlg,Il1r2", "Lyz2,C1qa,Csf1r", "Cd74,H2-Ab1,H2-Aa", "Lyz2,Ly6c2,Cd74,Csf1r", "C1qa,C1qb,Csf1r,Adgre1", "Gpnmb,C1qa,Adgre1,Csf1r", "Ly6c2,Ccr2,Lyz2,Plac8", "Lyz2,C1qa,Gpnmb,Isg15,Irf7", "Ly6c2,Ccr2,Lyz2,Plac8", "Ly6c2,Ccr2,Lyz2,Plac8", "Ly6c2,Ccr2,Lyz2,Plac8"],
    "negative_evidence": ["Ptprc,Csf1r,Pdgfra", "Ptprc,Pecam1,Pdgfra", "Pecam1,Cdh5,Csf1r", "Pecam1,Cdh5,Csf1r", "Pdgfra,Cdh5,Retnlg", "Pdgfra,Cdh5,Retnlg", "Pdgfra,Cdh5,Retnlg", "Pdgfra,Cdh5,Retnlg", "Pdgfra,Cdh5,Retnlg", "Csf1r,Adgre1,C1qa", "Pdgfra,Cdh5,Retnlg", "Adgre1,Gpnmb,Mrc1", "Pdgfra,Cdh5,Retnlg", "Pdgfra,Cdh5,Retnlg", "Pdgfra,Cdh5,Retnlg", "Pdgfra,Cdh5,Retnlg", "Pdgfra,Cdh5,Retnlg", "Pdgfra,Cdh5,Retnlg", "Pdgfra,Cdh5,Retnlg", "Pdgfra,Cdh5,Retnlg"]
})
annotation["source"] = "10.1002/advs.202507444;10.1016/j.celrep.2020.02.067"
annotation["interpretation"] = "Marker-supported identity/state; not a one-to-one M1/M2 mapping"
annotation.to_csv(out / "cluster_annotation.csv", index=False, encoding="utf-8")
meta = meta.merge(annotation[["cluster", "cell_type"]], on="cluster", how="left", validate="many_to_one", sort=False)
assert len(meta) == counts.shape[1]
meta.to_csv(out / "cell_metadata.csv.gz", index=False, encoding="utf-8")
mnp = meta.cell_type.isin(["C1qa macrophages", "Gpnmb macrophages", "MNP unresolved", "APC", "MNP activated", "Inflammatory MNP", "MNP IFN state"])
members = pd.concat([
    meta.assign(cell_index=np.arange(len(meta)), context=meta.cell_type)[["cell_index", "sample_id", "condition", "context"]],
    meta.loc[mnp].assign(cell_index=np.flatnonzero(mnp), context="All MNP")[["cell_index", "sample_id", "condition", "context"]]
], ignore_index=True)
members["context_sample"] = members.context + ":" + members.sample_id
units = members[["context_sample", "context", "sample_id", "condition"]].drop_duplicates().sort_values(["context", "sample_id"]).reset_index(drop=True)
unit_index = pd.Index(units.context_sample).get_indexer(members.context_sample)
membership = sparse.csc_matrix((np.ones(len(members)), (members.cell_index, unit_index)), shape=(len(meta), len(units)))
pb = (counts @ membership).toarray()
symbols = pd.Index(features[1].astype(str)).unique()
gene_map = sparse.csr_matrix((np.ones(len(features)), (symbols.get_indexer(features[1].astype(str)), np.arange(len(features)))), shape=(len(symbols), len(features)))
pb = gene_map @ pb
pd.DataFrame(pb, index=symbols, columns=units.context_sample).rename_axis("gene").to_csv(out / "pseudobulk_counts.csv.gz", encoding="utf-8")
units["n_cells"] = np.asarray(membership.sum(axis=0)).ravel().astype(int)
units.to_csv(out / "pseudobulk_samples.csv", index=False, encoding="utf-8")
print(units.to_string(index=False))
