#!/usr/bin/env python
"""NATMI ExtractEdges.py adapter for canonical h5ad input.

Usage: python run_natmi.py <canonical.h5ad> <output_edges.csv> <species>
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import _mem_guard
_mem_guard.cap_self(55)

CONFIG = {
    "paths": {
        "vendor": Path(os.environ["LLR_VENDOR_ROOT"]).resolve() / "NATMI",
    },
    "threads": os.cpu_count(),
    "seed": 12345,
    "plot": {
        "dpi": 300,
        "font": "Arial",
        "fontsize": 12,
        "palette": "viridis",
    },
    "species_map": {
        "human": "human",
        "mouse": "mouse",
        "rat": "rat",
        "zebrafish": "zebrafish",
        "caenorhabditis": "nematode",
        "cao(c.elegans)": "nematode",
    },
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])

if len(sys.argv) != 4:
    print("Usage: run_natmi.py <canonical.h5ad> <output_edges.csv> <species>", file=sys.stderr)
    sys.exit(2)

inpath, outpath, species_raw = sys.argv[1], sys.argv[2], sys.argv[3].lower()
species = CONFIG["species_map"].get(species_raw)
if species is None:
    print(f"NATMI requires a supported species for connectomeDB2020; got: {species_raw}", file=sys.stderr)
    sys.exit(1)

adata = sc.read_h5ad(inpath)
group_col = "cell_type" if "cell_type" in adata.obs.columns else "clusters"
if group_col not in adata.obs.columns and "cluster" in adata.obs.columns:
    group_col = "cluster"
if group_col not in adata.obs.columns:
    group_col = adata.obs.columns[0]

if "lognorm" in adata.layers:
    x = adata.layers["lognorm"]
elif "counts" in adata.layers:
    x = adata.layers["counts"].copy()
    adata.X = x.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    x = adata.X
else:
    x = adata.X.copy()

x = x.toarray() if sparse.issparse(x) else np.asarray(x)
x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)

with tempfile.TemporaryDirectory() as tmp:
    em_file = os.path.join(tmp, "natmi_expression.csv")
    ann_file = os.path.join(tmp, "natmi_annotation.csv")
    natmi_out = os.path.join(tmp, "natmi_out")

    pd.DataFrame(x.T, index=adata.var_names.astype(str), columns=adata.obs_names.astype(str)).to_csv(em_file)
    labels = adata.obs[group_col].astype(str).values
    labels = np.where(pd.Series(labels).str.fullmatch(r"\d+").fillna(False), "C" + labels, labels)
    pd.DataFrame({"cluster": labels}, index=adata.obs_names.astype(str)).to_csv(ann_file)

    cmd = [
        sys.executable,
        "-c",
        "import pandas as pd, runpy; pd.ExcelWriter.save = pd.ExcelWriter.close; runpy.run_path('ExtractEdges.py', run_name='__main__')",
        "--emFile",
        em_file,
        "--annFile",
        ann_file,
        "--species",
        species,
        "--idType",
        "symbol",
        "--interDB",
        "lrc2p",
        "--interSpecies",
        species,
        "--coreNum",
        str(CONFIG["threads"]),
        "--out",
        natmi_out,
    ]
    proc = subprocess.run(cmd, cwd=CONFIG["paths"]["vendor"], capture_output=True, text=True)
    if proc.returncode != 0:
        print(proc.stdout)
        print(proc.stderr, file=sys.stderr)
        sys.exit(proc.returncode)

    edge_path = os.path.join(natmi_out, "Edges_lrc2p.csv")
    if not os.path.exists(edge_path):
        print("NATMI produced no Edges_lrc2p.csv; no signaling interactions found", file=sys.stderr)
        sys.exit(1)

    edges = pd.read_csv(edge_path)
    if edges.empty:
        print("NATMI edge table is empty; no signaling interactions found", file=sys.stderr)
        sys.exit(1)
    edges.to_csv(outpath, index=False)

print(f"NATMI OK: {len(edges)} edges -> {outpath}")
