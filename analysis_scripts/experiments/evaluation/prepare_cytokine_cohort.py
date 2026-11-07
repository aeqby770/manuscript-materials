import os
from pathlib import Path
CONFIG={"paths":{"root":os.environ.get("LLR_WORK_ROOT", str(Path(__file__).resolve().parents[2])),"input":"inputs/external_validation/LARRY",
 "output":"results/evaluation/calendar_lineage_validation"},
 "threads":os.cpu_count(),"plot":{},"dataset":"larry_cytokine_perturbation",
 "scope":"All 65075 published cells and all cytokine conditions; no subsampling or seed expansion",
 "primary_object":"Day4 cells: author-assigned clone descendants at day6 within the same cytokine condition",
 "secondary_object":"Day2 shared starting cells: condition-specific clone descendant distributions at day6; retrospective fate inference with later RNA observed, not prospective prediction of an unseen intervention",
 "fit_information":"Published RNA matrix and experimental time/condition only; no clone IDs or clone embedding. Terminal RNA annotations are provided equally to all evaluators as fate label mapping.",
 "independent_truth":"Author LARRY clone matrix; mixed fates remain distributions",
 "unit":"Clone within condition. Sequencing libraries are technical partitions and not independent biological replicates.",
 "missing":"SCF has day4 RNA but no day6 sample; preserve as unavailable future truth rather than infer a fate label",
 "primary_metrics":["clone-balanced multiclass Brier score","cross entropy","coverage by time and condition","per-condition and complete pooled reports"],
 "controls":["current geometric branch weights","current local RNA priming weights","hard assignment","same-graph fate probabilities","condition-specific terminal RNA prior"],
 "selection":"Scope and metrics fixed from the actual design before computing any cytokine-cohort validation scores",
 "input_scale":"Author-provided L1-normalized transcript counts, not raw UMIs",
 "tracing_features":"No barcode identity features in the gene matrix; Gfpt1/Gfpt2 are endogenous genes and are retained"}
os.environ["OMP_NUM_THREADS"]=str(CONFIG["threads"])
import numpy as np
import pandas as pd
import scipy.io as sio
import anndata as ad
import matplotlib.pyplot as plt
plt.rcParams.update(CONFIG["plot"])
root=Path(CONFIG["paths"]["root"])
dest=root/CONFIG["paths"]["input"]
out=root/CONFIG["paths"]["output"]
meta=pd.read_csv(dest/"stateFate_cytokinePerturbation_metadata.txt.gz",sep="\t")
genes=pd.read_csv(dest/"stateFate_cytokinePerturbation_gene_names.txt.gz",header=None)[0]
matrix=sio.mmread(dest/"stateFate_cytokinePerturbation_normed_counts.mtx.gz").tocsr()
clone=sio.mmread(dest/"stateFate_cytokinePerturbation_clone_matrix.mtx.gz").tocsr().T.tocsr()
assert matrix.shape==(len(meta),len(genes))
assert clone.shape[0]==len(meta)
cells=meta["Library"]+":"+meta["Cell barcode"].str.replace("-","",regex=False)
assert cells.is_unique and genes.is_unique
obs=pd.DataFrame({"barcode_id":cells,"time_info":meta["Time point"],
 "condition":meta["Cytokine condition"],"sequencing_library":meta["Library"]})
obs.index=pd.Index(cells)
annotations=pd.DataFrame({"cell_id":cells,"terminal_rna_state":meta["Cell type annotation"],
 "time_days":meta["Time point"],"condition":meta["Cytokine condition"],"sequencing_library":meta["Library"]})
annotations.to_csv(dest/"larry_cytokine_independent_annotations.csv",index=False,encoding="utf-8")
rna=ad.AnnData(X=matrix,obs=obs,var=pd.DataFrame(index=pd.Index(genes,name="gene")))
rna.uns["input_scale"]="published L1-normalized transcript counts"
rna.write_h5ad(dest/(CONFIG["dataset"]+".h5ad"),compression="gzip")
coverage=annotations.groupby(["time_days","condition"],dropna=False).size().rename("n_cells").reset_index()
coverage.to_csv(out/"larry_cytokine_complete_scope.csv",index=False,encoding="utf-8")
print("RNA",rna.shape,"clones",clone.shape,"memberships",clone.nnz,flush=True)
print(coverage.to_string(index=False),flush=True)
