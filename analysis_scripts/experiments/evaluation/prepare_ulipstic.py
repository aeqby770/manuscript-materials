import os
from pathlib import Path
from io import StringIO
CONFIG={"paths":{"root":os.environ.get("LLR_WORK_ROOT", str(Path(__file__).resolve().parents[2])),"input":"inputs/external_validation/GSE253000",
 "output":"results/evaluation/calendar_lineage_validation"},
 "threads":os.cpu_count(),"plot":{},
 "hash_fraction":0.8,"donor_threshold":0.008,"contact_threshold":0.0085,
 "author_qc":{"min_genes":200,"max_mt_fraction":0.15,"max_counts":{"SNH1":80000,"SNH3":75000,"SNH4":75000}},
 "scope":"All deposited LCMV libraries SNH1/SNH3/SNH4 retained; longitudinal primary cohort is all author-QC eligible unsorted WT mediastinal LN cells",
 "primary_unit":"Six individual mice, two at each observed time 36/50/96 hours; enriched sorts not pooled into time-dependent contact fractions",
 "fit_features":"RNA only. Known experimental FLAG donor identity is supplied equally as source annotation; biotin is absent from model features and metadata.",
 "truth":"Independent biotin contact label at the fixed author threshold, plus continuous hashtag-normalized biotin",
 "metrics":["Per-mouse recipient-state AUROC/AUPRC","sample-level contact direction and peak interval","absolute calendar error and observation resolution","all planned scope coverage and unavailable outputs"],
 "time_resolution":"Two-hour labeling pulse and 36/50/96h sampling intervals",
 "limitations":"Physical contact is not a specific ligand-receptor causal effect; no true genealogy labels",
 "gene_identity":"Duplicate symbols resolve only through one-to-one official Entrez mapping. Same-gene loci are summed; unresolved loci retain unique Ensembl IDs. Every RNA count is retained.",
 "metadata_conflict":"SNH1 M6_biotin is called WT in readme but shares M6 identity with a mutant unsorted sample; retain as ambiguous, outside the unsorted WT primary cohort",
 "selection":"Experimental scope, thresholds and metrics fixed before computing model validation scores"}
os.environ["OMP_NUM_THREADS"]=str(CONFIG["threads"])
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
import anndata as ad
import matplotlib.pyplot as plt
plt.rcParams.update(CONFIG["plot"])
root=Path(CONFIG["paths"]["root"])
dest=root/CONFIG["paths"]["input"]
out=root/CONFIG["paths"]["output"]
f1=pd.read_csv(dest/"raw/GSM8012223_SNH1_features.tsv.gz",sep="\t",header=None,names=["id","symbol","feature_type"])
b1=pd.read_csv(dest/"raw/GSM8012223_SNH1_barcodes.tsv.gz",sep="\t",header=None)[0]
x1=sio.mmread(dest/"raw/GSM8012223_SNH1_matrix.mtx.gz").tocsr().T.tocsr()
p1=pd.DataFrame(x1[:,f1.feature_type.ne("Gene Expression").to_numpy()].toarray(),
 columns=f1.loc[f1.feature_type.ne("Gene Expression").to_numpy(),"symbol"])
h1=p1.columns.difference(["biotin_interaction","donor"])
o1=pd.DataFrame({"library":"SNH1","barcode":b1,"sample_hash":p1[h1].idxmax(axis=1),
 "winning_hash":p1[h1].max(axis=1),"total_hash":p1[h1].sum(axis=1),
 "biotin_counts":p1.biotin_interaction,"flag_counts":p1.donor})
a1=ad.AnnData(X=x1[:,f1.feature_type.eq("Gene Expression").to_numpy()],obs=o1,
 var=f1.loc[f1.feature_type.eq("Gene Expression").to_numpy(),["id","symbol"]].set_index("id"))
a1.obs_names="SNH1:"+b1
f3=pd.read_csv(dest/"raw/GSM8012225_SNH3_features.tsv.gz",sep="\t",header=None,names=["id","symbol","feature_type"])
b3=pd.read_csv(dest/"raw/GSM8012225_SNH3_barcodes.tsv.gz",sep="\t",header=None)[0]
x3=sio.mmread(dest/"raw/GSM8012225_SNH3_matrix.mtx.gz").tocsr().T.tocsr()
p3=pd.DataFrame(x3[:,f3.feature_type.ne("Gene Expression").to_numpy()].toarray(),
 columns=f3.loc[f3.feature_type.ne("Gene Expression").to_numpy(),"symbol"])
h3=p3.columns.difference(["biotin_interaction","donor"])
o3=pd.DataFrame({"library":"SNH3","barcode":b3,"sample_hash":p3[h3].idxmax(axis=1),
 "winning_hash":p3[h3].max(axis=1),"total_hash":p3[h3].sum(axis=1),
 "biotin_counts":p3.biotin_interaction,"flag_counts":p3.donor})
a3=ad.AnnData(X=x3[:,f3.feature_type.eq("Gene Expression").to_numpy()],obs=o3,
 var=f3.loc[f3.feature_type.eq("Gene Expression").to_numpy(),["id","symbol"]].set_index("id"))
a3.obs_names="SNH3:"+b3
f4=pd.read_csv(dest/"raw/GSM8012227_SNH4_features.tsv.gz",sep="\t",header=None,names=["id","symbol","feature_type"])
b4=pd.read_csv(dest/"raw/GSM8012227_SNH4_barcodes.tsv.gz",sep="\t",header=None)[0]
x4=sio.mmread(dest/"raw/GSM8012227_SNH4_matrix.mtx.gz").tocsr().T.tocsr()
p4=pd.DataFrame(x4[:,f4.feature_type.ne("Gene Expression").to_numpy()].toarray(),
 columns=f4.loc[f4.feature_type.ne("Gene Expression").to_numpy(),"symbol"])
h4=p4.columns.difference(["biotin_interaction","donor"])
o4=pd.DataFrame({"library":"SNH4","barcode":b4,"sample_hash":p4[h4].idxmax(axis=1),
 "winning_hash":p4[h4].max(axis=1),"total_hash":p4[h4].sum(axis=1),
 "biotin_counts":p4.biotin_interaction,"flag_counts":p4.donor})
a4=ad.AnnData(X=x4[:,f4.feature_type.eq("Gene Expression").to_numpy()],obs=o4,
 var=f4.loc[f4.feature_type.eq("Gene Expression").to_numpy(),["id","symbol"]].set_index("id"))
a4.obs_names="SNH4:"+b4
assert a1.var_names.equals(a3.var_names)
assert a1.var_names.equals(a4.var_names)
assert a1.var_names.is_unique
all_rna=ad.concat([a1,a3,a4],axis=0,merge="same")
original_n_genes=all_rna.X.getnnz(axis=1)
original_library_size=np.asarray(all_rna.X.sum(axis=1)).ravel()
gene_map=pd.read_csv(dest/"rna_feature_identity_mapping.csv",dtype=str)
gene_map.index=pd.Index(gene_map["id"])
gene_map=gene_map.reindex(all_rna.var_names)
assert gene_map.feature_key.notna().all()
codes,keys=pd.factorize(gene_map.feature_key,sort=False)
collapse=sp.csr_matrix((np.ones(len(codes),dtype=np.int64),(np.arange(len(codes)),codes)),
 shape=(len(codes),len(keys)))
var=gene_map.groupby("feature_key",sort=False,observed=True).agg(
 original_ids=("id",";".join),provided_symbol=("symbol","first"),n_input_loci=("id","size")).reindex(keys)
all_rna=ad.AnnData(X=(all_rna.X@collapse).tocsr(),obs=all_rna.obs.copy(),var=var)
assert np.array_equal(original_library_size,np.asarray(all_rna.X.sum(axis=1)).ravel())
obs=all_rna.obs.copy()
obs["barcode_id"]=obs.index
obs["sample_hash"]=obs.sample_hash.str.replace("M13A","M13",regex=False)
obs["mouse"]=obs.sample_hash.str.replace("_biotin","",regex=False)
obs["enriched_sort"]=obs.sample_hash.str.endswith("_biotin")
design=pd.read_csv(StringIO("""mouse,time_hours,background,tissue
M1,96,WT,MedLN
M3,96,WT,MedLN
M4,96,mut,MedLN
M6,50,mut,MedLN
M8,50,WT,MedLN
M9,96,mut,MedLN
M11,36,WT,MedLN
M13,36,WT,MedLN
M14,36,mut,MedLN
M16,50,WT,MedLN
Li1,96,WT,Liver
Li2,96,WT,Liver
Li3,96,WT,Liver
Li4,96,mut,Liver
Li5,96,mut,Liver
Lu1,96,WT,Lung
Lu2,96,WT,Lung
Lu3,96,WT,Lung
Lu4,96,mut,Lung
Lu5,96,mut,Lung
S1,96,WT,Spleen
S2,96,WT,Spleen
S3,96,WT,Spleen
S4,96,mut,Spleen
S5,96,mut,Spleen
""")).set_index("mouse")
obs=obs.join(design,on="mouse")
assert obs.time_hours.notna().all()
obs.loc[obs.library.eq("SNH1")&obs.sample_hash.eq("M6_biotin"),"background"]="ambiguous_readme"
obs["hash_fraction"]=np.divide(obs.winning_hash,obs.total_hash,
 out=np.zeros(len(obs)),where=obs.total_hash.gt(0))
biotin=np.divide(obs.biotin_counts,obs.winning_hash,out=np.zeros(len(obs)),where=obs.winning_hash.gt(0))
flag=np.divide(obs.flag_counts,obs.winning_hash,out=np.zeros(len(obs)),where=obs.winning_hash.gt(0))
labels=pd.DataFrame({"cell_id":obs.index,"biotin_counts":obs.biotin_counts,
 "sample_hash_count":obs.winning_hash,"normalized_biotin":biotin,
 "contact_positive":biotin>=CONFIG["contact_threshold"],"normalized_flag":flag})
obs["cell_type"]=np.where(flag>=CONFIG["donor_threshold"],"donor_P14","recipient")
obs["total_counts"]=np.asarray(all_rna.X.sum(axis=1)).ravel()
obs["n_genes"]=original_n_genes
obs["mt_fraction"]=np.divide(np.asarray(all_rna.X[:,all_rna.var_names.str.startswith("mt-")].sum(axis=1)).ravel(),
 obs.total_counts,out=np.zeros(len(obs)),where=obs.total_counts.gt(0))
obs["qc_pass"]=obs.n_genes.ge(CONFIG["author_qc"]["min_genes"])&obs.mt_fraction.le(CONFIG["author_qc"]["max_mt_fraction"])&obs.total_counts.le(obs.library.map(CONFIG["author_qc"]["max_counts"]))&obs.hash_fraction.ge(CONFIG["hash_fraction"])
obs["sample_id"]=obs.tissue.astype(str)+":"+obs.background.astype(str)+":"+obs.mouse.astype(str)
obs["time_info"]=obs.time_hours
obs["primary_calendar_eligible"]=obs.qc_pass&obs.tissue.eq("MedLN")&obs.background.eq("WT")&~obs.enriched_sort
labels.to_csv(dest/"GSE253000_independent_contact_labels.csv",index=False,encoding="utf-8")
obs.to_csv(dest/"GSE253000_all_called_cell_metadata.csv",index=True,encoding="utf-8")
all_rna.obs=obs.drop(columns=["biotin_counts"])
all_rna.write_h5ad(dest/"GSE253000_all_lcmv_rna.h5ad",compression="gzip")
primary=all_rna[obs.primary_calendar_eligible].copy()
primary.obs=primary.obs.drop(columns=["flag_counts","winning_hash","total_hash","hash_fraction"])
primary.write_h5ad(dest/"GSE253000_ulipstic_mLN_unsorted.h5ad",compression="gzip")
summary=obs.groupby(["library","sample_hash","mouse","tissue","time_hours","background","enriched_sort"],
 observed=True,dropna=False).agg(n_called=("barcode_id","size"),n_qc=("qc_pass","sum"),n_primary=("primary_calendar_eligible","sum")).reset_index()
summary.to_csv(out/"ulipstic_complete_sample_inventory.csv",index=False,encoding="utf-8")
primary_counts=primary.obs.groupby(["sample_id","time_hours","cell_type"],observed=True).size().rename("n_cells").reset_index()
primary_counts.to_csv(out/"ulipstic_primary_sample_inventory.csv",index=False,encoding="utf-8")
print("All LCMV RNA",all_rna.shape,"primary",primary.shape,flush=True)
print(primary_counts.to_string(index=False),flush=True)
