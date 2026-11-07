import os
import argparse
from pathlib import Path
CONFIG = {
    "paths": {
        "project": os.environ.get("LLR_WORK_ROOT", os.getcwd()),
        "predictions": "inputs/comparison/completed_real_ordering_predictions.csv.gz",
        "reference": "inputs/comparison/results/paper_tables/reference",
        "submission": "inputs/comparison/results/paper_tables/table1_main/table1_pseudotime_linealinkr064_submission.csv",
        "output": "inputs/comparison/results/paper_tables/completed_real_ordering.csv",
        "full_output": "results/evaluation/final065/saved_native/classic_full_native_ordering_metrics.csv",
        "native": "results/by_algorithm"
    },
    "threads": int(os.environ.get("LLR_NATIVE_THREADS", os.cpu_count())),
    "plot": {}
}
os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
os.environ["OPENBLAS_NUM_THREADS"] = str(CONFIG["threads"])
import numpy as np
import pandas as pd

project = Path(CONFIG["paths"]["project"])
reference = project / CONFIG["paths"]["reference"]
parser=argparse.ArgumentParser()
parser.add_argument("--native-full-results",action="store_true")
args=parser.parse_args()
match args.native_full_results:
    case False:
        predictions=pd.read_csv(project/CONFIG["paths"]["predictions"],encoding="utf-8",float_precision="round_trip")
    case True:
        native=project/CONFIG["paths"]["native"]
        files=list((native/"PHLOWER").glob("*/phlower_pseudotime.csv"))+list((native/"SLICER").glob("*/slicer_pseudotime.csv"))
        file_index=pd.DataFrame({"source_id":np.arange(len(files)),"source_file":pd.Series(files,dtype="string")})
        file_index["dataset"]=pd.Series(files).map(lambda p:p.parent.name)
        file_index["algorithm"]=pd.Series(files).map(lambda p:p.parent.parent.name)
        predictions=pd.concat(map(pd.read_csv,files),keys=file_index.source_id,names=["source_id","native_row"]).reset_index().merge(file_index,on="source_id",validate="many_to_one")
        coverage_files=pd.Series(files).map(lambda p:p.parent/"coverage.csv")
        coverage=pd.concat(map(pd.read_csv,coverage_files),ignore_index=True).rename(columns={"method":"algorithm"})
        predictions=predictions.merge(coverage[["dataset","algorithm","root_input_mode"]],on=["dataset","algorithm"],validate="many_to_one")
        predictions["host"]="server"
assert not predictions.duplicated(["dataset","algorithm","cell_id"]).any()
references = {
    "dentategyrus": pd.read_csv(reference / "ref_dentategyrus.csv", encoding="utf-8", float_precision="round_trip").rename(columns={"age(days)": "truth"})[["cell_id", "truth"]],
    "gastrulation_erythroid": pd.read_csv(reference / "ref_gastrulation_erythroid.csv", encoding="utf-8", float_precision="round_trip").rename(columns={"stage": "truth"})[["cell_id", "truth"]],
    "zebrafish_neural_dev": pd.read_csv(reference / "ref_zebrafish_neural_dev.csv", encoding="utf-8", float_precision="round_trip").rename(columns={"timepoint": "truth"})[["cell_id", "truth"]],
    "lung": pd.read_csv(reference / "ref_lung.csv", encoding="utf-8", float_precision="round_trip").rename(columns={"time_point": "truth"})[["cell_id", "truth"]],
    "larry_invitro": pd.read_csv(reference / "ref_larry_invitro.csv", encoding="utf-8", float_precision="round_trip").rename(columns={"time_info": "truth"})[["cell_id", "truth"]],
    "human_fetal_liver_organoid": pd.read_csv(reference / "ref_human_fetal_liver_organoid.csv", encoding="utf-8", float_precision="round_trip").rename(columns={"time_days": "truth"})[["cell_id", "truth"]],
    "human_fetal_pancreas": pd.read_csv(reference / "ref_human_fetal_pancreas.csv", encoding="utf-8", float_precision="round_trip").rename(columns={"pcw": "truth"})[["cell_id", "truth"]],
    "zebrafish_tracerseq": pd.read_csv(reference / "ref_zebrafish_tracerseq.csv", encoding="utf-8", float_precision="round_trip").rename(columns={"time_hours": "truth"})[["cell_id", "truth"]]
}
truth = pd.concat(references, names=["dataset", "row"]).reset_index().drop(columns="row")
truth["truth"] = truth.truth.astype("string").str.extract(r"([0-9]+(?:\.[0-9]+)?)", expand=False).astype(float)
keys = ["dataset", "algorithm"]
joined = predictions.merge(truth, on=["dataset", "cell_id"], how="left", validate="many_to_one", indicator=True)
joined["matched"] = joined["_merge"].eq("both")
joined["finite_prediction"] = np.isfinite(joined.pseudotime)
aligned = joined.groupby(keys).agg(
    n_prediction_cells=("cell_id", "size"), n_aligned=("matched", "sum"),
    n_finite_prediction=("finite_prediction", "sum"),
    source_host=("host", "first"), source_file=("source_file", "first")
)
aligned = aligned.join(truth.groupby("dataset").size().rename("n_reference"), on="dataset")
aligned=aligned.join(truth.assign(known_time=np.isfinite(truth.truth)).groupby("dataset").known_time.sum().rename("n_reference_known_time"),on="dataset")
valid = joined.loc[joined.matched & np.isfinite(joined.pseudotime) & np.isfinite(joined.truth)].copy()
grouped = valid.groupby(keys)
valid["p01"] = (valid.pseudotime - grouped.pseudotime.transform("min")) / (grouped.pseudotime.transform("max") - grouped.pseudotime.transform("min"))
valid["t01"] = (valid.truth - grouped.truth.transform("min")) / (grouped.truth.transform("max") - grouped.truth.transform("min"))
valid["squared_error"] = (valid.p01 - valid.t01) ** 2
correlation_index=pd.MultiIndex.from_frame(valid[keys].drop_duplicates().assign(metric_variable="p01"))
spearman=grouped[["p01","t01"]].corr(method="spearman").reindex(correlation_index)["t01"].droplevel(-1).rename("Spearman_signed")
pearson=grouped[["p01","t01"]].corr(method="pearson").reindex(correlation_index)["t01"].droplevel(-1).rename("Pearson")
metrics = aligned.join(grouped.size().rename("n_used")).join(spearman).join(pearson).join(valid.groupby(keys).squared_error.mean().pow(0.5).rename("RMSE")).reset_index()
metrics["n_used"]=metrics.n_used.fillna(0).astype(int)
metrics["cell_coverage"] = metrics.n_aligned / metrics.n_reference
metrics["Spearman"] = metrics.Spearman_signed
metrics["Spearman_absolute"] = metrics.Spearman_signed.abs()
metrics["direction_policy"] = "native_signed_no_truth_flip"
metrics["metric_source"]=np.where(args.native_full_results,"current_full_native_saved_results","reused_native_predictions_local_server_archive")
metrics["row_present"] = True
metrics["output_status"] = np.select(
    [metrics.n_reference.isna(),
     metrics.n_aligned.ne(metrics.n_reference) | metrics.n_prediction_cells.ne(metrics.n_reference),
     metrics.n_finite_prediction.ne(metrics.n_prediction_cells),
     np.isfinite(metrics.Spearman_signed) & np.isfinite(metrics.RMSE)],
    ["no_independent_observed_time_reference", "partial_cell_coverage", "nonfinite_prediction", "scored"],
    default="nonfinite_metric"
)
submission = pd.read_csv(project / CONFIG["paths"]["submission"], encoding="utf-8", float_precision="round_trip")
root_modes=pd.concat([submission[keys+["root_input_mode"]],
 predictions.reindex(columns=keys+["root_input_mode"]).dropna(subset=["root_input_mode"])],ignore_index=True).drop_duplicates(keys,keep="last")
metrics = metrics.merge(root_modes, on=keys, how="left", validate="one_to_one")
match args.native_full_results:
    case True:
        repair_meta=coverage.reindex(columns=keys+["implementation_repair","n_empty_suffix_sentinels","runtime_scope"])
        metrics=metrics.merge(repair_meta,on=keys,how="left",validate="one_to_one")
output_key={False:"output",True:"full_output"}[args.native_full_results]
metrics.to_csv(project/CONFIG["paths"][output_key],index=False,encoding="utf-8")
print(metrics.groupby("algorithm").agg(n_datasets=("dataset", "size"), mean_signed=("Spearman_signed", "mean"), mean_absolute=("Spearman_absolute", "mean"), mean_rmse=("RMSE", "mean"), min_cell_coverage=("cell_coverage", "min")).sort_values("mean_signed", ascending=False).to_string())
print(metrics.groupby("output_status").size().to_string())
print("aligned_cells", int(metrics.n_aligned.sum()), "used_cells", int(metrics.n_used.sum()))

