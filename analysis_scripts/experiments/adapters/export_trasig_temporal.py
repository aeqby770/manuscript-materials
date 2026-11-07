import os
from pathlib import Path

CONFIG = {
    "paths": {"revision": os.environ.get("LLR_WORK_ROOT", str(Path(__file__).resolve().parents[2]))},
    "threads": os.cpu_count(),
    "seed": 12345,
    "plot": {},
    "window": 20,
    "export_modes": {"GSE274389_signal_recorder": "wnt",
                     "GSE143437_muscle_regeneration": "fap"},
    "output_files": {"wnt": "trasig_wnt_temporal.csv", "fap": "trasig_fap_temporal.csv",
                     "all": "trasig_temporal.csv"},
    "fap_targets": ["Monocytes/Macrophages/Platelets", "Pro-inflammatory macrophages",
                    "Anti-inflammatory macrophages", "Resident Macrophages/APCs",
                    "MuSCs and progenitors"],
}

from pathlib import Path
import pickle
import sys

import bottleneck as bn
import numpy as np
import pandas as pd

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])

dataset = sys.argv[1]
algorithm_directory = (sys.argv[2:] + ["TraSig"])[0]
export_stage = (sys.argv[3:] + ["complete"])[0]
directory = (Path(sys.argv[4]) if len(sys.argv)>4 else Path(CONFIG["paths"]["revision"]) / "results/by_algorithm" / algorithm_directory / dataset)
expression = pd.read_csv(directory / f"{dataset}_log1p_lr.txt", index_col=0)
inputs = pd.read_csv(directory / "trasig_cell_inputs.csv").set_index("cell_id")
inputs = inputs.loc[expression.index]
labels = pd.read_csv(directory / "trasig_path_labels.csv", dtype={"cell_type": str}).set_index("path")
match export_stage:
    case "complete":
        native = pd.read_csv(directory / "trasig_scores.csv")
    case "observed":
        observed = pickle.loads((directory / "results" /
            f"{dataset}_lr_log1p_slingshot_smallerWindow_nlap_20_metrics_0.pickle").read_bytes())
        lr_pairs = np.asarray(pickle.loads((directory / f"lr_{dataset}_log1p.pickle").read_bytes()))
        native = pd.DataFrame({
            "path_pair": np.repeat(list(observed), len(lr_pairs)),
            "ligand": np.tile(lr_pairs[:, 0], len(observed)),
            "receptor": np.tile(lr_pairs[:, 1], len(observed)),
            "observed_score": np.concatenate(pd.DataFrame.from_dict(observed, orient="index")["dot"].to_numpy()),
            "permutation_p_value": np.nan,
        })
native["inference_status"] = export_stage
assert expression.index.equals(inputs.index)
assert not native.duplicated(["path_pair", "ligand", "receptor"]).any()

# Exact native smallerWindow/nLap=20 path; absent time bins are NaN before
# rolling sums, and absent windows become zero after the native nan2zero step.
paths = np.sort(inputs["cell_path"].unique())
times = np.round(np.arange(0, 1.01, 0.01), 2)
bin_index = pd.MultiIndex.from_product([times, paths], names=["time", "path"])
grouped = expression.assign(time=np.round(inputs["cell_time"], 2), path=inputs["cell_path"]).groupby(["time", "path"])
summed = grouped.sum().reindex(bin_index).to_numpy().reshape(len(times), len(paths), expression.shape[1])
counts = grouped.size().reindex(bin_index).to_numpy().reshape(len(times), len(paths))
half = CONFIG["window"] // 2
small_count = bn.move_sum(counts, window=half, min_count=1, axis=0)
full_count = bn.move_sum(counts, window=CONFIG["window"], min_count=1, axis=0)
small = np.divide(
    bn.move_sum(summed, window=half, min_count=1, axis=0),
    small_count[:, :, None],
    out=np.zeros_like(summed), where=np.isfinite(small_count[:, :, None]) & (small_count[:, :, None] > 0),
)
full = np.divide(
    bn.move_sum(summed, window=CONFIG["window"], min_count=1, axis=0),
    full_count[:, :, None],
    out=np.zeros_like(summed), where=np.isfinite(full_count[:, :, None]) & (full_count[:, :, None] > 0),
)
smooth = np.concatenate([small[half:CONFIG["window"]], full[CONFIG["window"]:]], axis=0)
support = np.nan_to_num(np.concatenate([small_count[half:CONFIG["window"]], full_count[CONFIG["window"]:]], axis=0))
curve_times = times[half:]
path_pairs = native["path_pair"].str.split("_", n=1, expand=True, regex=False).astype(int)
source_index = pd.Index(paths).get_indexer(path_pairs[0])
target_index = pd.Index(paths).get_indexer(path_pairs[1])
ligand_index = expression.columns.get_indexer(native["ligand"])
receptor_index = expression.columns.get_indexer(native["receptor"])
assert min(source_index.min(), target_index.min(), ligand_index.min(), receptor_index.min()) >= 0
curves = smooth[:, source_index, ligand_index] * smooth[:, target_index, receptor_index]
recovered = curves.mean(axis=0)
error = np.abs(recovered - native["observed_score"].to_numpy())
# The author implementation writes scores rounded to ten decimal places.
assert error.max() <= 5.1e-11, error.max()
agreement = native[["path_pair", "ligand", "receptor", "observed_score"]].copy()
agreement["recovered_mean"] = recovered
agreement["absolute_error"] = error
agreement.to_csv(directory / "trasig_temporal_agreement.csv", index=False, encoding="utf-8", lineterminator="\n")

# Export the prespecified biological task, retaining every LR/source/target.
native["source_type"] = labels["cell_type"].reindex(path_pairs[0]).to_numpy()
native["target_type"] = labels["cell_type"].reindex(path_pairs[1]).to_numpy()
assert not native[["source_type", "target_type"]].isna().any().any()
selection = {
    "wnt": native["ligand"].str.upper().str.startswith("WNT"),
    "fap": native["source_type"].eq("FAPs") & native["target_type"].isin(CONFIG["fap_targets"]),
    "all": np.ones(len(native), dtype=bool),
}
export_mode = CONFIG["export_modes"].get(dataset, "all")
row_index = np.flatnonzero(selection[export_mode])
selected = native.iloc[row_index].copy()
selected["rank_native_score"] = selected.groupby(
    ["source_type", "target_type"], observed=True
)["observed_score"].rank(method="min", ascending=False)
selected.to_csv(directory / "trasig_selected_task_pairs.csv", index=False,
                encoding="utf-8", lineterminator="\n")
rows = selected.iloc[np.tile(np.arange(len(selected)), len(curve_times))].reset_index(drop=True)
rows["time"] = np.repeat(curve_times, len(selected))
rows["score"] = curves[:, row_index].ravel()
rows["source_cells_in_window"] = support[:, source_index[row_index]].ravel()
rows["target_cells_in_window"] = support[:, target_index[row_index]].ravel()
rows.to_csv(directory / CONFIG["output_files"][export_mode], index=False, encoding="utf-8", lineterminator="\n")
print(f"native_rows={len(native)} paths={len(paths)} exported_pairs={selected[['ligand', 'receptor']].drop_duplicates().shape[0]} max_error={error.max():.3g}")
