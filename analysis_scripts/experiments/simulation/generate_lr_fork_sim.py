# -*- coding: utf-8 -*-
# Same-source counts -> LR -> fork simulation generator (revision 0.6.5).
# One dataset per invocation: --scenario=<name>. See SCENARIOS for truth design.
import os
import sys
import json
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from pathlib import Path

CONFIG = {
    "paths": {
        "out_dir": str(Path(__file__).resolve().parents[2] / "inputs" / "revision" / "sim"),
    },
    "threads": os.cpu_count(),
    "seed": 12345,
    "plot": {},
    "n_cells": 1200,
    "n_donors": 4,
    "n_traj_genes": 200,
    "n_noise_genes": 976,
    "nb_theta": 8.0,
    "effect_amp": 1.2,
}
os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])

LR_PAIRS = [
    ("EGF", "EGFR"), ("FGF1", "FGFR1"), ("FGF2", "FGFR2"), ("PDGFA", "PDGFRA"),
    ("PDGFB", "PDGFRB"), ("IGF1", "IGF1R"), ("CCL5", "CCR1"), ("CXCL12", "CXCR4"),
    ("IL4", "IL4R"), ("EPO", "EPOR"), ("KITL", "KIT"), ("GDF15", "TGFBR2"),
]
N_PAIRS = len(LR_PAIRS)

# Segment layout per topology: (label, pt_lo, pt_hi, cell fraction, beta_a, beta_b)
# Original (hard) benchmark layout: short stem, strongly skewed branch-1
# density, moderate branch separation. Restored after round-1 debug tuning.
SEGMENTS_2BR = [
    ("stem",    0.03, 0.50, 0.28, 1.0, 1.0),
    ("branch0", 0.50, 1.00, 0.38, 1.0, 1.0),
    ("branch1", 0.50, 1.00, 0.34, 0.55, 1.25),  # density skew toward fork
]
SEGMENTS_NESTED = [
    ("stem",          0.03, 0.50, 0.24, 1.0, 1.0),
    ("branch0",       0.50, 1.00, 0.30, 1.0, 1.0),
    ("branch1_shared",0.50, 0.75, 0.16, 1.0, 1.0),
    ("branch1a",      0.75, 1.00, 0.15, 1.0, 1.0),
    ("branch1b",      0.75, 1.00, 0.15, 1.0, 1.0),
]
# LR pair roles: 'focal' = rise post-fork on focal branch only;
# 'shared' = identical rise on both branches (null for D(t));
# 'flat' = no post-fork change (null); 'nested' = rise only on branch1a post fork2.
ROLES = {
    "shared_trend":     ["shared"] * 6 + ["flat"] * 6,
    "focal_branch":     ["focal"] * 6 + ["shared"] * 3 + ["flat"] * 3,
    "donor_confounded": ["focal"] * 2 + ["shared"] * 6 + ["flat"] * 4,
    "dropout":          ["focal"] * 6 + ["shared"] * 3 + ["flat"] * 3,
    "sampling_thinning":["focal"] * 6 + ["shared"] * 3 + ["flat"] * 3,
    "nested_fork":      ["nested"] * 3 + ["shared"] * 4 + ["flat"] * 5,
    "transient_focal":  ["transient"] * 6 + ["shared"] * 3 + ["flat"] * 3,
    "donor_ligand_only":["focal"] * 2 + ["shared"] * 6 + ["flat"] * 4,
    "donor_receptor_only":["focal"] * 2 + ["shared"] * 6 + ["flat"] * 4,
    "donor_nested_fork": ["nested"] * 3 + ["shared"] * 4 + ["flat"] * 5,
    "donor_transient_focal": ["transient"] * 6 + ["shared"] * 3 + ["flat"] * 3,
}
SCENARIOS = {
    "shared_trend":     {"segments": SEGMENTS_2BR,    "nested": False, "donor_confound": False, "dropout": 0.0, "thin_q": None},
    "focal_branch":     {"segments": SEGMENTS_2BR,    "nested": False, "donor_confound": False, "dropout": 0.0, "thin_q": None},
    "donor_confounded": {"segments": SEGMENTS_2BR,    "nested": False, "donor_confound": True,  "dropout": 0.0, "thin_q": None},
    "dropout":          {"segments": SEGMENTS_2BR,    "nested": False, "donor_confound": False, "dropout": 0.25, "thin_q": None},
    "sampling_thinning":{"segments": SEGMENTS_2BR,    "nested": False, "donor_confound": False, "dropout": 0.0, "thin_q": 0.25},
    "nested_fork":      {"segments": SEGMENTS_NESTED, "nested": True,  "donor_confound": False, "dropout": 0.0, "thin_q": None},
    "transient_focal":  {"segments": SEGMENTS_2BR,    "nested": False, "donor_confound": False, "dropout": 0.0, "thin_q": None},
    "donor_ligand_only":{"segments": SEGMENTS_2BR,    "nested": False, "donor_confound": True,  "dropout": 0.0, "thin_q": None},
    "donor_receptor_only":{"segments": SEGMENTS_2BR,   "nested": False, "donor_confound": True,  "dropout": 0.0, "thin_q": None},
    "donor_nested_fork": {"segments": SEGMENTS_NESTED, "nested": True, "donor_confound": True, "dropout": 0.0, "thin_q": None},
    "donor_transient_focal": {"segments": SEGMENTS_2BR, "nested": False, "donor_confound": True, "dropout": 0.0, "thin_q": None},
}
FORK1_PT = 0.50
FORK2_PT = 0.75
FOCAL_BRANCH_2BR = "branch0"
FOCAL_BRANCH_NESTED = "branch1a"
PULSE_PLATEAU_START, PULSE_PLATEAU_END, PULSE_END = 0.62, 0.68, 0.85

if sys.argv[1:] == ["--describe"]:
    print(json.dumps({"fork_times": [FORK1_PT, FORK2_PT], "two_branch_segments": SEGMENTS_2BR,
        "nested_segments": SEGMENTS_NESTED, "roles": ROLES,
        "pulse": {"start": FORK1_PT, "plateau_start": PULSE_PLATEAU_START,
                  "plateau_end": PULSE_PLATEAU_END, "end": PULSE_END},
        "pairs": LR_PAIRS, "seed": CONFIG["seed"]}))
    raise SystemExit(0)


np.random.seed(CONFIG["seed"])
scenario = sys.argv[1].split("=", 1)[1]
seed_override = [a.split("=", 1)[1] for a in sys.argv[2:]
                 if a.startswith("--seed=")]
if seed_override:
    CONFIG["seed"] = int(seed_override[0])
    np.random.seed(CONFIG["seed"])
spec = SCENARIOS[scenario]
rng = np.random.default_rng(CONFIG["seed"] + list(SCENARIOS).index(scenario) * 7919)

# ---------------------------------------------------------------------------
# Cells: segment membership, pseudotime, donors, cell types.
# ---------------------------------------------------------------------------
segments = spec["segments"]
seg_labels = np.array([s[0] for s in segments])
seg_frac = np.array([s[3] for s in segments])
seg_idx = rng.choice(len(segments), size=CONFIG["n_cells"], p=seg_frac / seg_frac.sum())
pt_lo = np.array([s[1] for s in segments])[seg_idx]
pt_hi = np.array([s[2] for s in segments])[seg_idx]
beta_a = np.array([s[4] for s in segments])[seg_idx]
beta_b = np.array([s[5] for s in segments])[seg_idx]
pt = pt_lo + (pt_hi - pt_lo) * rng.beta(beta_a, beta_b)
branch = seg_labels[seg_idx]

donor = rng.integers(0, CONFIG["n_donors"], size=CONFIG["n_cells"])
br0_post = (branch == FOCAL_BRANCH_NESTED if spec["nested"] else branch == FOCAL_BRANCH_2BR)
override = spec["donor_confound"] & br0_post & (rng.random(CONFIG["n_cells"]) < 0.65)
donor = np.where(override, CONFIG["n_donors"] - 1, donor)
is_confounded_donor = donor == (CONFIG["n_donors"] - 1)
donor_label = np.array([f"donor_{chr(65 + d)}" for d in donor])

# Original two-state design: SRC expresses ligands, TGT receptors. The runner
# now accepts two-state metadata, so no third category is invented.
ct_code = (rng.random(CONFIG["n_cells"]) >= 0.5).astype(int)
cell_type = np.array(["SRC", "TGT"])[ct_code]
is_src = ct_code == 0
is_tgt = ct_code == 1
cell_id = np.array([f"{scenario}_cell_{i}" for i in range(CONFIG["n_cells"])])

# ---------------------------------------------------------------------------
# Gene programs on log-mean scale.
# ---------------------------------------------------------------------------
n_traj = CONFIG["n_traj_genes"]
n_noise = CONFIG["n_noise_genes"]
n_genes = n_traj + 2 * N_PAIRS + n_noise

ramp = lambda x, lo, hi: np.clip((x - lo) / (hi - lo), 0.0, 1.0)
on = lambda *labels: np.isin(branch, list(labels))

traj_specs = [("early", 60, -pt / 0.12),
              ("shared", 40, 1.2 * ramp(pt, FORK1_PT, 0.85) * on("branch0", "branch1",
                                                                  "branch1_shared", "branch1a", "branch1b"))]
traj_specs += [("branch0", 35, 1.4 * ramp(pt, FORK1_PT, 0.80) * (branch == "branch0"))]
traj_specs += [("branch1", 35, 1.4 * ramp(pt, FORK1_PT, 0.80) * on("branch1", "branch1_shared",
                                                                 "branch1a", "branch1b"))] if not spec["nested"] else \
              [("branch1_shared", 20, 1.3 * ramp(pt, FORK1_PT, 0.72) * on("branch1_shared", "branch1a", "branch1b")),
               ("branch1a", 20, 1.4 * ramp(pt, FORK2_PT, 0.95) * (branch == "branch1a")),
               ("branch1b", 20, 1.4 * ramp(pt, FORK2_PT, 0.95) * (branch == "branch1b"))]
traj_specs += [("noise_traj", n_traj - sum(n for _, n, _ in traj_specs), np.zeros(CONFIG["n_cells"]))]

traj_names = [f"TRAJ_{p}_{j}" for p, n, _ in traj_specs for j in range(n)]
traj_profiles = np.vstack([np.tile(prof, (n, 1)) for _, n, prof in traj_specs])

lig_names = [l for l, _ in LR_PAIRS]
rec_names = [r for _, r in LR_PAIRS]
roles = ROLES[scenario]
amp = CONFIG["effect_amp"] * rng.normal(1.0, 0.15, N_PAIRS)
post_all = ramp(pt, FORK1_PT, 0.85) * on("branch0", "branch1", "branch1_shared", "branch1a", "branch1b")
focal_sel = (branch == FOCAL_BRANCH_NESTED if spec["nested"] else branch == FOCAL_BRANCH_2BR)
nested_sel = branch == FOCAL_BRANCH_NESTED
fork_pt = FORK2_PT if spec["nested"] else FORK1_PT
prof_by_role = np.vstack([
    ramp(pt, FORK1_PT, 0.85) * focal_sel,
    ramp(pt, FORK2_PT, 0.95) * nested_sel,
    post_all,
    np.zeros(CONFIG["n_cells"]),
    ramp(pt, FORK1_PT, PULSE_PLATEAU_START) * (1.0 - ramp(pt, PULSE_PLATEAU_END, PULSE_END)) * focal_sel,
])
role_id = np.array([{"focal": 0, "nested": 1, "shared": 2, "flat": 3,
                     "transient": 4}[r] for r in roles])
role_prof = amp[:, None] * prof_by_role[role_id]
lr_side_lig = np.where(is_src[None, :], role_prof, 0.0)
lr_side_rec = np.where(is_tgt[None, :], role_prof, 0.0)

gene_names = traj_names + lig_names + rec_names + [f"NOISE_{j}" for j in range(n_noise)]
gene_class = ([n for n, cnt, _ in traj_specs for _ in range(cnt)] +
              ["lr_ligand"] * N_PAIRS + ["lr_receptor"] * N_PAIRS + ["noise"] * n_noise)

log_base = np.concatenate([
    rng.uniform(np.log(3.0), np.log(8.0), n_traj),
    rng.uniform(np.log(2.0), np.log(5.0), 2 * N_PAIRS),
    rng.uniform(np.log(0.5), np.log(3.0), n_noise),
])
log_mu = log_base[:, None] + np.vstack([traj_profiles, lr_side_lig, lr_side_rec,
                                        np.zeros((n_noise, CONFIG["n_cells"]))])
donor_eff = rng.normal(0.0, 0.12, CONFIG["n_donors"])
log_mu += donor_eff[donor][None, :]
conf_extra = np.zeros_like(log_mu)
conf_extra[n_traj:n_traj + 2 * N_PAIRS, :] = (
    spec["donor_confound"] * is_confounded_donor[None, :] * np.log(3.0) *
    np.r_[np.full(N_PAIRS, float(scenario != "donor_receptor_only")),
          np.full(N_PAIRS, float(scenario != "donor_ligand_only"))][:, None]
)
log_mu += conf_extra
mu = np.exp(np.clip(log_mu, -8.0, np.log(200.0)))

# NB counts via Gamma-Poisson. The historical dropout stress below RETAINS
# counts with exp(-strength * lambda): high expression is more often erased.
# E[Yobs|mu] = mu * (1 + strength*mu/nb_theta)^(-nb_theta-1), which decreases
# above mu=1/strength. Preserve this legacy input, but do not describe it as
# ordinary sparse sampling or infer generic dropout robustness from its sign.
# sampling_thinning below is the existing direction-preserving count stress.
lam = rng.gamma(CONFIG["nb_theta"], mu / CONFIG["nb_theta"])
counts = rng.poisson(lam).astype(np.float32)
keep = rng.random(mu.shape) < np.exp(-lam * spec["dropout"])
counts[n_traj:n_traj + 2 * N_PAIRS, :] = np.where(
    keep[n_traj:n_traj + 2 * N_PAIRS, :], counts[n_traj:n_traj + 2 * N_PAIRS, :], 0.0
)

# Sampling thinning: uniform binomial Yobs ~ Binomial(Y, q) on ALL genes.
# q is the retention probability (not "1-q dropout"); the direction relation
# E[Yobs] = q * E[Y] is preserved while counts get sparser.
if spec["thin_q"] is not None:
    counts = rng.binomial(
        counts.astype(np.int64), spec["thin_q"]
    ).astype(np.float32)

# ---------------------------------------------------------------------------
# Assemble canonical h5ad.
# ---------------------------------------------------------------------------
lib = counts.sum(axis=0)
lognorm = np.log1p(counts / np.maximum(lib, 1.0) * 1e4).astype(np.float32)

umap_x = pt.copy()
umap_y = np.where(branch == "branch0", 1.6 * ramp(pt, FORK1_PT, 1.0),
          np.where(on("branch1", "branch1_shared"), -1.6 * ramp(pt, FORK1_PT, 0.75),
          np.where(branch == "branch1a", -0.55 - 1.0 * ramp(pt, FORK2_PT, 1.0),
          np.where(branch == "branch1b", -2.6 - 1.0 * ramp(pt, FORK2_PT, 1.0), 0.0))))
umap = np.column_stack([umap_x, umap_y]) + rng.normal(0, 0.04, (CONFIG["n_cells"], 2))

bif_level = np.select(
    [branch == "stem", on("branch0", "branch1", "branch1_shared"), on("branch1a", "branch1b")],
    [0, 1, 2],
)
obs = pd.DataFrame({
    "true_pseudotime": pt,
    "true_branch": pd.Categorical(branch),
    "true_bifurcation_level": bif_level.astype(np.int64),
    "clusters": pd.Categorical(cell_type),
    "cell_type": pd.Categorical(cell_type),
    "batch": pd.Categorical(donor_label),
    "donor": pd.Categorical(donor_label),
    "cell_id": cell_id,
}, index=pd.Index(cell_id, name="cell_barcode"))

var = pd.DataFrame({
    "gene_class": pd.Categorical(gene_class),
    "highly_variable": np.isin(gene_class, ["early", "shared", "branch0", "branch1",
                                            "branch1_shared", "branch1a", "branch1b"]),
    "means": mu.mean(axis=1),
    "dispersions": np.full(n_genes, np.nan),
    "dispersions_norm": np.full(n_genes, np.nan),
}, index=pd.Index(gene_names, name="gene_symbol"))

adata = ad.AnnData(
    X=sparse.csr_matrix(counts.T),
    obs=obs,
    var=var,
)
adata.layers["counts"] = sparse.csr_matrix(counts.T)
adata.layers["lognorm"] = sparse.csr_matrix(lognorm.T)
adata.obsm["X_umap"] = umap.astype(np.float32)
adata.obsm["X_pca"] = umap.astype(np.float32)
adata.uns["log1p"] = {}
adata.uns["hvg"] = {}
adata.uns["sim_scenario"] = scenario
adata.uns["fork_pt"] = float(fork_pt)

out_dir = Path(CONFIG["paths"]["out_dir"])
out_dir.mkdir(parents=True, exist_ok=True)
ds_name = (f"lrsim_{scenario}" if CONFIG["seed"] == 12345
           else f"lrsim_{scenario}_s{CONFIG['seed']}")
adata.write_h5ad(out_dir / f"{ds_name}_canonical.h5ad")

truth = pd.DataFrame({
    "dataset": ds_name,
    "scenario": scenario,
    "pair_idx": np.arange(1, N_PAIRS + 1),
    "ligand": lig_names,
    "receptor": rec_names,
    "role": roles,
    "true_positive": np.isin(roles, ["focal", "nested", "transient"]),
    "focal_branch": np.where(np.array(roles) == "nested", FOCAL_BRANCH_NESTED, FOCAL_BRANCH_2BR),
    "other_branch": np.where(np.array(roles) == "nested", "branch1b", "branch1"),
    "fork_pt": np.where(np.array(roles) == "nested", FORK2_PT, FORK1_PT),
    "donor_confounded": spec["donor_confound"],
    "amp_log": amp,
})
truth.to_csv(out_dir / f"{ds_name}_truth.csv", index=False, encoding="utf-8", lineterminator="\n")

cells = pd.DataFrame({
    "cell_id": cell_id,
    "true_pseudotime": pt,
    "true_branch": branch,
    "cell_type": cell_type,
    "donor": donor_label,
})
cells.to_csv(out_dir / f"{ds_name}_cells.csv", index=False,
             encoding="utf-8", lineterminator="\n")

manifest = {
    "dataset": ds_name, "scenario": scenario, "seed": int(CONFIG["seed"]),
    "n_cells": CONFIG["n_cells"], "n_genes": int(n_genes), "n_donors": CONFIG["n_donors"],
    "fork_pt": float(fork_pt), "nested": spec["nested"],
    "donor_confounded": spec["donor_confound"], "dropout_strength": spec["dropout"],
    "thin_q": spec["thin_q"],
    "segments": [list(map(float, s[1:3])) + [s[0], float(s[3])] for s in segments],
    "species_for_lr_db": "human",
}
(out_dir / f"{ds_name}_manifest.json").write_text(
    json.dumps(manifest, indent=2) + "\n", encoding="utf-8", newline="\n")

print(f"{ds_name}: {adata.shape[0]} cells x {adata.shape[1]} genes, "
      f"pos={int(truth['true_positive'].sum())}/{N_PAIRS}, "
      f"counts_mean={counts.mean():.3f}, lr_nonzero={float((counts[n_traj:n_traj+2*N_PAIRS] > 0).mean()):.3f}")
