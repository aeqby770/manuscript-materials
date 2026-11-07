import os
CONFIG = {'paths': {'root': 'E:/科研项目/LinkR', 'data': 'inputs/external_validation/GSE48968/GSE48968_BMDC_WT.h5ad',
    'pairs': 'results/validation/gse48968/external_validation_time_root/all_lr_pairs.csv',
    'output': 'results/by_algorithm/GSE48968_BMDC_WT', 'vendor': 'inputs/vendor/TraSig'},
    'threads': os.cpu_count(), 'plot': {}, 'tpm_factor': 0.01, 'time_window': 1}
os.environ['OMP_NUM_THREADS'] = str(CONFIG['threads'])
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from scipy.io import mmwrite
import matplotlib.pyplot as plt
plt.rcParams.update(CONFIG['plot'])
root = Path(CONFIG['paths']['root'])
out = root / CONFIG['paths']['output']
out.mkdir(parents=True, exist_ok=True)
data = ad.read_h5ad(root / CONFIG['paths']['data'])
pairs = pd.read_csv(root / CONFIG['paths']['pairs'], encoding='utf-8-sig')
E = data.X.T.tocsc().copy()
E.data = np.log1p(E.data * CONFIG['tpm_factor'])
expression = pd.DataFrame(E.toarray(), index=data.var_names, columns=data.obs_names)
partners = pd.Index(pd.unique(pairs[['ligand_genes', 'receptor_genes']].to_numpy().ravel()), name='partner')
parts = partners.to_series().str.split('_').explode().rename('gene').reset_index()
assert parts.gene.isin(expression.index).all()
# Match production chemistry: average each measured subunit within a population
# before taking the minimum across a complex's subunits.
parts.to_csv(out / 'complex_feature_gene_map.csv', index=False, encoding='utf-8')
gene_pools = expression.T.assign(stimulation=data.obs.stimulation.to_numpy(), time_hours=data.obs.time_hours.to_numpy()).groupby(['time_hours', 'stimulation'], observed=True).mean()
gene_index = expression.index.get_indexer(parts.gene)
partner_profiles = pd.DataFrame(gene_pools.to_numpy()[:, gene_index].T,
    index=pd.Index(parts.partner, name='partner')).groupby(level=0, sort=False).min().reindex(partners).T
partner_profiles.index = gene_pools.index
pools = partner_profiles
times = np.sort(data.obs.time_hours.unique())
conditions = pd.Index(np.sort(data.obs.stimulation.unique()), name='stimulation')
full = pools.reindex(pd.MultiIndex.from_product([times, conditions], names=['time_hours', 'stimulation']))
available = full.notna().all(axis=1).to_numpy().reshape(len(times), len(conditions))
cube = full.fillna(0).to_numpy().reshape(len(times), len(conditions), len(partners))
li = partners.get_indexer(pairs.ligand_genes)
ri = partners.get_indexer(pairs.receptor_genes)
S = cube[:, :, li]
R = cube[:, :, ri]
sys.path.insert(0, str(root / CONFIG['paths']['vendor']))
from trasig.trasig import trasig
native = trasig(['dot'], {0: conditions.tolist()}, dict(zip(conditions, np.zeros(len(conditions), dtype=int))),
    data.obs.time_hours.to_numpy(), data.obs.stimulation.to_numpy(), conditions.to_numpy())
native_dot = native.paired_dot_paths(S, R)
# No random permutations or pseudo p values: this is the published observable dot score.
dot = pd.DataFrame(native_dot.T, index=pd.Index(pairs.pair_idx, name='pair_idx'), columns=conditions).reset_index().melt(
    id_vars='pair_idx', var_name='stimulation', value_name='score')
dot['method'] = 'TraSig_native_dot_same_resource'
dot['time_window_measurements'] = CONFIG['time_window']
product = S * R
max_product = np.max(product, axis=0)
peak = pd.DataFrame(max_product.T, index=pd.Index(pairs.pair_idx, name='pair_idx'), columns=conditions).reset_index().melt(
    id_vars='pair_idx', var_name='stimulation', value_name='score')
peak['method'] = 'maximum_observed_mean_LR_product'
global_gene_mean = expression.mean(axis=1).to_numpy()
pooled = pd.Series(global_gene_mean[gene_index], index=pd.Index(parts.partner, name='partner')).groupby(level=0, sort=False).min().reindex(partners).to_numpy()
static = pairs[['pair_idx']].copy()
static['score'] = pooled[li] * pooled[ri]
static['stimulation'] = 'all_primary_conditions'
static['method'] = 'pooled_mean_LR_product'
dot_global = dot.groupby('pair_idx', as_index=False).agg(score=('score', 'max'))
dot_global['method'] = 'TraSig_native_dot_same_resource'
dot_global['stimulation'] = 'all_primary_conditions'
peak_global = peak.groupby('pair_idx', as_index=False).agg(score=('score', 'max'))
peak_global['method'] = 'maximum_observed_mean_LR_product'
peak_global['stimulation'] = 'all_primary_conditions'
current = pd.read_csv(root / 'results/validation/gse48968/external_validation_time_root/native_pair_ranking.csv', encoding='utf-8')
previous = pd.read_csv(root / 'results/validation/gse48968/external_validation_time_root_legacy/native_pair_ranking.csv', encoding='utf-8')
current['method'] = 'LineaLinkR_G_point_Q_curve'
previous['method'] = 'LineaLinkR_G_point_G_curve'
own = pd.concat([current, previous], ignore_index=True)[['pair_idx', 'score', 'method']]
own['stimulation'] = 'all_primary_conditions'
rankings = pd.concat([dot, peak, static, dot_global, peak_global, own], ignore_index=True).merge(pairs, on='pair_idx', validate='many_to_one')
common = current.loc[current.covered, 'pair_idx'].to_numpy()
common = np.intersect1d(common, previous.loc[previous.covered, 'pair_idx'].to_numpy())
rankings['in_common_evaluable_panel'] = rankings.pair_idx.isin(common)
rankings['score_common'] = rankings.score.where(rankings.in_common_evaluable_panel)
rankings['rank_common'] = rankings.groupby(['method', 'stimulation'], observed=True).score_common.rank(method='min', ascending=False)
rankings['n_common_candidates'] = len(common)
rankings['rank'] = rankings.groupby(['method', 'stimulation'], observed=True).score.rank(method='min', ascending=False)
rankings['n_candidates'] = rankings.groupby(['method', 'stimulation'], observed=True).score.transform('size')
rankings.to_csv(out / 'same_resource_pair_rankings.csv', index=False, encoding='utf-8')
pd.DataFrame({'input_cells': [data.n_obs], 'input_genes': [data.n_vars], 'candidate_pairs': [len(pairs)],
    'native_api': ['trasig.trasig.paired_dot_paths'], 'complex_rule': ['minimum of measured subunit population means'],
    'calendar_time': ['actual hours; no ordering fit'], 'window_measurements': [CONFIG['time_window']],
    'permutations': [0], 'missing_condition_times': [int((~available).sum())],
    'missing_profile_rule': ['TraSig default nan2zero; unobserved conditions retained in availability table']}).to_csv(
    out / 'controlled_comparison_coverage.csv', index=False, encoding='utf-8')
pd.DataFrame(available, index=times, columns=conditions).to_csv(out / 'observed_profile_availability.csv', encoding='utf-8')
prepared = out / 'log_tpm_input'
prepared.mkdir(exist_ok=True)
with (prepared / 'expr.mtx').open('wb') as stream:
    mmwrite(stream, E)
(prepared / 'genes.txt').write_text('\n'.join(data.var_names) + '\n', encoding='utf-8')
(prepared / 'cells.txt').write_text('\n'.join(data.obs_names) + '\n', encoding='utf-8')
data.obs.to_csv(prepared / 'meta.csv', encoding='utf-8')
pairs.to_csv(out / 'all_lr_pairs.csv', index=False, encoding='utf-8')
print(rankings.loc[rankings.ligand_genes.eq('IFNB1') & rankings.receptor_genes.eq('IFNAR1_IFNAR2'),
    ['method', 'stimulation', 'score', 'rank', 'n_candidates']].to_string(index=False))
