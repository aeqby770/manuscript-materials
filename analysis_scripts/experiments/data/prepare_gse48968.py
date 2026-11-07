import os
from pathlib import Path
CONFIG = {'paths': {'root': os.environ.get("LLR_WORK_ROOT", str(Path(__file__).resolve().parents[2])),
    'input': 'inputs/external_validation/GSE48968',
    'matrix_a': 'GSE48968_allgenesTPM_GSM1189042_GSM1190902.txt.gz',
    'matrix_b': 'GSE48968_allgenesTPM_GSM1406531_GSM1407094.txt.gz',
    'series': 'GSE48968_family.soft.gz',
    'metrics': 'author_sequencing_metrics.csv',
    'metadata': 'GSE48968_all_library_metadata.csv',
    'training': 'GSE48968_BMDC_WT.h5ad',
    'roots': 'GSE48968_BMDC_WT_root_cells.txt'},
    'threads': os.cpu_count(), 'plot': {}}
os.environ['OMP_NUM_THREADS'] = str(CONFIG['threads'])
import gzip
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import matplotlib.pyplot as plt
plt.rcParams.update(CONFIG['plot'])
root = Path(CONFIG['paths']['root'])
dest = root / CONFIG['paths']['input']
a = pd.read_csv(dest / CONFIG['paths']['matrix_a'], sep='\t', index_col=0, encoding='utf-8-sig')
b = pd.read_csv(dest / CONFIG['paths']['matrix_b'], sep='\t', index_col=0, encoding='utf-8-sig')
assert a.index.is_unique and b.index.is_unique and a.columns.is_unique and b.columns.is_unique
assert a.columns.intersection(b.columns).empty
with gzip.open(dest / CONFIG['paths']['series'], 'rt', encoding='utf-8-sig') as stream:
    lines = pd.Series(stream.read().splitlines(), name='line')
gsm = lines.str.extract(r'^\^SAMPLE = (GSM[0-9]+)$', expand=False).ffill()
fields = lines.str.extract(r'^!Sample_(title|source_name_ch1|description|geo_accession) = (.*)$')
fields.columns = ['field', 'value']
fields['gsm'] = gsm
fields = fields.loc[fields.field.notna()].copy()
assert not fields[['gsm', 'field']].duplicated().any()
geo = fields.pivot(index='gsm', columns='field', values='value')
assert geo.index.equals(pd.Index(geo.geo_accession, name='gsm'))
geo['gsm_number'] = geo.index.str.replace('GSM', '', regex=False).astype(int)
geo_a = geo.loc[geo.gsm_number.between(1189042, 1190902)].copy()
geo_b = geo.loc[geo.gsm_number.between(1406531, 1407094)].sort_values('gsm_number').copy()
assert geo_a.title.is_unique and len(geo_a) == len(a.columns)
assert len(geo_b) == len(b.columns)
ma = geo_a.reset_index().set_index('title').loc[a.columns].reset_index().rename(columns={'index': 'title'})
ma['raw_cell_id'] = a.columns.to_numpy()
ma['matrix_source'] = 'A'
raw_b_well = b.columns.to_series(index=geo_b.index).str.extract(r'_S([0-9]+)_rsem$', expand=False)
geo_b_well = geo_b.title.str.extract(r'_S([0-9]+)$', expand=False)
assert raw_b_well.notna().all() and geo_b_well.notna().all()
assert raw_b_well.equals(geo_b_well)
mb = geo_b.reset_index()
mb['raw_cell_id'] = b.columns.to_numpy()
mb['matrix_source'] = 'B'
mb['raw_group'] = mb.raw_cell_id.str.replace(r'_S[0-9]+_rsem$', '', regex=True)
mb['geo_group'] = mb.title.str.replace(r'_S[0-9]+$', '', regex=True)
assert mb.groupby('raw_group').geo_group.nunique().eq(1).all()
assert mb.groupby('geo_group').raw_group.nunique().eq(1).all()
metrics = pd.read_csv(dest / CONFIG['paths']['metrics'], encoding='utf-8-sig')
qc_names = metrics.Experiment.str.replace(r'_rsem$', '', regex=True)
meta = pd.concat([ma, mb], ignore_index=True)
assert len(meta) == len(a.columns) + len(b.columns) and meta.gsm.is_unique
meta['author_quality_passed'] = meta.title.isin(qc_names)
meta['experiment_id'] = meta.title.str.replace(r'_S[0-9]+$', '', regex=True)
meta['single_cell'] = meta.title.str.contains(r'_S[0-9]+$', regex=True)
meta['biological_replicate_flag'] = meta.experiment_id.isin(['Unstimulated_Replicate', 'LPS_4h_Replicate'])
meta['technical_replicate_flag'] = meta.experiment_id.str.contains('Technical_Replicate', regex=False)
meta['primary_wt_timecourse'] = meta.experiment_id.str.fullmatch(r'Unstimulated|(?:LPS|PAM|PIC)_[1246]h')
meta['model_input'] = meta.matrix_source.eq('A') & meta.author_quality_passed & meta.single_cell & meta.primary_wt_timecourse
meta['stimulation'] = meta.experiment_id.str.extract(r'(LPS|PAM|PIC|IFNB)', expand=False).fillna('none')
meta['time_hours'] = pd.to_numeric(meta.experiment_id.str.extract(r'_([0-9]+)h(?:_|$)', expand=False))
meta.loc[meta.experiment_id.str.startswith('Unstimulated'), 'time_hours'] = 0
meta.loc[meta.experiment_id.eq('On_Chip_Unstim'), 'time_hours'] = 0
meta['experimental_role'] = np.select(
    [~meta.single_cell, ~meta.author_quality_passed, meta.model_input,
     meta.biological_replicate_flag, meta.technical_replicate_flag,
     meta.experiment_id.str.contains('GolgiPlug', regex=False),
     meta.experiment_id.str.contains('_KO_', regex=False),
     meta.experiment_id.str.startswith('On_Chip'),
     meta.experiment_id.str.startswith('IFNB')],
    ['population_control', 'author_quality_excluded', 'primary_wt_timecourse',
     'held_out_biological_replicate', 'held_out_technical_replicate',
     'held_out_secretion_block', 'held_out_genetic_perturbation',
     'held_out_isolation', 'held_out_cytokine_stimulation'],
    default='held_out_other_perturbation')
meta.to_csv(dest / CONFIG['paths']['metadata'], index=False, encoding='utf-8')
fit = meta.loc[meta.model_input].copy()
fit.index = pd.Index(fit.gsm.to_numpy(), name='barcode_id')
fit['cell_type'] = 'BMDC'
fit['clusters'] = 'BMDC'
fit['experiment_id'] = fit.experiment_id.astype(str)
expression = a.loc[:, fit.raw_cell_id].to_numpy(dtype=np.float64).T
assert np.isfinite(expression).all() and (expression >= 0).all()
assert fit.time_hours.notna().all() and fit.index.is_unique
assert fit.experiment_id.nunique() == 13
obs = fit.loc[:, ['cell_type', 'clusters', 'experiment_id', 'stimulation', 'time_hours',
                  'raw_cell_id', 'title', 'description']].copy()
var = pd.DataFrame(index=pd.Index(a.index.to_numpy(), name='gene_symbol'))
data = ad.AnnData(sp.csr_matrix(expression), obs=obs, var=var)
data.uns['source'] = 'Shalek et al. Nature 2014; GSE48968 and official Supplementary Table 2'
data.uns['expression_unit'] = 'Author RSEM linear TPM; not UMI counts'
data.uns['fit_scope'] = 'All author-quality-passed single cells from the original WT LPS/PAM/PIC time courses and unstimulated control; no new quality filter or gene filter'
data.uns['held_out'] = 'All genetic, secretion-block, isolation, direct-cytokine, biological-replicate and technical-replicate experiments are withheld from the fitted state graph'
data.uns['replication'] = 'Primary dates each represent one published experimental condition; held-out biological and technical replicates retain their original separate roles'
data.uns['gene_universes'] = 'Matrix A has 27723 genes; matrix B has 10972 genes. Missing genes across matrices are unmeasured, never zero-filled.'
data.write_h5ad(dest / CONFIG['paths']['training'], compression='gzip')
(dest / CONFIG['paths']['roots']).write_text('\n'.join(fit.index[fit.time_hours.eq(0)]) + '\n', encoding='utf-8')
check = ad.read_h5ad(dest / CONFIG['paths']['training'], backed='r')
assert check.obs_names.equals(data.obs_names) and check.var_names.equals(data.var_names)
assert check.shape == data.shape and check.uns['expression_unit'] == data.uns['expression_unit']
check.file.close()
print('TRAINING', data.shape, 'conditions', fit.experiment_id.nunique())
print('B_MAPPING', mb.groupby(['raw_group', 'geo_group']).size().to_dict())
print('ROLES', meta.groupby(['matrix_source', 'experimental_role']).size().to_dict())
print('AUTHOR_QC', meta.author_quality_passed.sum(), 'SINGLE_CELLS', (meta.author_quality_passed & meta.single_cell).sum())
print('TPM_TOTAL_RANGE', expression.sum(axis=1).min(), expression.sum(axis=1).max())
print('TIMECOURSE', fit.groupby(['stimulation', 'time_hours']).size().to_dict())
