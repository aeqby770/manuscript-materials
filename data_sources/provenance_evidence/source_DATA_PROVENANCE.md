# Dataset provenance evidence index

The authoritative 0.6.4 dataset register is `../dataset_provenance_064.tsv`; its scope and evidence rules are documented in `../DATA_PROVENANCE_064.md`. This directory holds supporting evidence and is not a separate dataset list.

## Active registry evidence

- Active dataset IDs and canonical paths: `../dataset_inventory_audit.csv`
- Registry construction and species labels: `../../analysis_scripts/config/00_config.R`
- Curated object metadata: `../benchmark_manifests/dataset_manifest.json`
- GEO candidate metadata: `../benchmark_manifests/geo_candidate_manifest.json`
- Current simulation metadata: `../benchmark_manifests/simulation_manifest.json`
- Exact simulation generators: `../simulation_generators/generate_simulated_data.py` and `../simulation_generators/generate_complex_simulated.py`
- CellRank STEMNET derivation: `cellrank_stemnet_source.ipynb`
- Extension-dataset reconstruction metadata: `../derived_dataset_metadata/*_dataset_info.json`

## Verified corrections

The `cellagentchat_scc` source is Ji et al., Cell 2020, PMID 32579974, doi:10.1016/j.cell.2020.05.039. PMID 32946785 and doi:10.1016/j.cell.2020.08.043 identify the correction notice. The exact GEO sample is GSM4565826 (`P6_ST_vis_rep2`, skin, cSCC).

All eight `sim_3br_*` H5AD files contain four `true_branch` values (`0`-`3`). The generator returns `n_branches=4`; the outer manifest rows record 3. `simulated_complex_3branch` also contains four branch labels. The formal TSV marks these rows as label inconsistencies.

## Open record

GSE329822 is the verified source for `zebrafish_neural_dev`, but its official GEO record has no linked publication, DOI or PMID as of 2026-07-24. That field remains unresolved.
