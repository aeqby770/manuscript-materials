# LineaLinkR 0.6.4 dataset provenance

`dataset_provenance_064.tsv` is the row-level source register for the 37 datasets marked `in_build_list=TRUE` in `dataset_inventory_audit.csv` and the five external validation or narrative datasets used by the manuscript. It records the dataset identifier, display name, organism, tissue or process, repository accession, original publication, DOI and PMID, official data URL, evidence status, notes, analysis role, source class and evidence path.

## Coverage

| Source class | Rows | Interpretation |
|---|---:|---|
| Independent public biological data | 25 | Distinct biological studies or public cohorts, including five external validation or narrative datasets |
| Derived public object | 1 | `cellrank_stemnet`, derived from the same pancreas experiment as `pancreas` |
| Project-generated simulation | 16 | `simulated_branch` plus 15 canonical simulation conditions |
| **Total** | **42** | 37 active-registry rows plus five external validation or narrative rows |

Twenty-five rows have an original publication verified against PubMed or Crossref. They represent 24 distinct publications because `pancreas` and `cellrank_stemnet` share GSE132188 and the Bastidas-Ponce et al. study. Twenty-four rows resolve to a GEO series or sample URL, one resolves to BioStudies/ArrayExpress, and the Palantir CD34 object resolves to the author-hosted Rep1 file and its HCA project. The simulation rows have no external accession or publication and are identified as project-generated rather than unresolved literature.

## Source hierarchy

Mappings were checked in this order:

1. The exact local H5AD, adjacent metadata, or generator manifest.
2. The official GEO, BioStudies/ArrayExpress, HCA, or author data record.
3. PubMed article metadata and DOI.
4. The upstream package notebook or repository when the stored object is a processed derivative.

Reviews and secondary dataset catalogues were not used in place of the original data paper. File hashes and modification times were not used to infer provenance.

## Corrected publication mapping

The previous `cellagentchat_scc` row cited PMID 32946785 and doi:10.1016/j.cell.2020.08.043. Those identifiers belong to the two-page correction notice. The original Ji et al. article associated with the cutaneous SCC spatial sample is PMID 32579974 and doi:10.1016/j.cell.2020.05.039. GSM4565826 identifies `P6_ST_vis_rep2`, with `Skin` and `tissue: cSCC` in the official GEO sample record.

## Derived and packaged objects

`cellrank_stemnet` is a computational derivative of the embryonic pancreas experiment GSE132188. It is not an independent cohort. The local source notebook loads a processed GSE132188 pancreas object and writes the STEMNET benchmark object.

`human_cd34` is the 5,780-cell Palantir replicate-1 bone-marrow object. Its official author URL is `https://s3.amazonaws.com/dp-lab-data-public/palantir/human_cd34_bm_rep1.h5ad`; the same project is registered in the Human Cell Atlas as `091cf39b-01bc-42e5-9437-f419a66c8a45`. No GEO accession is asserted because the author and HCA records do not assign one to this object.

## Simulations

The simulation evidence consists of the generated H5AD metadata, `benchmark_manifests/simulation_manifest.json`, and the two scripts in `simulation_generators/`. Both generators set seed 42. The dropout variants use deterministic seeds 3042, 5042, 7042 and 9042; the batch variant uses seed 10041.

The branch labels were also read directly from every canonical simulation H5AD. All eight `sim_3br_*` objects contain four distinct `true_branch` values (`0`, `1`, `2`, `3`). The generator's `topology_3branch()` function returns `n_branches=4`, although the outer manifest rows for these files record `n_branches=3`. These rows are marked `project_generated_verified_label_inconsistency`; neither the `3br` identifier nor the manifest value should be treated as evidence for three observed branch labels.

The legacy identifier `simulated_complex_3branch` has the same four-label topology. Its manifest and generator both record four branches. The identifier remains in the active registry for path compatibility, while the provenance row reports the verified four-label content.

## Unresolved item

`zebrafish_neural_dev` is verified as GEO GSE329822 and BioProject PRJNA1461894. The official GEO record describes ventral spinal-cord pMN progenitors at 3, 5 and 7 days post-fertilization but does not link a publication, DOI or PMID as of 2026-07-24. The publication fields therefore remain unresolved. No publication from a similar zebrafish dataset was substituted.

GSE291091 is linked to the Advanced Science article through the matching title, authors and Crossref DOI `10.1002/advs.202507444`. The official GEO record still displays `Citation missing`, so no PMID is asserted.

## Registry boundary

The benchmark boundary remains the active 37-row registry, not earlier directory counts. Five manuscript datasets outside that boundary are retained with explicit roles: GSE274389 and GSE216606 for external validation, GSE143437 and GSE291091 for biological validation, and GSE171524 for the disease narrative. Legacy aliases `smoke_human_cd34_500c`, `cellrank_bone_marrow` and `bonemarrow_scvelo` are excluded because they are a smoke subset or duplicate representations of `human_cd34`.
