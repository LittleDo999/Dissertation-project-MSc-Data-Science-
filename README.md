# CES1, Macrophage Heterogeneity and Cell–Cell Communication in Colorectal Cancer

MSc Data Science dissertation project — Jingyuan Du, University of Bristol.

## Project overview

This project reanalyses published colorectal cancer single-cell RNA-sequencing data to investigate the immune microenvironment from the perspective of CES1 expression and tumour-associated macrophage (TAM) heterogeneity. CellPhoneDB is used to infer candidate cell–cell interactions, with particular attention to C1QC+ and SPP1+ TAM populations and their interactions with T-cell subsets.

The analysis is exploratory. Transcript-based interaction predictions do not establish physical cell contact, protein activity or causal effects of CES1.

## Repository contents

| File | Purpose |
| --- | --- |
| `Zhang_data.ipynb` | Data loading and alignment, log transformation, highly variable gene selection, PCA/UMAP, CES1 and myeloid-marker visualisation, Leiden clustering, tissue-composition summaries and global metadata export. |
| `CellPhoneDB V5.ipynb` | Preparation of CellPhoneDB inputs, Method 2 analysis at subcluster and global-cluster levels, inspection of deconvoluted outputs, and interaction visualisation. |
| `README.md` | Project description, analysis scope and reproduction notes. |

This README documents the supplied snapshots `Zhang_data(1).ipynb` and `CellPhoneDB V5(2).ipynb`, corresponding to the repository names above. Findings are based on code and saved outputs; a fresh end-to-end execution has not been performed.

## Data source and required inputs

The project uses data associated with:

Zhang et al. (2020). *Single-Cell Analyses Inform Mechanisms of Myeloid-Targeted Therapies in Colon Cancer*. **Cell**, 181(2), 442–459.e29. [Publication record](https://pubmed.ncbi.nlm.nih.gov/32302573/).

The upstream notebook reads these human 10x leukocyte files associated with [GEO accession GSE146771](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146771):

- `GSE146771_CRC.Leukocyte.10x.TPM.txt.gz`
- `GSE146771_CRC.Leukocyte.10x.Metadata.txt.gz`

The successful loading block uses a space separator for the expression file and a tab separator for metadata, aligns cells through `CellName`, and produces **43,817 cells × 13,538 genes**. An earlier exploratory block uses an incorrect tab separator for the expression file and returns zero cells; use the corrected loading block.

The code applies `sc.pp.log1p`, stores the pre-subsetting expression object in `adata.raw`, and restricts `adata.X` to **2,000 highly variable genes**. PCA and the global UMAP use this restricted matrix (`n_neighbors=15`, `n_pcs=30`). No additional cell-level QC filtering or batch-correction step is implemented in the supplied code. Published annotations are retained rather than newly assigned for all cell populations.

The CellPhoneDB notebook requires:

- `zhang_processed.h5ad`: an annotated expression object with `Sub_ClusterID` in `adata.obs` and gene identifiers compatible with HGNC symbols.
- `cpdb_inputs/metadata_global.tsv`: a separate global-cluster metadata file used by the global analysis. It is exported in the final non-empty cell of `Zhang_data.ipynb` using `adata.obs_names` and `adata.obs['Global_Cluster']`. Run this export before the global CellPhoneDB analysis.
- `cpdb_db/cellphonedb.zip`: the CellPhoneDB database. The notebook requests database release `v5.0.0` if this file is absent.

The notebook creates `cpdb_inputs/metadata.tsv` and `cpdb_inputs/counts_log1p.h5ad` from the processed object. The latter copies `adata.X`, not `adata.raw.X`. The upstream code log-transforms the supplied TPM-labelled matrix; it does not calculate TPM from raw counts. A subsequent user-side check of the current `zhang_processed.h5ad` reports both `.X` and `.raw` as **43,817 cells × 13,538 genes**. These dimensions do not establish the expression transformation or the identity of the input used for earlier runs.

The uploaded CellPhoneDB notebook contains a saved input-loading output of **43,817 cells × 2,000 genes**, whereas the current processed file checked by the user contains **13,538 genes** in both the main and raw matrices. This discrepancy means that the saved notebook output and current file cannot be treated as a single verified execution state. Check the actual `cpdb_inputs/counts_log1p.h5ad` and its relationship to the result files used in the dissertation before deciding whether a rerun is needed. Even the current exported input may have changed since earlier runs. If a run used only the 2,000-gene input, interactions requiring excluded genes could not be evaluated in that run; this restriction must not be assumed for all existing results. The current processed object does not require recovery from `.raw` merely to regain the 13,538-gene dimension.

The supplied upstream notebook reads `zhang_processed.h5ad` in a later cell but contains no command to create it. Add an explicit export after global preprocessing before attempting a clean end-to-end run. A later write in the CellPhoneDB notebook does not resolve its earlier dependency on that same file.

## Software environment

Recorded in the CellPhoneDB notebook:

- Python **3.8.20**
- CellPhoneDB Python package **5.0.1**
- Requested CellPhoneDB database release **v5.0.0**

The upstream notebook uses a separate Python **3.10.20** environment; its installation output also records scikit-misc **0.5.2** and NumPy **2.2.6**. Do not assume these are the CellPhoneDB environment versions. Other packages used across the workflows include Scanpy, AnnData, pandas, Matplotlib, seaborn, ktplotspy and Leiden clustering dependencies. Jupyter is needed to run the notebooks interactively. Their exact versions are not recorded here; export the original working environment before claiming a fully specified reproducible installation.

The package version and database release are separate version identifiers.

## Analysis workflow

1. Load and align the source expression matrix and metadata; log-transform expression and select highly variable genes for global PCA/UMAP.
2. Examine CES1 and marker expression; inspect **12,625 myeloid cells**. The notebook includes exploratory myeloid embeddings and Leiden clustering at resolution **0.8**, plus a cross-tabulation against the original `Sub_ClusterID` annotations. The saved Leiden output contains **12 clusters (0–11)** despite a plot title referring to 0–9. Tissue-composition plots pool cells within tissue categories rather than estimating patient-level differences.
3. Export the processed object and global metadata before running the CellPhoneDB notebook.
4. Export expression data and subcluster metadata for CellPhoneDB.
5. Run **Method 2**, using `cpdb_statistical_analysis_method.call`, for subclusters.
6. Inspect predicted interaction patterns and compare `hM12` (C1QC+ TAM) and `hM13` (SPP1+ TAM), including interactions with T-cell subsets.
7. Read the Method 2 `deconvoluted` output to inspect participant-level expression information.
8. Run Method 2 with global-cluster metadata and compare selected myeloid, B-cell and T-cell interaction patterns.

### Method 2 parameters

Both recorded calls use:

| Parameter | Value |
| --- | --- |
| `counts_data` | `hgnc_symbol` |
| `iterations` | `1000` |
| `threshold` | `0.1` |
| `pvalue` | `0.05` |
| `debug_seed` | `42` |
| `threads` | `8` |
| `result_precision` | `3` |
| `separator` | `\|` |
| `score_interactions` | `True` |

Saved execution logs show completion of both subcluster and global-cluster Method 2 runs. Most downstream plots filter interactions using `p < 0.05`.

### Method 3 clarification

The inspected notebook **does not implement CellPhoneDB Method 3** (`cpdb_degs_analysis_method`). It contains no DEG input file supplied to that method, no invocation of the method, and no recorded Method 3 results.

The cell commented `# Method 3` reads `statistical_analysis_deconvoluted*.txt`, which is an output of Method 2. This step should be described as **inspection of Method 2 deconvoluted outputs**. It is neither an independent DEG-based interaction analysis nor spatial validation.

See the [official CellPhoneDB methods documentation](https://cellphonedb.readthedocs.io/en/stable/RESULTS-DOCUMENTATION.html).

## Outputs and interpretation

Subcluster results are written to `cpdb_outputs/`; global-cluster results are written to `cpdb_outputs_global/`. Recorded outputs include:

- `statistical_analysis_means_*.txt`
- `statistical_analysis_pvalues_*.txt`
- `statistical_analysis_significant_means_*.txt`
- `statistical_analysis_interaction_scores_*.txt`
- `statistical_analysis_deconvoluted_*.txt`
- `statistical_analysis_deconvoluted_percents_*.txt`

Visualisations include a communication heatmap and selected interaction bubble plots. In several bubble plots, colour represents the CellPhoneDB interaction mean and size represents `-log10(p + 0.001)`. Other plots use different size encodings; consult their plotting code. Interaction means and specificity scores are distinct quantities. Neither should be labelled as a CellChat communication probability.

A significant prediction in one cell pair and a non-significant prediction in another does not, by itself, establish a statistically significant difference between those pairs. These comparisons are descriptive.

## Reproduction status and required preparation

The current notebook retains saved exploratory outputs but requires the following preparation for a clean rerun:

1. Obtain the two source files listed above, use the corrected loading block, and add the missing processed-object export. Remove or fix the early violin call using the non-existent `Detailed Sub Clusters` observation column; the later code uses `Sub_Cluster`.
2. Replace machine-specific absolute paths with paths appropriate to the local checkout.
3. Run the global metadata export at the end of the upstream notebook; ensure cell identifiers match the expression input.
4. Select all output tables from the same intended run. Some cells use hard-coded July output filenames, while others select the latest available files; saved main-analysis logs are from September.
5. Remove or revise the exploratory block referencing `Tumor_cells|hM12` and `Tumor_cells|hM13`: these columns were absent in the recorded output, so that block did not produce the intended comparison.
6. Rename the misleading `# Method 3` comment and change the final plot legend from `Communication Prob. (Mean)` to `CellPhoneDB interaction mean`.
7. Verify which expression input produced the reported communication results. The current processed object has 13,538 genes, but the supplied notebook retains a 2,000-gene input-loading output. Check the exported CellPhoneDB input and match it to the selected run; do not infer historical input dimensions from a current file alone. The upstream code assigns `.raw` after log transformation, so its name does not imply untransformed counts. Rerun and compare results if an input correction is necessary.
8. Restart the kernel and execute in order after resolving these dependencies. Keep the selected input and output versions together.

The CellPhoneDB preparation shown does not filter `Tissue` to tumour-only cells. The upstream data include normal and tumour tissue labels, so these calls should not be described as tumour-only or as separate tumour-versus-normal communication analyses. The analysed immune-cell data do not represent the complete tumour and stromal microenvironment. The results should be interpreted within the cell populations and genes retained in the input.

## Supplementary full-gene analysis

An additional CellPhoneDB Method 2 analysis was completed on 23 September 2026 using 43,817 cells and all 13,538 available genes. Separate runs used global and subcluster annotations.

The supplementary result tables are stored in `supplementary/fullgenes_20260923/`. No additional figures were generated from these tables for the submitted dissertation.
The dissertation's main communication figures retain the original 2,000-highly-variable-gene analysis.

The full expression input is not included in this repository.

## References and acknowledgements

- Zhang et al. (2020), original single-cell study: [publication record](https://pubmed.ncbi.nlm.nih.gov/32302573/).
- CellPhoneDB: [official repository](https://github.com/ventolab/CellphoneDB) and [methods documentation](https://cellphonedb.readthedocs.io/en/stable/RESULTS-DOCUMENTATION.html).

The original study authors are acknowledged for generating the source data. CellChat comparator figures supplied by the supervisor and discussed in the dissertation should be attributed to their source; they are not outputs generated by the CellPhoneDB notebook.
