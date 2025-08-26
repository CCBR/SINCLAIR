# Expected Outputs

The following folders are created under the working directory containing the user's input files (e.g. fastqs) and where SINCLAIR was initialized in:

- **batch_correct**: contains the various RDS files containing scRNA data that has been batch corrected using the previously listed methods. A .html report is included that compares the

    - A detailed walkthrough for each method can be found via the following links:
      [CCA](https://direct.mit.edu/neco/article/16/12/2639/6880/Canonical-Correlation-Analysis-An-Overview-with), [HARMONY](https://portals.broadinstitute.org/harmony/advanced.html), [scVCI](https://docs.scvi-tools.org/en/0.6.8/tutorials/basic_tutorial.html), [RPCA](https://satijalab.org/seurat/articles/integration_rpca), [LIGER](https://welch-lab.github.io/liger/index.html)

- **cellranger_counts**: contains the h5 files, if `cellranger count` is deployed
- **pipeline_info**: contains execution_reports, execution_trace and pipeline_dag files from NextFlow
- **samplesheets**: contains the manifests used to identify samples, contrasts, and sample:contrast groupings
- **seurat**: contains multiple seurat-generated directories:
    - **preprocess**: contains sample level data; both RDS and PDF files from pre-processing
    - **merge**: contrast grouped, sample level data; both RDS and PDf files of pre-processed, merging

The results directory structure will have the following format, where group1 and group2 are the terms specified in the contrast_manifest.csv file.

```
results
├── batch_correct
│   ├── group1-group2_batch_correction_cca.rds
│   ├── group1-group2_batch_correction_harmony.rds
│   ├── group1-group2_batch_correction_integration.html
│   ├── group1-group2_batch_correction_liger.rds
│   └── group1-group2_batch_correction_rpca.rds
├── cellranger_counts
│   ├── sample1
│   │   └── outs
│   │       └── filtered_feature_bc_matrix.h5
│   ├── ...
│   └── sampleN
│       └── outs
│           └── filtered_feature_bc_matrix.h5
├── samplesheets
│   ├── project_contrast_samplesheet.csv
│   ├── project_gex_samplesheet.csv
│   └── project_groups_samplesheet.csv
└── seurat
    ├── merge
    │   ├── group1-group2_seurat_merged.pdf
    │   └── group1-group2_seurat_merged.rds
    └── preprocess
        ├── sample1_seurat_preprocess.pdf
        ├── sample1_seurat_preprocess.rds
        ├── ...
        ├── sampleN_seurat_preprocess.pdf
        └── sampleN_seurat_preprocess.rds

```

For single-sample analsyses, the output directory will be similar in structure to the above, but missing batch_correct results:

```
results
├── cellranger_counts
│   ├── sample1
│   │   └── outs
│   │       └── filtered_feature_bc_matrix.h5
├── samplesheets
│   ├── project_contrast_samplesheet.csv
│   ├── project_gex_samplesheet.csv
│   └── project_groups_samplesheet.csv
└── seurat
    ├── merge
    │   ├── group1-group2_seurat_merged.pdf
    └── preprocess
        ├── sample1_seurat_preprocess.pdf

```

The following list contains a description of each file and their potential use(s) cases for single-cell analysis.

- `batch_correct/group1-group2_batch_correction_<method>.rds`

    - A Seurat object containing batch-corrected data using a specific integration method.
    - Use case:
        - Input for downstream analyses such as differential gene expression or pathway analysis.
        - Refer to the `batch_correction_integration.html` report to determine which corrected dataset is most appropriate.

- `batch_correct/group1-group2_batch_correction_integration.html`

    - Interactive HTML report summarizing batch correction and integration diagnostics.
    - Includes UMAP plots colored by cell type and silhouette plots comparing integration quality across methods.
    - Use case:
        - Helps identify the best-performing integration method for downstream analyses by visually evaluating technical correction.

- `cellranger_counts/sampleN/outs/filtered_feature_bc_matrix.h5`

    - Filtered gene expression matrix output from **Cell Ranger**.
    - Use case:
        - Input for downstream Seurat-based analysis workflows.

- `samplesheets/project_gex_samplesheet.csv`

    - Sample sheet containing metadata for gene expression analysis (e.g., sample names, input paths).
    - Use case:
        - Defines sample identities and source locations during preprocessing.

- `samplesheets/project_groups_samplesheet.csv`

    - Metadata file that maps samples to experimental or biological groups.
    - Use case:
        - Used during batch correction or merging steps to define group relationships.

- `samplesheets/project_contrast_samplesheet.csv`

    - Specifies contrast group definitions for differential gene expression testing.
    - Use case:
        - Informs statistical comparisons between defined conditions or groups.

- `seurat/merge/group1-group2_seurat_merged.rds`

    - Seurat object combining two or more samples or groups into a single dataset.
    - Use case:
        - Serves as a foundation for integrated downstream analysis before or after batch correction.

- `seurat/preprocess/sampleN_seurat_preprocess.rds`
    - Preprocessed Seurat object for an individual sample, including quality control and normalization.
    - Use case:
        - Used in early-stage analysis or as an input for merging with other samples.
