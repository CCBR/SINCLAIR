# 4. Expected Outputs

The following directories are created under the WORKDIR directory:

- batch_correct: contains the various RDS files for batch correction methods ([CCA](https://direct.mit.edu/neco/article/16/12/2639/6880/Canonical-Correlation-Analysis-An-Overview-with), HARMONY, ScVI, RPCA, LIGER) and batch correction report (HTML)
- cellranger_counts: contains the h5 files, if `cellranger count` is deployed
- pipeline_info: contains execution_reports, execution_trace and pipeline_dag files from NextFlow
- samplesheets: contains the manifests used to identify samples, contrasts, and sample:contrast groupings
- seurat: contains multiple seurat-generated directories:
  - preprocess: contains sample level data; both RDS and PDF files from pre-processing
  - merge: contrast grouped, sample level data; both RDS and PDf files of pre-processed, merging

```
─ batch_correct
│   ├── group1-group2_batch_correction_cca.rds
│   ├── group1-group2_batch_correction_harmony.rds
│   ├── group1-group2_batch_correction_integration.html
│   ├── group1-group2_batch_correction_liger.rds
│   ├── group1-group2_batch_correction_rpca.rds
│   ├── group1-group2_batch_correction_scvi.rds
├── cellranger_counts
│   ├── sample1
│   │   └── outs
│   │       └── filtered_feature_bc_matrix.h5
├── pipeline_info
│   ├── execution_report_2023-09-20_12-45-47.html
│   ├── execution_timeline_2023-09-20_12-45-47.html
│   ├── execution_trace_2023-09-20_12-45-47.txt
│   └── pipeline_dag_2023-09-20_12-45-47.svg
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
```
