# CCBR/SINCLAIR pipeline parameters

SINgle CelL AnalysIs Resource

## Input/output options

Define where the pipeline should find input data and save output data.

| Parameter  | Description                                                                                                                                                                                                                                                                                                                                                                      | Type     | Default                      | Required | Hidden |
| ---------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------- | ---------------------------- | -------- | ------ |
| `input`    | Path to CSV file containing information about the samples in the experiment. <details><summary>Help</summary><small>You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 5 columns, and a header row.</small></details> | `string` | assets/input_manifest.csv    | True     |        |
| `contrast` | Path to CSV file with contrast specification                                                                                                                                                                                                                                                                                                                                     | `string` | assets/contrast_manifest.csv | True     |        |
| `outdir`   | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.                                                                                                                                                                                                                                                         | `string` | output                       | True     |        |

## Main options

| Parameter            | Description                                                                             | Type      | Default                   | Required | Hidden |
| -------------------- | --------------------------------------------------------------------------------------- | --------- | ------------------------- | -------- | ------ |
| `species`            |                                                                                         | `string`  | hg19                      |          |        |
| `run_cellranger`     | whether to run cellranger                                                               | `string`  | Y                         |          |        |
| `vars_to_regress`    |                                                                                         | `string`  | percent.mt,nFeature_RNA   |          |        |
| `qc_filtering`       |                                                                                         | `string`  | manual                    |          |        |
| `nCount_RNA_max`     |                                                                                         | `integer` | 500000                    |          |        |
| `nCount_RNA_min`     |                                                                                         | `integer` | 1000                      |          |        |
| `nFeature_RNA_max`   |                                                                                         | `integer` | 5000                      |          |        |
| `nFeature_RNA_min`   |                                                                                         | `integer` | 200                       |          |        |
| `percent_mt_max`     |                                                                                         | `integer` | 10                        |          |        |
| `percent_mt_min`     |                                                                                         | `integer` | 0                         |          |        |
| `run_doublet_finder` |                                                                                         | `string`  | Y                         |          |        |
| `seurat_resolution`  |                                                                                         | `string`  | 0.1,0.2,0.3,0.5,0.6,0.8,1 |          |        |
| `npcs`               |                                                                                         | `integer` | 50                        |          |        |
| `resolution_list`    |                                                                                         | `string`  | 0.1,0.2,0.3,0.5,0.6,0.8,1 |          |        |
| `genome_dir`         | Path to the genome references. Overridden by platform configs, e.g. conf/biowulf.config | `string`  |                           |          |        |

## Institutional config options

| Parameter                    | Description | Type     | Default | Required | Hidden |
| ---------------------------- | ----------- | -------- | ------- | -------- | ------ |
| `config_profile_name`        |             | `string` |         |          |        |
| `config_profile_description` |             | `string` |         |          |        |
| `config_profile_contact`     |             | `string` |         |          |        |
| `config_profile_url`         |             | `string` |         |          |        |

## Generic options

| Parameter          | Description | Type      | Default                        | Required | Hidden |
| ------------------ | ----------- | --------- | ------------------------------ | -------- | ------ |
| `publish_dir_mode` |             | `string`  | link                           | True     |        |
| `tracedir`         |             | `string`  | ${params.outdir}/pipeline_info | True     |        |
| `max_memory`       |             | `string`  | 128.GB                         | True     |        |
| `max_cpus`         |             | `integer` | 48                             | True     |        |
| `max_time`         |             | `string`  | 240.h                          | True     |        |

## Containers

Docker/Singularity containers to use for processes. Must be available in dockerhub

| Parameter         | Description | Type     | Default                             | Required | Hidden |
| ----------------- | ----------- | -------- | ----------------------------------- | -------- | ------ |
| `base_container`  |             | `string` | nciccbr/ccbr_ubuntu_base_20.04:v6.1 |          | True   |
| `baser_container` |             | `string` | nciccbr/sinclair_baser:0.1.0        |          | True   |

## Hidden options

| Parameter               | Description | Type     | Default                                                              | Required | Hidden |
| ----------------------- | ----------- | -------- | -------------------------------------------------------------------- | -------- | ------ |
| `Rlib_dir`              |             | `string` | /gpfs/gsfs10/users/CCBR_Pipeliner/db/PipeDB/Rlibrary_4.3_scRNA_RHEL8 |          | True   |
| `conda_path`            |             | `string` | /gpfs/gsfs10/users/CCBR_Pipeliner/db/PipeDB/Conda/envs/scvi-env      |          | True   |
| `python_path`           |             | `string` | /gpfs/gsfs10/users/CCBR_Pipeliner/db/PipeDB/Conda/envs/scvi-env/bin  |          | True   |
| `Rpkg`                  |             | `string` | ${projectDir}/conf/Rpack.config                                      |          | True   |
| `script_functions`      |             | `string` | ${projectDir}/bin/scRNA_functions.R                                  |          | True   |
| `script_preprocess`     |             | `string` | ${projectDir}/bin/seurat_preprocess.Rmd                              |          | True   |
| `script_merge`          |             | `string` | ${projectDir}/bin/seurat_merge.Rmd                                   |          | True   |
| `script_bc_harmony`     |             | `string` | ${projectDir}/bin/batch_correction_harmony.Rmd                       |          | True   |
| `script_bc_rpca`        |             | `string` | ${projectDir}/bin/batch_correction_rpca.Rmd                          |          | True   |
| `script_bc_cca`         |             | `string` | ${projectDir}/bin/batch_correction_cca.Rmd                           |          | True   |
| `script_liger`          |             | `string` | ${projectDir}/bin/batch_correction_liger.Rmd                         |          | True   |
| `script_bc_integration` |             | `string` | ${projectDir}/bin/batch_correction_integration.Rmd                   |          | True   |

<!-- this doc is generated by: nf-core pipelines schema docs -->
