# SINCLAIR Quickstart

## Overview

**SINCLAIR** is an end-to-end NextFlow pipeline for processing single cell RNASeq (scRNASeq) data from either raw fastq.gz files or count files in `.h5` format, as produced by the 10X CellRanger pipeline, and primarily uses [Seurat v5](https://www.satijalab.org/Seurat) as its backbone for downstream processing.

In short, SINCLAIR performs the following functions:

- Alignment from FASTQ (optional)
- Initial quality control and cell filtering per sample
- Sample combination
- Batch correction
- Preliminary cell type annotation
- Preliminary clustering

The final outputs are a set of Seurat `.rds` files that contain all provided samples with and without batch correction, with the latter evaluated using different algorithms (HARMONY, scVI, RPCA, LIGER).

## Initialization

Sinclair can be run either on SLURM or locally. The default — and in most cases, preferable method is running Sinclair on SLURM. To start an interactive SLURM session, a simple sinteractive command is sufficient:

```
sinteractive
```

To run Sinclair locally, start an interactive session with a minimum of 64GB memory, 16 CPUs, 128GB of local storage space, and 8 hours wall-time:

```
sinteractive --mem=64g --cpus-per-task=16 --time=8:00:00 --gres=lscratch:128
```

As of CCBR Pipeliner release 8, instantiate Pipeliner as a module:

```
module load ccbrpipeliner
```

Navigate to your working directory and initialize SINCLAIR:

```
sinclair init
```

## Setting up input files

All input files should follow nomenclature as if generated via [CellRanger](https://www.10xgenomics.com/support/jp/software/cell-ranger/8.0/tutorials/inputs/cr-specifying-fastqs). When starting from fastq files, each sample should have its own directory containing at least `R1` and `R2` (i.e. forward and reverse reads). Additional files that may be included include `I1` index files and reads from multiple lanes. Example minimum data structure for two samples:

```
`/path/to/sample1/sample1_S1_L0001_R1_001.fastq.gz`
`/path/to/sample1/sample1_S1_L0001_R2_001.fastq.gz`

`/path/to/sample2/sample2_S1_L0001_R1_001.fastq.gz`
`/path/to/sample2/sample2_S1_L0001_R2_001.fastq.gz`
```

When starting from `.h5` files that are generated from CellRanger alignment, the directory structure is simpler:

`/path/to/sample1/outputs/filtered_feature_bc_matrix.h5`
`/path/to/sample2/outputs/filtered_feature_bc_matrix.h5`

The `.h5` matrix files should be indicated as `filtered_feature_bc_matrix.h5`, with the sample name indicated in the directory path.

### Converting to .h5 from Matrix Market Exchange (Mtx) Triplet Files

Older versions of CellRanger, as well as some other workflows (e.g. DropSeq, Smart-Seq, PipSeq, etc.) may output matrix.mtx, features.tsv, and barcodes.tsv files. For these instances, refer to [Addendum: Starting from matrix files](preparing-files.md#addendum-starting-from-matrix-files) section for a tutorial on how to convert these files into .h5.

A brief description of each file is:

Matrix Market Exchange (.mtx) files can also be used to store sparse matrices of gene expression data. For each non-zero entry, its row index, column index, and value are stored within a mtx file.

For example, a matrix containing 20000 genes across 5000 barcoded cells, with 10000 non-zero entries can be represented in a mtx file like:

**Matrix dimensions:** 20000 rows × 5000 columns
**Non-zero entries:** 10000

| Row | Column | Value |
| --- | ------ | ----- |
| 1   | 1      | 3     |
| 1   | 2      | 5     |
| 2   | 1      | 1     |

The features.tsv file typically represents the genes (or other molecule types) represented in the rows of the matrix.mtx file. It is a tab-deliminated file with three columns (Feature ID, Feature Name, and Feature Type), like:

| Feature ID      | Gene  | Feature Type    |
| --------------- | ----- | --------------- |
| ENSG00000121410 | CD3D  | Gene Expression |
| ENSG00000198786 | MS4A1 | Gene Expression |
| ENSG00000160791 | LYZ   | Gene Expression |

The barcodes.tsv file list the cell barcodes, eaching representing a droplet or cell that was detected during single-cell sequencing. It is typically a single column, plain text file with one barcode per line:

| Barcode            |
| ------------------ |
| AAACCTGAGATAGGAA-1 |
| AAACCTGAGATAGGGA-1 |
| AAACCTGAGATAGGTT-1 |

## Setting Up Manifest Files

Manifest files are comma-separated variable (.csv) files in the `assets` folder of the SINCLAIR working directory. These contain the filepaths for the input sample files and the contrasts to be included as group identities for downstream differential expression.

Two options exist for sample inputs: `input_manifest.csv` and `input_manifest_cellranger.csv`. Usage depends on the entry point, i.e. whether the samples have already been aligned using CellRanger.

### input_manifest.csv

This file is used when starting from FASTQ files and requires alignment via CellRanger. A .csv template is available in the assets directory and can be modified as shown below:

| `masterID` | `uniqueID`       | `groupID` | `dataType` | `input_dir`      |
| ---------- | ---------------- | --------- | ---------- | ---------------- |
| parentID_1 | uniqueSampleID_1 | groupID_1 | gex        | path/to/fastqs/1 |
| parentID_2 | uniqueSampleID_2 | groupID_2 | gex        | path/to/fastqs/2 |

The `masterID `column can be used to indicate if samples are replicates of the same sample and are not required to be unique. The `groupID` indicates the contrast group for each sample. The `input_dir` must point to a series of fastq files as generated by the 10X Chromium pipeline, or fastq files that follow the [same naming convention](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-specifying-fastqs#file-naming-convention).

### input_manifest_cellranger.csv

This file is used if the alignment and read counting has already been performed, such as through 10X CellRanger software or a similar tool. The input type is expected to be in `.h5` format. Tools will be made available to convert `mtx` triplet files into `.h5`

| masterID     | uniqueID           | groupID     | dataType | input_dir            |
| ------------ | ------------------ | ----------- | -------- | -------------------- |
| `parentID_1` | `uniqueSampleID_1` | `groupID_1` | `gex`    | `path/to/h5Counts/1` |
| `parentID_2` | `uniqueSampleID_2` | `groupID_2` | `gex`    | `path/to/h5Counts/2` |

The primary difference here is that the input_dir points to the directory containing `.h5` file rather than uncounted `.fastq.gz` files.

### contrast_manifest.csv

This file contains the comparisons to be generated, and will indicate samples that should be included in different combinations. For each contrast indicated, only the samples within the specified groups will be processed and combined.

| contrast1 | contrast2 | contrast3 |
| --------- | --------- | --------- |
| group1    | group2    |
| group1    | group2    | group3    |

As many contrasts as there exists groups can be included, with a minimum of 2 groups, as specified both in the `input_manifest.csv`/`input_manifest_cellranger.csv` and the `contrast_manifest.csv` files. If running SINCLAIR on a single sample, the above can be formatted as:

| contrast1 |
| --------- |
| group1    |

## Starting the Run

These instructions will start a basic run of SINCLAIR. For more detailed instructions, please refer to [`3. Running the Pipeline.`](./run.md). When running SINCLAIR from CCBRPipeliner, the following commands are used. When running from a GitHub installation, `sinclair` should be replaced with `bin/sinclair`.

To start a local instance with CellRanger alignment (which is also the default setting):

```
sinclair run --mode local --species=<genome> --run_cellranger true
```

To start a slurm run:

```
sinclair run --mode slurm --species <genome> --run_cellranger true
```

By default, the genome is `hg19`; other options include `mm10` and `hg38`. In order to run SINCLAIR without CellRanger alignment, the parameter `--run_cellranger false` needs to be set and SINCLAIR will now look at the `input_manifest_cellranger.csv` manifest.

## Expected Outputs

Upon initialization, the following folders and files will be loaded into your directory:

```
├── assets
│   ├── contrast_manifest.csv
│   ├── img
│   │   ├── scRNA.jpeg
│   │   └── scRNA_process.jpeg
│   ├── input_manifest_cellranger.csv
│   ├── input_manifest.csv
│   ├── params.yml
│   ├── README.md
│   ├── schema_input.json
│   ├── slurm_header_biowulf.sh
│   └── slurm_header_frce.sh
├── conf
│   ├── base.config
│   ├── base_stub.config
│   ├── biowulf.config
│   ├── ci_stub.config
│   ├── debug.config
│   ├── frce.config
│   ├── interactive.config
│   ├── modules.config
│   ├── Rpack.config
│   ├── slurm.config
│   ├── test.config
│   ├── test_h5.config
│   └── test_pbmc.config
├── log
└── nextflow.config
```

- `assets/` – Contains input files, images, configuration templates, and documentation.

- `conf/` – Nextflow configuration files for different environments and use cases.

- `log/` – Directory where pipeline logs and trace files are stored.

- `nextflow.config` – Main configuration file for the Nextflow pipeline.

During execution, the SINCLAIR workflow stores all temporary outputs in the `work` directory. This directory also supports workflow recovery: if the run fails, intermediate files in work allow the pipeline to resume from the point of failure when the user re-runs the pipeline.

Final results are saved in the `results` directory unless a different output directory was specified in the parameters. The `results` directory will contain 4 subdirectories:

- `batch_correct` contains the combined Seurat `.rds` files for each of the contrasts, with a separate file for each batch correction method, as well as a summary `.html` file.
- `cellranger_counts` contains the `.h5` counts files for each sample produced by the CellRanger software.
- `samplesheets` contains the parsed sample sheets based on the manifest files, as interpreted by NextFlow and SINCLAIR.
- `seurat` contains two subdirectories:
    > - `merge` contains the combined sample Seurat `.rds` files for each set of contrasts prior to batch correction (which can otherwise be referred to as the "uncorrected" object).
    > - `preprocess` contains the individual sample `.rds` files.

When proceeding to downstream secondary analysis, such as [differential expression](./differentialExpression.md), please utilize the `batch_correction_integration.html` files to determine which batch correction method, or even lack thereof, best fits the data. The appropriate file can then be analyzed in R through the Seurat workflow.

For multi-sample analsyses, the output directory will have the following structure:

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
