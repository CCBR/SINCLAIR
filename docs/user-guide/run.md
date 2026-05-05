# Running Pipeline

## Setting Biowulf Interactive Session

Before running SINCLAIR, open the terminal and log in to Biowulf using your NIH credentials. Then, create an interactive session and navigate to your project / working directory. A guide to navigating Biowulf can be found in (https://hpc.nih.gov/docs/userguide.html).

```sh
# login to Biowulf
ssh -Y $USER@biowulf.nih.gov

# create interactive session with desired job specifications
sinteractive --mem=<ram> --cpus-per-task=<cores> --time=<wall_time> --gres=lscratch:<local_scratch_space>

# if using slurm to run sinclair, a simple sinteractive command will suffice
sinteractive --time=<wall_time>

# load ccbrpipeliner containing SINCLAIR and other analytical tools
module load ccbrpipeliner
```

## Running the SINCLAIR command

As of ccbrpipeliner version 8, sinclair can be run with the command:

```sh
# initialize the pipeline (only needs to be done once)
sinclair init --output <output_dir>

# run the pipeline
sinclair run --output <output_dir> [OPTIONS]
```

Various options can be controlled in the command line call and pipeline [parameters](../params.md# CCBR/SINCLAIR pipeline parameters) can be set in the `params.yml` file.

The most commonly used options are described below.

_Default values indicated with \*_

### General CLI arguments

- `--help` Prints the help statement
- `--output` The pipeline output directory (same as the nextflow `launchDir`)
- `--mode` Determines if the workflow runs on the current system or is submitted as a slurm job
  - `local`\*
  - `slurm`
- `--forceall` Forces all steps of the workflow to be run

### Nextflow arguments

_Note that [nextflow arguments](https://www.nextflow.io/docs/latest/reference/cli.html#run) are prepended by a **single hyphen** rather than a double hyphen_

- `-params-file assets/params.yml` Specify the pipeline parameters in a YAML file
- `-profile` Uses pre-defined profiles to determine particular run configurations
  - `test` Applies samples and manifests for the test dataset run
- `-preview` Preview the pipeline without executing it

### Pipeline parameters

These are parameters used within the nextflow workflow.
They can be passed in via the command line or set in the `params.yml` file.
View the full list of pipeline parameters [here](../params.md).

#### Input and output parameters

| Parameter          | Description                                                                        | Example(s) / Default(s)\*                                                                                        |
| ------------------ | ---------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------- |
| `--input`          | Input manifest `.csv` file                                                         | `./assets/input_manifest_cellranger.csv`\*<br>`./assets/input_manifest.csv`<br>`other/user-defined/manifest.csv` |
| `--contrast`       | Contrast manifest `.csv` file                                                      | `./assets/contrasts.csv`\*                                                                                       |
| `--outdir`         | Results directory inside the pipeline output directory                             | `./output`\*                                                                                                     |
| `--species`        | Species/genome reference for alignment (optional: also for cell type annotation)   | `hg19`\*<br>`hg38`<br>`mm10`                                                                                     |
| `--run_cellranger` | Whether to run Cell Ranger for alignment; determines which input manifest to parse | `true`<br>`false`                                                                                                |

#### Seurat parameters

The following is a list containing parameters that can be used for downstream Seurat analysis.

- vars_to_regress – Optional variables to regress out, as a way of eliminating technical noise:

  - percent.mt: Percentage of mitochondrial reads; high values may indicate dead/stressed cells.
  - nFeature_RNA: Number of detected features per cell.
  - S.Score: S-phase cell cycle score.
  - G2M.Score: G2/M-phase cell cycle score.
  - nCount_RNA: Total RNA molecule count per cell.

- qc_filtering – Quality control filtering method:

  - miqc: Uses MiQC parameters.
  - manual: Uses manually specified thresholds.
  - nCount_RNA_max – Maximum total RNA count allowed per cell (cells above this are removed):

    - Default: 50000\*

  - nCount_RNA_min – Minimum total RNA count allowed per cell (cells below this are removed):

    - Default: 1000\*

  - nFeatures_RNA_max – Maximum number of features (e.g., genes) per cell:

    - Default: 5000\*

  - nFeature_RNA_min – Minimum number of features per cell:

    - Default: 200\*

  - percent_mt_max – Maximum mitochondrial read percentage allowed per cell:

    - Default: 10\*

  - percent_mt_min – Minimum mitochondrial read percentage allowed per cell:
    - Default: 0\*

- run_doublet_finder – Boolean flag to run the DoubletFinder tool:

  - Default: true

  - seurat_resolution – Comma-separated list of clustering resolutions for Seurat:

    - Default: "0.1,0.2,0.3,0.5,0.6,0.8,1"\*

  - npcs – Number of principal components used in downstream analyses (e.g., UMAP, clustering):
    - Default: 50\*

## Examples

This run will operate on the `slurm` workflow manager, perform CellRanger alignment to the `mm10` mouse genome, and cluster the cells at the specified resolutions:

```
sinclair run --mode slurm --run_cellranger true --species mm10 --seurat_resolution 0.2,0.4,0.6,0.8,1
```

Users will receive autogenerated emails from slurm@hpc.nih.gov on the status of submitted slurm jobs.

![alt text](image.png)

This run will operate locally, starting from pre-aligned .h5 files generated from CellRanger and take the human `hg38` genome as its cue for downstream cell type annotation, while forcing the run to start from the beginning.

```
sinclair run --mode local --run_cellranger false --species hg38 --forceall
```

Specify [pipeline parameters](../params.md) in the `params.yml` file and show a preview of the pipeline run (without actually running it):

```
sinclair run -params-file assets/params.yml -preview
```
