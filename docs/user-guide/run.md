# 3. Running the Pipeline

## Running the SINCLAIR command

As of ccbrpipeliner version 8, sinclair can be run with the command:

```sh
# initialize the pipeline (only needs to be done once)
sinclair init --output <output_dir>
# run the pipeline
sinclair run --output <output_dir> [OPTIONS]
```

Various options can be controlled in the command line call and pipeline parameters can be set in the `params.yml` file.

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

- `--input` The input manifest `.csv` file
  - `./assets/input_manifest_cellranger.csv`\*
  - `./assets/input_manifest.csv`
  - `other/user-defined/manifest.csv`
- `--contrast` The contrast manifest `.csv` file
  - `./assets/contrasts.csv`\*
- `--outdir` The nextflow results directory inside the pipeline output directory. Can be manually set
  - `./output`\*
- `--species` Which species and genome is to be used for reference in alignment (option) cell type annotation
  - `hg19`\*
  - `hg38`
  - `mm10`
- `--run_cellranger` Whether to run CellRanger for alignment. Also indicates which input manifest file to parse
  - `true`
  - `false`

#### Seurat parameters

<details>
<summary></summary>

- `vars_to_regress` Variables whose effects should be regressed to eliminate potential noise
  - `percent.mt`
  - `nFeature_RNA`
  - `S.Score`
  - `G2M.Score`
  - `nCount_RNA`
- `qc_filtering` Filtering method
  - `miqc`\* Uses the MiQC parameters
  - `manual` Uses the
- `nCount_RNA_max` Maximum number of reads allowed per cell. Cells exceeding the threshold are removed
  - 50000\*
- `nCount_RNA_min` Minimum number of reads allowed per cell. Cells below the threshold are removed
  - 1000\*
- `nFeatures_RNA_max` Maximum number of features (e.g. genes) allowed per cell
  - 5000\*
- `nFeature_RNA_min` Minimum number of features (e.g. genes) allowed per cell
  - 200\*
- `percent_mt_max` Maximum mitochondrial percentage allowed per cell
  - 10\*
- `percent_mt_min` Minimum mitochondrial percentage allowed per cell
  - 0\*
- `run_doublet_finder` Boolean for running the DoubletFinder tool (default T)
- `seurat_resolution` Comma-separated string for resolutions to use when finding unsupervised clusters
  - "0.1,0.2,0.3,0.5,0.6,0.8,1"\*
- `npcs` Number of principal components calculated and used downstream in neighbor-identification, dimensionality reduction (e.g. UMAP/T-SNE), and unsupervised clustering
  - 50\*
  </details>

## Examples

This run will operate on the `slurm` workflow manager, perform CellRanger alignment to the `mm10` mouse genome, and cluster the cells at the specified resolutions:

```
sinclair run --mode slurm --run_cellranger true --species mm10 --seurat_resolution 0.2,0.4,0.6,0.8,1
```

This run will operate locally, starting from pre-aligned .h5 files generated from CellRanger and take the human `hg38` genome as its cue for downstream cell type annotation, while forcing the run to start from the beginning.

```
sinclair run --mode local --run_cellranger false --species hg38 --forceall
```

Specify [pipeline parameters](../params.md) in the `params.yml` file and show a preview of the pipeline run (without actually running it):

```
sinclair run -params-file assets/params.yml -preview
```
