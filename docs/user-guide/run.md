# 3. Running the Pipeline

## Running with SINCLAIR command:

As of ccbrpipeliner version 8, sinclair can be run with the command:

```
sinclair run [OPTIONS]
```

Various options can be controlled in the command line call and pipeline parameters can be set in the `params.yml` file.
View the full list of pipeline parameters that can be set in `params.yml` [here](../params.md).

The most commonly used options are described below.

_Default values indicated with \*_

### General CLI arguments:

- `--help` Returns the help statement
- `--mode` Determines if the workflow runs on the current system or is submitted as a slurm job
  - `local`\*
  - `slurm`
- `--forceall` Forces all steps of the workflow to be run
- `-profile` Uses pre-defined profiles to determine particular run configurations
  _Note that this parameter is prepended by a single hyphen, not a double hyphen_

  - `biowulf` Uses the configuration optimized for the Biowulf NIH HPC
  - `test` Applies samples and manifests for the test dataset run

### Input and output parameters:

- `--input` Points to the input manifest `.csv` file
  - `./assets/input_manifest_cellranger.csv`\*
  - `./assets/input_manifest.csv`
  - `other/user-defined/manifest.csv`
- `--contrast` Points to the contrast manifest `.csv` file
  - `./assets/contrasts.csv`\*
- `--outdir` Points to the output directory. Can be manually set
  - `./outputs`\*
- `--species` Which species and genome is to be used for reference in alignment (option) cell type annotation
  - `hg19`\*
  - `hg38`
  - `mm10`
- `--run_cellranger` Whether to run CellRanger for alignment. Also indicates which input manifest file to parse
  - `Y`\*
  - `N`

### Seurat parameters:

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

## Examples:

This run will operate on the `slurm` workflow manager, perform CellRanger alignment to the `mm10` mouse genome, and cluster the cells at the specified resolutions:

```
sinclair run --mode slurm --run_cellranger Y --species mm10 --seurat_resolution 0.2,0.4,0.6,0.8,1
```

This run will operate locally, starting from pre-aligned .h5 files generated from CellRanger and take the human `hg38` genome as its cue for downstream cell type annotation, while forcing the run to start from the beginning.

```
sinclair run --mode local --run_cellranger N --species hg38 --forceall
```
