# Advanced Usage

## Loading and Running an Installation from GitHub

### Installation

Clone the repository from GitHub:

```sh
git clone https://github.com/CCBR/SINCLAIR
cd SINCLAIR
```

Install the CLI dependencies in `pyproject.toml` with:

```sh
pip install .
```

Additional required dependencies:

- Nextflow
- Singularity

If you're running sinclair on biowulf, the nextflow and singularity modules will be loaded automatically.

### Preparation and Running

Prepare files as described in the (quickstart)[quickstart.md] and the (preparation)[preparing-files.md] guides.

To start a local instance with CellRanger alignment (which is also the default setting):

```
bin/sinclair run --mode local --species=<genome> --run_cellranger=Y
```

To start a slurm run:

```
bin/sinclair run --mode slurm --species <genome> --run_cellranger Y
```

## Manually adjusting config files

## 2.1 Configs

The configuration files control parameters and software of the pipeline. These files are listed below:

- nextflow.config
- conf/base.config
- conf/modules.config
- conf/Rpack.config

### 2.1.1 NextFlow Config

The configuration file dictates the global information to be used during the pipeline. Users can alter the default values, as needed. View the full list of parameters [here](../params.md).

- input: path to input manifest; example manifests with (`input_manifest_cellranger.csv`) and without (`input_manifest.csv`) `cellranger` are included in assets
- contrast: path to contrast manifest; example manifest (`contrast_manifest.csv`) is included in assets
- outdir: path to output dir
- species: species [options: hg19, mm10]
- run_cellranger: determines whether to run cell ranger; if "Y" is selected, expects FQ inputs, if "N", expects .h5 inputs [options: "Y", "N"]
- vars_to_regress: a comma separated list of any variables to regress during `SCTransform` process; [options: "", "percent.mt,nFeature_RNA,S.Score,G2M.Score,nCount_RNA"]

### 2.1.2 Base Config

The configuration file dictates submission to Biowulf HPC. There are two different ways to control these parameters - first, to control the default settings, and second, to create or edit individual rules. These parameters should be edited with caution, after significant testing.

### 2.1.3 Modules Config

The configuration file dictates process-specific processing parameters, including:

- the version of each software or program that is being used in the pipeline
- output location and file names
- additional arguments to be passed to the process

### 2.1.4 R Package Config

The configuration file dictates which R libraries, and which versions, are loaded into the accompanying R script

## Running via NextFlow commands

### Running via Nextflow

The Nextflow workflow can be also run as follows:

```
nextflow run main.nf \
    --input assets/input_manifest.csv \
    --contrast assets/contrast_manifest.csv \
    -params-file assets/params.yml
```

## 3.2 Commands explained

The following explains each of the command options:

- `-entry`: accepts the datatype to be used; IE gex
- `-profile`: how to run the processes; IE biowulf singularity, docker
- `--input`: input_manifest.csv location
- `--contrast`: contrast_manifest.csv location
- `--species`: species to be used
- `--run_cellranger`: whether or not to run cellranger on dataset; IE Y, N
- args: any additional arguments; IE --stub-run

## 3.3 Typical Workflow

A typical command workflow, running the pipeline for the first time locally, is as follows:

```
nextflow run main.nf \
    \
    -profile biowulf \
    --input assets/input_manifest.csv \
    --contrast assets/contrast_manifest.csv \
    --species hg19
```

A typical command workflow, running the pipeline for a repeated time locally, running cellranger, is as follows:

```
nextflow run main.nf -resume \
    \
    -profile biowulf \
    --run_cellranger Y \
    --input assets/input_manifest.csv \
    --contrast assets/contrast_manifest.csv \
    --run_cellranger Y \
    --species hg19
```

A typical command workflow, running the pipeline in a `dryrun mode`, without running cellranger, is as follows:

```
nextflow run main.nf \
    \
    -profile biowulf \
    --run_cellranger Y \
    --input assets/input_manifest.csv \
    --contrast assets/contrast_manifest.csv \
    --species hg19 \
    --run_cellranger N \
    --stub-run
```
