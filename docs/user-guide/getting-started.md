# Getting Started

The scRNA github repository is stored locally, and will be used for project deployment. Multiple projects can be deployed from this one point simultaneously, without concern.

The scRNA Pipeline begins at various stages, depending on the users needs. The pipeline can begin with GEX FASTQ files, performing cell counting, with 10X Genomics [CellRanger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger). Then, normalization and pre-processing occurs, using custom [R](https://www.r-project.org/) scripts with packages like [Seurat](https://satijalab.org/seurat/). Alternatively the user can begin with .h5 files, beginning the pipeline post-preprocessing.

## Login to the cluster

scRNA has been exclusively tested on Biowulf HPC. Log in to the cluster's head node and move to / create the working directory on which to store the input and outputs files of SINCLAIR.

```
# ssh into cluster's head node
ssh -Y $USER@biowulf.nih.gov
```

```
# (If not done) create working directory
mkdir /path/to/WORKDIR
```

```
# navigate to working directory
cd /path/to/WORKDIR
```

## Load an interactive session

An interactive session should be started before performing any of the pipeline sub-commands, even if the pipeline is to be executed on the cluster.

To start an interactive SLURM session, a simple sinteractive command is sufficient:

```
sinteractive
```

To run SINCLAIR locally, start an interactive session with a minimum of 64GB memory, 16 CPUs, 128GB of local storage space, and 8 hours wall-time:

```
sinteractive --mem=64g --cpus-per-task=16 --time=8:00:00 --gres=lscratch:128
```

## Loading SINCLAIR

The ccbrpipeliner module on Biowulf also loads module dependencies, and should be loaded prior to running SINCLAIR:

```
module load ccbrpipeliner/8
```

Initialize the output directory for SINCLAIR. This will create a new directory where you can copy the necessary files to run the pipeline:

```
sinclair init --output /path/to/output/dir
```

From here, proceed to [preparing the files](./preparing-files.md) and then [running the pipeline](./run.md).
