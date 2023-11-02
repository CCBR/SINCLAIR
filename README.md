# SINCLAIR - **SIN**gle **C**el**L** **A**nalys**I**s **R**esource

## Table of contents

- [1. Introduction](#1-Introduction)
- [2. Overview](#2-Overview)
- [3. Documentation](#3-Documentation)
- [4. Feedback](#4-Feedback)

### 1. Introduction

The [SINCLAIR - **SIN**gle **C**el**L** **A**nalys**I**s **R**esource](#sinclair---single-cell-analysis-resource) was developed by the CCR Collaborative Bioinformatics Resource as an open-source, reproducible solution for multiple single cell next-generation modalities. It has been developed solely on [Biowulf](https://hpc.nih.gov/) using [Nextflow](https://www.nextflow.io/).

### 2. Overview

The pipeline currently begins with either sample FASTQ file or h5 Aligned reads, completing per sample quality control, and per-contrast integration. Quality control reports are generated, as are per-contrast integration reports.

![Single cell RNA-Seq GEX pipeline](./resources/scRNA.svg) <sup>**Overview of Single Cell RNASeq Gene Expression Pipeline**</sup>

### 3. Basic Deployment

```
  USAGE:
    bash sinclair -m/--runmode=<RUNMODE> -w/--workdir=<WORKDIR>
  Required Arguments:
    1.  RUNMODE: [Type: String] Valid options:
      *) init : initialize workdir
      *) run : run with slurm
      *) resume : continue pipeline
      *) stubrun : nextflow stubrun; will run locally with test data
      *) runlocal : run without submitting to sbatch
      *) testrun: run on cluster with test dataset
    2.  WORKDIR: [Type: String]: Absolute or relative path to the output folder with write permissions.
    3.  ENTRY: [Type: String] Valid options:
      *) GEX
  Optional Arguments:
    1.  RESUME: [Type: String:] Valid options: Y, N; default: N
```

Example workflow

```
# 1) run initialization
bash sinclair --runmode=init --workdir=/path/to/output/dir --entry=GEX

# 2) update the config files as needed
## can change whether cellranger is deployed, species, names of manifest files (default locations listed below)
/path/to/output/dirnextflow.config
/path/to/output/dir/assets/contrast_manifest.csv /path/to/output/dir/assets/input_manifest.csv

# 3) deploy the pipeline
## A) STUBRUN
bash sinclair --runmode=stubrun --workdir=/path/to/output/dir --entry=GEX

## B) local run
bash sinclair --runmode=runlocal --workdir=/path/to/output/dir --entry=GEX

## C) submit to slurm
bash sinclair --runmode=run --workdir=/path/to/output/dir --entry=GEX

# 4) OPTIONAL resume
bash sinclair --runmode=run --workdir=/path/to/output/dir --entry=GEX --resume=Y
```

### 4. Detailed Documentation

Please view the repositories [documentation](https://symmetrical-adventure-ovjq9gl.pages.github.io/) for full details on deploying the pipeline, features, testing, and expected outputs.

### 5. Feedback

For comments/suggestions/advice please reach out to [Samantha Chill](mailto:samantha.sevilla@nih.gov) or [Nathan Wong](mailto:nathan.wong@nih.gov).
