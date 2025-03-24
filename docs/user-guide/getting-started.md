# Overview

The scRNA github repository is stored locally, and will be used for project deployment. Multiple projects can be deployed from this one point simultaneously, without concern.

## 1. Getting Started

## 1.1 Introduction

The scRNA Pipelie beings at various stages, depending on the users needs. The pipeline can begin with GEX FASTQ files, performing cell counting, with 10X Genomic's [CellRanger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger). Then, normalization and pre-processing occurs, using custom [R](https://www.r-project.org/) scripts with packages like [Seurat](https://satijalab.org/seurat/). Alternatively the user can begin with .h5 files, beginning the pipeline post-preprocessing. [TODO: add as we go!]

## 1.2 Setup Dependencies

scRNA has several dependencies listed below. These dependencies can be installed by a sysadmin. All dependencies will be automatically loaded if running from Biowulf.

- nextflow: "nextflow/23.04.1"
- cellranger "cellranger:7.1.0"
- R: "R/4.3"

Docker containers to run the pipeline are currently in development.

## 1.3 Login to the cluster

scRNA has been exclusively tested on Biowulf HPC. Login to the cluster's head node and move into the pipeline location.

```
# ssh into cluster's head node
ssh -Y $USER@biowulf.nih.gov
```

## 1.4 Load an interactive session

An interactive session should be started before performing any of the pipeline sub-commands, even if the pipeline is to be executed on the cluster.

```
# Grab an interactive node
sinteractive --time=12:00:00 --mem=8gb  --cpus-per-task=4 --pty bash
```
