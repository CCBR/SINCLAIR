# 3. Running the Pipeline

## 3.1 Pipeline Overview

The Nextflow workflow can be run as follows:

```
nextflow run main.nf \
    -entry $datatype \
    -profile biowulf \
    --input assets/input_manifest.csv \
    --contrast assets/contrast_manifest.csv \
    --outdir /data/sevillas2/scRNA_test \
    --species $species \
    $args
```

## 3.2 Commands explained

The following explains each of the command options:

- entry: accepts the datatype to be used; IE gex
- profile: how to run the processes; IE biowulf singularity, docker
- input: input_manifest.csv location
- contrast: contrast_manifest.csv location
- outdir: complete path to the output dir
- species: species to be used
- run_cellranger: whether or not to run cellranger on dataset; IE Y, N
- args: any additional arguments; IE --stub-run

## 3.3 Typical Workflow

A typical command workflow, running the pipeline for the first time locally, is as follows:

```
nextflow run main.nf \
    -entry gex \
    -profile biowulf \
    --input assets/input_manifest.csv \
    --contrast assets/contrast_manifest.csv \
    --outdir /path/to/scRNA_test \
    --species hg19
```

A typical command workflow, running the pipeline for a repeated time locally, running cellranger, is as follows:

```
nextflow run main.nf -resume \
    -entry gex \
    -profile biowulf \
    --run_cellranger Y \
    --input assets/input_manifest.csv \
    --contrast assets/contrast_manifest.csv \
    --outdir /path/to/scRNA_test \
    --run_cellranger Y \
    --species hg19
```

A typical command workflow, running the pipeline in a `dryrun mode`, without running cellranger, is as follows:

```
nextflow run main.nf \
    -entry gex \
    -profile biowulf \
    --run_cellranger Y \
    --input assets/input_manifest.csv \
    --contrast assets/contrast_manifest.csv \
    --outdir /path/to/scRNA_test \
    --species hg19 \
    --run_cellranger N \
    --stub-run
```

Alternatively a script was created to run the pipeline, which takes the following flags:

- species: hg19
- datatype: GEX
- outDir: /path/to/output/dir
- resume: Y or N
- run_cellranger: Y or N
- stubrun: Y or N

Examples:

```
# run GEX on test data, for the first time
sh run_scRNA.sh hg19 GEX N /path/to/output/dir

# first pass, with and without cellranger
sh run_scRNA.sh hg19 GEX /path/to/output/dir N Y N
sh run_scRNA.sh hg19 GEX /path/to/output/dir N Y N

# resume, with and without cellranger
sh run_scRNA.sh hg19 GEX /path/to/output/dir Y Y N
sh run_scRNA.sh hg19 GEX /path/to/output/dir Y N N

# first pass, with cellranger, with and without a dryrun
sh run_scRNA.sh hg19 GEX /path/to/output/dir N Y Y
sh run_scRNA.sh hg19 GEX /path/to/output/dir N Y N
```
