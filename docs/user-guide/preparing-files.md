# 2. Preparing Files
The pipeline is controlled through editing configuration and manifest files. Defaults are found in the /PIPELINEDIR/conf and /PIPELINEDIR/ directories

![SINCLAIR Process Overview](https://https://github.com/CCBR/SINCLAIR/blob/dev/resources/img/scRNA.jpeg?raw=true) <sup>**Overview of Single Cell RNASeq Gene Expression Process**</sup>

## 2.1 Configs
The configuration files control parameters and software of the pipeline. These files are listed below:

- nextflow.config
- conf/base.config
- conf/modules.config
- conf/process_params.config
- conf/Rpack.config

### 2.1.1 NextFlow Config
The configuration file dictates the global information to be used during the pipeline. 

### 2.1.2 Base Config
The configuration file dictates submission to Biowulf HPC. There are two different ways to control these parameters - first, to control the default settings, and second, to create or edit individual rules. These parameters should be edited with caution, after significant testing.

### 2.1.3 Modules Config
The configuration file dictates process-specific processing parameters, including:

- the version of each software or program that is being used in the pipeline
- output location and file names
- additional arguments to be passed to the process

### 2.1.4 R Package Config
The configuration file dictates which R libraries, and which versions, are loaded into the accompanying R script

### 2.1.3 Process Parameters
The configuration file dictates process-specific user parameters, which varies for each process. Users can choose varied resolution values or QC methods, for example.

## 2.2 Preparing Manifests
There are two manifests, which are required. These files describe information on the samples and desired contrasts. These files are:

- /assets/input_manifest.csv
- /assets/contrast_manifest.csv

### 2.2.1 Input Manifest
This manifest will include information to sample level information. It includes the following column headers:

- masterID: This is the biological sample ID; duplicates are allowed in this column
- uniqueID: This is a unique sample level ID; duplicates are not allowed in this column
- groupID: This is the groupID which should match to the `contrast_manifest`; duplicates are allowed in this column
- dataType: This is the datatype for the input sample; options are 'gex' 'atac' 'vdj'
- input_dir: This is the input directory for the data files of the sample type (IE "/path/to/sample1/fastq")

An example sampleManifest file is shown below:

| masterID | uniqueID | groupID | dataType | input_dir |
| --- |--- |--- |--- |--- |
| WB_Lysis_1 | sample1 | group1 | gex | /data/CCBR_Pipeliner/Pipelines/TechDev_scRNASeq_Dev2023/test_dir/| WB_Lysis_Granulocytes_3p_Introns_8kCells_fastqs/sample1
| WB_Lysis_1 | sample2 | group1 | gex | /data/CCBR_Pipeliner/Pipelines/TechDev_scRNASeq_Dev2023/test_dir/| WB_Lysis_Granulocytes_3p_Introns_8kCells_fastqs/sample2
| WB_Lysis_2 | sample3 | group2 | gex | /data/CCBR_Pipeliner/Pipelines/TechDev_scRNASeq_Dev2023/test_dir/| WB_Lysis_Granulocytes_3p_Introns_8kCells_fastqs/sample3
| WB_Lysis_2 | sample4 | group2 | gex | /data/CCBR_Pipeliner/Pipelines/TechDev_scRNASeq_Dev2023/test_dir/| WB_Lysis_Granulocytes_3p_Introns_8kCells_fastqs/sample4
| WB_Lysis_3 | sample5 | group3 | gex,/data/CCBR_Pipeliner/Pipelines/TechDev_scRNASeq_Dev2023/test_dir/| WB_Lysis_Granulocytes_3p_Introns_8kCells_fastqs/sample5
| WB_Lysis_1 | sample6 | group1 | atac | /data/CCBR_Pipeliner/Pipelines/TechDev_scRNASeq_Dev2023/test_dir/| WB_Lysis_Granulocytes_3p_Introns_8kCells_fastqs/sample1

### 2.2.2 Contrast Manifest
This manifest will include sample information to performed differential comparisons. A few requirements:

- groups listed must match groups within the `input_manifest` groupID column
- headers should be included for the max number of contrasts. In the example below, the second contrast contains 3 groups, and so the header includes contrast1-contrast3
- multiple groups can be added by increasing the header and adding additional contrasts, as needed

An example contrast file:

| contrast1 | contrast2 | contrast3
| --- | --- |--- |
| group1 | group2
| group1 | group2 | group3