FILE: WB_Lysis_Granulocytes_3p_Introns_8kCells_fastqs.tar
TYPE: GEX
NAME: Whole Blood RBC Lysis for PBMCs and Neutrophils, Granulocytes, 3'
PRODUCT: Single Cell Gene Expression Dataset by Cell Ranger 6.1.0
INFO: Donor Information: healthy female.
Isolation protocol: CG000392 RevA: Isolation of Leukocytes, Bone Marrow and Peripheral Blood Mononuclear Cells for Single Cell RNA Sequencing - Whole Blood Lysis for Granuloctyes track.
Whole transcriptome/Gene Expression libraries were generated as described in the Chromium Next GEM Single Cell 3' Reagent Kits v3.1 (Dual Index) User Guide (CG000204 Rev D).
TUTORIAL: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/tutorials/neutrophils
SOURCE: https://www.10xgenomics.com/resources/datasets/whole-blood-rbc-lysis-for-pbmcs-neutrophils-granulocytes-3-3-1-standard

## WB_1

- WB_Lysis_Granulocytes_3p_Introns_8kCells_S1_L001_I1_001.fastq.gz
- WB_Lysis_Granulocytes_3p_Introns_8kCells_S1_L001_I2_001.fastq.gz
- WB_Lysis_Granulocytes_3p_Introns_8kCells_S1_L001_R1_001.fastq.gz
- WB_Lysis_Granulocytes_3p_Introns_8kCells_S1_L001_R2_001.fastq.gz

## WB_2

- WB_Lysis_Granulocytes_3p_Introns_8kCells_S1_L002_I1_001.fastq.gz
- WB_Lysis_Granulocytes_3p_Introns_8kCells_S1_L002_I2_001.fastq.gz
- WB_Lysis_Granulocytes_3p_Introns_8kCells_S1_L002_R1_001.fastq.gz
- WB_Lysis_Granulocytes_3p_Introns_8kCells_S1_L002_R2_001.fastq.gz

## WB_3 and WB_4 are copies

- WB_1 --> WB_3
- WB_2 --> WB_4

## samples were then subsampled with a set seed

- WB_1 --> sample1,sample2 --> s100
- WB_2 --> sample3,sample4 --> s101
- WB_3 --> sample5,sample6 --> s102
- WB_4 --> sample7,sample8 --> s103

## grouping tests should be

- group1,group2
- group3,group4
- group1+group3,group2+group4
