import os
import re
import sys
import pandas as pd
import numpy as np

sample = sys.argv[1]
pipeline = sys.argv[2]
rnaRef = sys.argv[3]
vdjRef = sys.argv[4]
atacRef = sys.argv[5]
# output = sys.argv[6]


data = pd.read_excel(
    "/data/khanlab/scrnaPipeliner/SequencingMasterFile_OutsidePatients_singlecell_Copy_revised_final5.xlsx",
    engine="openpyxl",
)

data = data[
    [
        "Project",
        "custom ID",
        "Type of sequencing",
        "Enrichment step",
        "FCID",
        "Library ID",
        "SampleRef",
        "sc cite-seq feature ref",
        "Case Name",
    ]
]

data["finalName"] = "Sample_" + data["Library ID"] + "_" + data["FCID"]
data[
    "fastq"
] = "/data/khanlab/ref/scRNAseq/fastqs/"  # +  "Sample_" + data["Library ID"] + "_" + data["FCID"]

# fastqs = data["fastq"].str.split("/",expand=True).iloc[:,5]
# print(fastqs)

tempData = data.loc[data["Case Name"] == sample]
fcid = tempData.FCID.str.split("_", expand=True)
sc = tempData[tempData["Enrichment step"].str.contains("GEX")]
sc = sc[["fastq", "finalName"]]
sc["type"] = "Gene Expression"
sc.columns = ["fastqs", "sample", "library_type"]
cite = tempData[tempData["Enrichment step"].str.contains("cite")]
cite = cite[["fastq", "finalName"]]
cite["type"] = "Antibody Capture"
cite.columns = ["fastqs", "sample", "library_type"]
scCite = pd.concat([sc, cite])
scCite.to_csv(str(sample + "_libraries.csv"), index=False)


if "rna" in pipeline and "atac" not in pipeline:
    if "cite" not in pipeline:
        toRun = str(
            "cellranger count "
            + "--id="
            + sample
            + "_gex"
            + " "
            + "--libraries="
            + str(sample + "_libraries.csv")
            + " "
            + "--transcriptome="
            + rnaRef
            + " "
            + "--localcores=16 --localmem=64"
        )
        print(toRun)
        os.system(toRun)

    else:
        print(tempData)
        citeRef = tempData[tempData["Enrichment step"].str.contains("cite")]
        citeRef = citeRef["sc cite-seq feature ref"].to_string(index=False)
        citeRef = citeRef.strip(" ")
        print(citeRef)
        toRun = str(
            "cellranger count "
            + "--id="
            + sample
            + "_gex_cite"
            + " "
            + "--libraries="
            + str(sample + "_libraries.csv")
            + " "
            + "--transcriptome="
            + rnaRef
            + " "
            + "--feature-ref="
            + citeRef
            + " "
            + "--localcores=16 --localmem=64"
        )
        print(toRun)
        os.system(toRun)

if "tcr" in pipeline or "bcr" in pipeline:
    vdj = tempData[tempData["Enrichment step"].str.contains("VDJ")]
    vdjFQ = vdj["fastq"].to_string(index=False)
    vdjID = vdj["finalName"].to_string(index=False)

    toRun = str(
        "cellranger vdj "
        + "--id="
        + sample
        + "_tcr"
        + " "
        + "--fastqs="
        + vdjFQ
        + " "
        + "--reference="
        + vdjRef
        + " "
        + "--sample="
        + vdjID
        + " "
        + "--localcores=16 --localmem=64"
    )
    os.system(toRun)


if "atac" in pipeline:
    toRun = str(
        "cellranger-arc count "
        + "--id="
        + sample
        + " "
        + "--libraries="
        + lib
        + " "
        + "--reference="
        + atacRef
        + " "
        + "--localcores=16 --localmem=64"
    )
    print(toRun)
    os.system(toRun)
