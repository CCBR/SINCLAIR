import os
import re
import sys
import pandas as pd
import numpy as np

sample = sys.argv[1]
workpath = sys.argv[2]
pipeline = sys.argv[3]
ref = sys.argv[4]
output = sys.argv[5]


# data = pd.read_csv(str(workpath+"groups.tab"), sep='\t', low_memory=False, encoding= 'unicode_escape',header = None)

# data.columns = ["sample","group","name","type"]


if "rna" in pipeline and "atac" not in pipeline:
    if "cite" not in pipeline:
        toRun = str(
            "Rscript "
            + "workflow/scripts/scRNA.R "
            + "cellrangerOut/"
            + sample
            + "/outs/filtered_feature_bc_matrix.h5 "
            + ref
            + " "
            + output
        )
        # 		toRun = str("Rscript "+"--id="+sample+ "_gex"  +" "  + "--libraries=" + str(sample + "_libraries.csv") + " " + "--transcriptome=" + rnaRef + " "  + "--localcores=16 --localmem=64")
        print(toRun)
        os.system(toRun)

    else:
        print(tempData)
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


# if "rna" in pipeline and "atac" not in pipeline:
# 	if "cite" not in pipeline:
# 		toRun = str("cellranger count "+"--id="+sample+ "_gex"  +" "  + "--libraries=" + str(sample + "_libraries.csv") + " " + "--transcriptome=" + rnaRef + " "  + "--localcores=16 --localmem=64")
# 		print(toRun)
# 		os.system(toRun)

# 	else:
# 		print(tempData)
# 		toRun = str("cellranger count "+"--id="+sample+  "_gex_cite"  +" " + "--libraries=" + str(sample + "_libraries.csv") + " " + "--transcriptome=" + rnaRef + " " + "--feature-ref=" + citeRef + " "  + "--localcores=16 --localmem=64")
# 		print(toRun)
# 		os.system(toRun)
