loc = "/data/CCBR_Pipeliner/db/PipeDB/scrna4.2Rlibs"
library(htmltools,lib.loc=loc)
library(Seurat,lib.loc=loc)
library(stringr)#,lib.loc=loc)
library("DoubletFinder",lib.loc=loc)
library(SingleR,lib.loc=loc)
library(scRNAseq,lib.loc=loc)
library(SingleCellExperiment,lib.loc=loc)
library(celldex,lib.loc=loc)
library(Orthology.eg.db,lib.loc=loc)
library(org.Mm.eg.db,lib.loc=loc)
library(org.Hs.eg.db,lib.loc=loc)
library(flexmix,lib.loc=loc)
library(SeuratWrappers,lib.loc=loc)
library(djvdj,lib.loc=loc)

args <- commandArgs(trailingOnly = T)


h5 = as.character(args[1])
ref =  as.character(args[2])
outFile = as.character(args[3])
rnaCounts = Read10X_h5(h5)

if(class(rnaCounts) == "list") {rnaCounts = rnaCounts$'Gene Expression'}

so <- CreateSeuratObject(rnaCounts)


groupFile = read.delim("groups.tab",header=F,stringsAsFactors = F)

sample = groupFile$V3[groupFile$V1 == tail(strsplit(h5,"/")[[1]],3)[1] & groupFile$V4 == "gex"]

print(groupFile[groupFile$V3 == sample & groupFile$V4 == "vdj",])
#VDJ INPUT STARTS HERE
if (nrow(groupFile[groupFile$V3 == sample & groupFile$V4 == "vdj",]) > 0 ) {

	tcrSamples = groupFile$V1[groupFile$V3 == sample & groupFile$V4 == "vdj"]
	tcrSamples = paste0("cellrangerOut/",tcrSamples,"/outs")

	so = import_vdj(input = so, vdj_dir = tcrSamples,  filter_paired = FALSE  )

}
