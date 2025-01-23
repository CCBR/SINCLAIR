.libPaths( c("/data/CCBR_Pipeliner/db/PipeDB/scrna5",.libPaths()))


library(Seurat)
library(stringr)
library(SingleR)
library(scRNAseq)
library(SingleCellExperiment)
library(celldex)
library(Orthology.eg.db)

options(future.globals.maxSize = 1e12)
#options(Seurat.object.assay.version = "v5")


args <- commandArgs(trailingOnly = TRUE)

sampleDir <- as.character(args[1])
outDirMerge = as.character(args[2])
ref = as.character(args[3])
contrasts = as.character(args[4])




file.names <- dir(path = sampleDir,pattern ="rds")
groupFile = read.delim("groups.tab",header=F,stringsAsFactors = F)
groupFile=groupFile[groupFile$V2 %in% stringr::str_split_fixed(contrasts,pattern = "-",n = Inf)[1,],]

splitFiles = gsub(".rds","",file.names)#str_split_fixed(file.names,pattern = "[.rd]",n = 2)
#splitFiles = stringr::str_split_fixed(splitFiles,pattern = "__",n = Inf)[,1]
file.names=file.names[match(groupFile$V1,splitFiles,nomatch = F)]
print(groupFile$V1)
print(splitFiles)
print(file.names)


readObj = list()
for (obj in file.names) {
  Name=strsplit(obj,".rds")[[1]][1]
  assign(paste0("S_",Name),readRDS(paste0(sampleDir,"/",obj)))
  readObj = append(readObj,paste0("S_",Name))
 }



combinedObj.list=list()
i=1
for (p in readObj){
  combinedObj.list[[p]] <- eval(parse(text = readObj[[i]]))
  combinedObj.list[[p]]$Sample = names(combinedObj.list)[i]
  i <- i + 1
 }

reference.list <- combinedObj.list[unlist(readObj)]

print(length(reference.list))
print(reference.list)

selectFeatures <- SelectIntegrationFeatures(object.list = reference.list, nfeatures = 3000)

combinedObj.integratedRNA = merge(reference.list[[1]],reference.list[2:length(reference.list)])
VariableFeatures(combinedObj.integratedRNA) = selectFeatures

saveRDS(combinedObj.integratedRNA,outDirMerge)
