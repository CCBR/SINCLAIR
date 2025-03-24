loc = "/data/CCBR_Pipeliner/db/PipeDB/scrna4.2Rlibs"
library(htmltools,lib.loc=loc)
library(rlang,lib.loc=loc)
library(globals,lib.loc=loc)
library(stringr,lib.loc=loc)
library(harmony,lib.loc=loc)
library("SeuratWrappers",lib.loc=loc)
library(Seurat,lib.loc=loc)
library("SingleR",lib.loc=loc)
library(scRNAseq,lib.loc=loc)
library(SingleCellExperiment,lib.loc=loc)
library(dplyr)
library(Matrix)
library(tools)


args <- commandArgs(trailingOnly = TRUE)

matrix <- as.character(args[1])
outDirMerge = as.character(args[2])
outDirSample = as.character(args[3])
outDirGroup = as.character(args[4])
ref = as.character(args[5])
contrasts = as.character(args[6])



options(future.globals.maxSize = 64000 * 1024^2)


file.names <- dir(path = matrix,pattern ="rds")

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
  assign(paste0("S_",Name),readRDS(paste0(matrix,"/",obj)))
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
combinedObj.integratedRNA$Sample = stringr::str_split_fixed(combinedObj.integratedRNA$Sample,pattern = "__",n = Inf)[,1]
combinedObj.integratedRNA$groups = groupFile$V2[match(combinedObj.integratedRNA$Sample,  paste0("S_",groupFile$V1),nomatch = F)]

combinedObj.integratedRNA <- RunPCA(object = combinedObj.integratedRNA, npcs = 50, verbose = FALSE)
harmSample <- RunHarmony(combinedObj.integratedRNA, group.by.vars = "Sample",assay.use = "SCT")
harmGroup <- RunHarmony(combinedObj.integratedRNA, group.by.vars = "groups",assay.use = "SCT")





runInt = function(obj,npcs,harm){
obj <- RunPCA(object = obj, npcs = npcs, verbose = FALSE)

if(harm=="Yes") {obj <- FindNeighbors(obj,dims = 1:npcs,reduction = "harmony")}
if(harm=="No") {obj <- FindNeighbors(obj,dims = 1:npcs)}

obj <- FindClusters(obj,dims = 1:npcs, print.output = 0, resolution = 0.1,algorithm = 3)
obj <- FindClusters(obj,dims = 1:npcs, print.output = 0, resolution = 0.2,algorithm = 3)
obj <- FindClusters(obj,dims = 1:npcs, print.output = 0, resolution = 0.5,algorithm = 3)
obj <- FindClusters(obj,dims = 1:npcs, print.output = 0, resolution = 0.6,algorithm = 3)
obj <- FindClusters(obj,dims = 1:npcs, print.output = 0, resolution = 0.8,algorithm = 3)
colnames(obj@meta.data) = gsub("integrated_snn_res","SLM_int_snn_res",colnames(obj@meta.data))
colnames(obj@meta.data) = gsub("SCT_snn_res","SLM_SCT_snn_res",colnames(obj@meta.data))


obj <- FindClusters(obj,dims = 1:npcs, print.output = 0, resolution = 0.1,algorithm = 4)#
obj <- FindClusters(obj,dims = 1:npcs, print.output = 0, resolution = 0.2,algorithm = 4)
obj <- FindClusters(obj,dims = 1:npcs, print.output = 0, resolution = 0.5,algorithm = 4)
obj <- FindClusters(obj,dims = 1:npcs, print.output = 0, resolution = 0.6,algorithm = 4)
obj <- FindClusters(obj,dims = 1:npcs, print.output = 0, resolution = 0.8,algorithm = 4)
colnames(obj@meta.data) = gsub("integrated_snn_res","Leiden_int_snn_res",colnames(obj@meta.data))
colnames(obj@meta.data) = gsub("SCT_snn_res","Leiden_SCT_snn_res",colnames(obj@meta.data))


if(harm=="Yes") {obj <- RunUMAP(object = obj,dims = 1:npcs,n.components = 3,reduction = "harmony")}
if(harm=="No") {obj <- RunUMAP(object = obj, reduction = "pca",dims = 1:npcs,n.components = 3)}


#obj$groups = groupFile$V2[match(obj$Sample,  groupFile$V1,nomatch = F)]

runSingleR = function(obj,refFile,fineORmain){
avg = AverageExpression(obj,assays = "SCT")
avg = as.data.frame(avg)
ref = refFile
s = SingleR(test = as.matrix(avg),ref = ref,labels = ref[[fineORmain]])

clustAnnot = s$labels
names(clustAnnot) = colnames(avg)
names(clustAnnot) = gsub("SCT.","",names(clustAnnot))

obj$clustAnnot = clustAnnot[match(obj$seurat_clusters,names(clustAnnot))]
return(obj$clustAnnot)
}

if(ref == "hg38"){
obj$clustAnnot_HPCA_main <- runSingleR(obj,celldex::HumanPrimaryCellAtlasData(),"label.main")
obj$clustAnnot_HPCA <-  runSingleR(obj,celldex::HumanPrimaryCellAtlasData(),"label.fine")
obj$clustAnnot_BP_encode_main <-  runSingleR(obj,celldex::BlueprintEncodeData(),"label.main")
obj$clustAnnot_BP_encode <-  runSingleR(obj,celldex::BlueprintEncodeData(),"label.fine")
obj$clustAnnot_monaco_main <-  runSingleR(obj,celldex::MonacoImmuneData(),"label.main")
obj$clustAnnot_monaco <- runSingleR(obj,celldex::MonacoImmuneData(),"label.fine")
obj$clustAnnot_immu_cell_exp_main <-  runSingleR(obj,celldex::DatabaseImmuneCellExpressionData(),"label.main")
obj$clustAnnot_immu_cell_exp <- runSingleR(obj,celldex::DatabaseImmuneCellExpressionData(),"label.fine")
}

if(ref == "mm10"){

obj$clustAnnot_immgen_main <-  runSingleR(obj,celldex::ImmGenData(),"label.main")
obj$clustAnnot_immgen <- runSingleR(obj,celldex::ImmGenData(),"label.fine")
obj$clustAnnot_mouseRNAseq_main <-  runSingleR(obj,celldex::MouseRNAseqData(),"label.main")
obj$clustAnnot_mouseRNAseq <- runSingleR(obj,celldex::MouseRNAseqData(),"label.fine")

}
return(obj)
}


combinedObj.integratedRNA = runInt(combinedObj.integratedRNA,50,"No")
saveRDS(combinedObj.integratedRNA,outDirMerge)

harmGroup = runInt(harmGroup,50,"Yes")
harmSample = runInt(harmSample,50,"Yes")


saveRDS(combinedObj.integratedRNA,outDirMerge)
saveRDS(harmGroup,outDirGroup)
saveRDS(harmSample,outDirSample)
