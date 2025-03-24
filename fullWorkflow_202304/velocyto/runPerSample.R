loc = "/data/CCBR_Pipeliner/db/PipeDB/scrna4.2Rlibs"

library(Seurat,lib.loc = loc)
library(SeuratDisk,lib.loc = loc)
library(SeuratWrappers,lib.loc = loc)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
rdsFile = as.character(args[1])
so = readRDS(rdsFile)
Idents(so) = "Sample"

file.names = unique(so$Sample)
file.names = gsub("S_","",file.names)
for (file in file.names){
so_sub = subset(so,idents = paste0("S_",file))
bcs = str_split_fixed(colnames(so_sub),"-",2)[,1]
bcs = paste0(file,":",bcs,"x")
so_sub = RenameCells(so_sub, new.names = bcs)

  if (exists("dataset")){
	ldat <- ReadVelocity(file = paste0("cellrangerOut/",file,"/velocyto/",file,".loom"))
	bm <- as.Seurat(x = ldat)
	bm[["RNA"]] <- bm[["spliced"]]
	bm <- SCTransform(bm)
	bm = subset(bm,cells = bcs)
	bm = AddMetaData(bm, so_sub@meta.data[,9:length(colnames(so_sub@meta.data))])
	umap = rbind(umap, so_sub@reductions$umap@cell.embeddings)
	DefaultAssay(bm) <- "RNA"
    dataset<-merge(dataset, bm)
  }

  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
ldat <- ReadVelocity(file = paste0("cellrangerOut/",file,"/velocyto/",file,".loom"))
bm <- as.Seurat(x = ldat)
bm[["RNA"]] <- bm[["spliced"]]
bm <- SCTransform(bm)
bm = subset(bm,cells = bcs)
bm = AddMetaData(bm, so_sub@meta.data[,9:length(colnames(so_sub@meta.data))])
umap = so_sub@reductions$umap@cell.embeddings

DefaultAssay(bm) <- "RNA"
dataset = bm
  }

}

rownames(umap) = colnames(dataset)
umap = umap[,1:2]
dataset[['UMAP']] = CreateDimReducObject(embeddings = umap, key = "umap", global = T, assay = "RNA")
#dataset@reductions$umap@cell.embeddings = umap

#file.names <- dir(path = "/data/khanlab/projects/abdalla/reRunKathrineSCRNA/velo",pattern =".rds")

rdsFile = tail(str_split(rdsFile,pattern = "/")[[1]],1)
rdsFile = gsub(".rds","",rdsFile)
SaveH5Seurat(dataset, filename = paste0("velocyto/",rdsFile,".h5Seurat"))
Convert(paste0("velocyto/",rdsFile,".h5Seurat"), dest = "h5ad")
#saveRDS(dataset,"x.rds")
