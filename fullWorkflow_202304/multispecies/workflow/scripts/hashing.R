sc4Lib = "/data/CCBR_Pipeliner/db/PipeDB/Rlibrary_4.1.0_scRNA"
library('rlang', lib.loc=sc4Lib)
library(Seurat, lib.loc=sc4Lib)
library(SingleR, lib.loc=sc4Lib)
library(scRNAseq)
library(SingleCellExperiment, lib.loc=sc4Lib)

args <- commandArgs(trailingOnly = T)


h5 = as.character(args[1])
ref =  as.character(args[2])
outFile = as.character(args[3])
rnaCounts = Read10X_h5(h5)


so = CreateSeuratObject(h5$'Gene Expression')

so = NormalizeData(so)

so = FindVariableFeatures(so, selection.method = "mean.var.plot")

so = ScaleData(so, features = VariableFeatures(so))

so[["HTO"]] <- CreateAssayObject(counts = h5$'Antibody Capture'[grep("HTO",rownames(h5$'Antibody Capture')),] )

so = NormalizeData(so, assay = "HTO", normalization.method = "CLR")

so <- HTODemux(so, assay = "HTO", positive.quantile = 0.99)

Idents(so) = "HTO_classification.global"

so = subset(so, idents = "Singlet")


runSingleR = function(obj,refFile,fineORmain){
sce = as.SingleCellExperiment(obj,assay = "RNA")
ref = refFile
s = SingleR(test = sce, ref = ref,labels = ref[[fineORmain]])
return(s$pruned.labels)
print(head(s$pruned.labels))
}


if(ref == "hg38"){
so$HPCA_main <- runSingleR(so,celldex::HumanPrimaryCellAtlasData(),"label.main")
so$HPCA <-  runSingleR(so,celldex::HumanPrimaryCellAtlasData(),"label.fine")
so$BP_encode_main <-  runSingleR(so,celldex::BlueprintEncodeData(),"label.main")
so$BP_encode <-  runSingleR(so,celldex::BlueprintEncodeData(),"label.fine")
so$monaco_main <-  runSingleR(so,celldex::MonacoImmuneData(),"label.main")
so$monaco <-     runSingleR(so,celldex::MonacoImmuneData(),"label.fine")
so$immu_cell_exp_main <-  runSingleR(so,celldex::DatabaseImmuneCellExpressionData(),"label.main")
so$immu_cell_exp <- runSingleR(so,celldex::DatabaseImmuneCellExpressionData(),"label.fine")
so$annot = so$HPCA_main
}

if(ref == "mm10"){

so$immgen_main <-  runSingleR(so,celldex::ImmGenData(),"label.main")
so$immgen <- runSingleR(so,celldex::ImmGenData(),"label.fine")
so$mouseRNAseq_main <-  runSingleR(so,celldex::MouseRNAseqData(),"label.main")
so$mouseRNAseq <- runSingleR(so,celldex::MouseRNAseqData(),"label.fine")

so$annot = so$immgen_main

}



DF = as.data.frame(as.matrix(so@assays$HTO@data))
so$assignHTO = rownames(DF)[apply(DF,2,which.max)]
so.list = SplitObject(so,split.by="assignHTO")

#remove samples with less than 100 cells
ncells = vector()
for(n in 1:length(so.list) ){ncells = append(ncells,ncol(so.list[[n]])) }

so.list = so.list[ncells > 100]

so.list <- lapply(X = so.list, FUN = NormalizeData,normalization.method = "LogNormalize", scale.factor = 10000,assay = "RNA")
so.list <- lapply(X = so.list, FUN = SCTransform,variable.features.rv.th = 1.3,variable.features.n = NULL)
so.list <- lapply(X = so.list, FUN = RunPCA,npcs = 30)
so.list <- lapply(X = so.list, FUN = RunUMAP, dims = 1:30)
so.list <- lapply(X = so.list, FUN = FindNeighbors, dims = 1:30)
so.list <- lapply(X = so.list, FUN = FindClusters, algorithm = 3)
so.list <- lapply(X = so.list, FUN = saveRDS, paste0("x"))

for(n in 1:length(so.list)){}

for(n in 1:length(so.list)){saveRDS(so.list[[n]] , paste0("x_",names(so.list)[n]),".rds")}
