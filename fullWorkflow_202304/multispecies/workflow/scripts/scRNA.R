loc = "/data/CCBR_Pipeliner/db/PipeDB/scrna4.2Rlibs"
library(Seurat,lib.loc=loc)
#library(Routliers,lib.loc=loc)
library(stringr)#,lib.loc=loc)
library("DoubletFinder",lib.loc=loc)
library(SingleR,lib.loc=loc)
library(scRNAseq,lib.loc=loc)
library(SingleCellExperiment,lib.loc=loc)
#library(scater,lib.loc=loc)
library(celldex,lib.loc=loc)
library(Orthology.eg.db,lib.loc=loc)
library(org.Mm.eg.db,lib.loc=loc)
library(org.Hs.eg.db,lib.loc=loc)



args <- commandArgs(trailingOnly = T)


h5 = as.character(args[1])
ref =  as.character(args[2])
outFile = as.character(args[3])
rnaCounts = Read10X_h5(h5)

so <- CreateSeuratObject(counts = rnaCounts)


###Run Seurat Clustering
seuratClustering = function(so){



so$Sample = tail(strsplit(h5,"/")[[1]],1)



so <- RunPCA(object = so, features = VariableFeatures(object = so), do.print = TRUE, pcs.print = 1:5,genes.print = 0,verbose=F,npcs = 30)

npcs = 30


so <- FindNeighbors(so,dims = 1:npcs)
so <- FindClusters(so,dims = 1:npcs, print.output = 0, resolution = 0.1,algorithm = 3)
so <- FindClusters(so,dims = 1:npcs, print.output = 0, resolution = 0.2,algorithm = 3)
so <- FindClusters(so,dims = 1:npcs, print.output = 0, resolution = 0.5,algorithm = 3)
so <- FindClusters(so,dims = 1:npcs, print.output = 0, resolution = 0.6,algorithm = 3)
so <- FindClusters(so,dims = 1:npcs, print.output = 0, resolution = 0.8,algorithm = 3)


so <- RunUMAP(so,dims = 1:npcs,n.components = 3L)


return(so)
}



doublets <-function(dfso){
  dfso <- SCTransform(dfso)
  dfso <- RunPCA(dfso, pc.genes = dfso@var.genes, pcs.print = 0,verbose = F,npcs =10)
  npcs = 10
  dfso <- RunUMAP(dfso, verbose=TRUE,dims = 1:npcs)


  sweep.res.list_kidney <- paramSweep_v3(dfso,PCs = 1:10, sct = T)
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
  bcmvn_kidney <- find.pK(sweep.stats_kidney)
  ## pK Identification (ground-truth) ------------------------------------------------------------------------------------------


  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  homotypic.prop <- modelHomotypic(dfso$annot)
  perc = 0.005 * (length(colnames(dfso))/1000)
  nExp_poi <- round(perc*length(colnames(dfso)))#dfso@cell.names
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  dfso <- doubletFinder_v3(dfso, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE,PCs = 1:10,sct = T)
  pAAN=tail(names(dfso@meta.data),2)[1]
  dfso <- doubletFinder_v3(dfso, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = pAAN,PCs = 1:10,sct = T)

  return(dfso)
}

#convertHumanGeneList <- function(x){

 # require("biomaRt")
  #human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  #mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

 # genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)

 # humanx <- unique(genesV2[, 2])

 # return(humanx)
#}




convertHumanGeneList  <- function(gns){
  egs <- mapIds(org.Hs.eg.db, gns, "ENTREZID","SYMBOL")
  mapped <- select(Orthology.eg.db, egs, "Mus.musculus","Homo.sapiens")
  mapped$MUS <- mapIds(org.Mm.eg.db, as.character(mapped$Mus.musculus), "SYMBOL", "ENTREZID")
  return(as.character(unlist(mapped$MUS )))
}



if(ref=="hg38"){so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")}
if(ref=="mm10"){so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^mt-")}

#so <- subset(so, subset = nFeature_RNA > 200)
#nCount_out = outliers_mad(so$nCount_RNA,threshold = 3)$LL_CI_MAD#
#nFeature_out = outliers_mad(so$nFeature_RNA,threshold = 3)$LL_CI_MAD
#mt_out = outliers_mad(so$percent.mt,threshold = 3)$UL_CI_MAD

#so <- subset(so, subset = nFeature_RNA > nFeature_out & nCount_RNA > nFeature_out & percent.mt < mt_out)


if(ref=="hg38"){
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
}

if(ref=="mm10"){
s.genes <- convertHumanGeneList(cc.genes$s.genes)
g2m.genes <- convertHumanGeneList(cc.genes$g2m.genes)
}


so = SCTransform(so)
so = NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000,assay = "RNA")
so = ScaleData(so,assay = "RNA")

so = CellCycleScoring(so,s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


so = seuratClustering(so)

print("so")
runSingleR = function(obj,refFile,fineORmain){
obj = DietSeurat(obj, graphs = "umap")
sce = as.SingleCellExperiment(obj,assay = "SCT")
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

print("anot")



#dfso = doublets(so)
#so$DF_hi.lo = dfso[[tail(names(dfso@meta.data),1)]]
#so=subset(so,cells=names(so$DF_hi.lo)[so$DF_hi.lo =="Singlet"])

saveRDS(so,outFile)
