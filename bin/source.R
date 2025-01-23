




convertHumanGeneList  <- function(gns){
  egs <- mapIds(org.Hs.eg.db, gns, "ENTREZID","SYMBOL")
  mapped <- select(Orthology.eg.db, egs, "Mus.musculus","Homo.sapiens")
  mapped$MUS <- mapIds(org.Mm.eg.db, as.character(mapped$Mus.musculus), "SYMBOL", "ENTREZID")
  return(as.character(unlist(mapped$MUS )))
}



seuratSample = function(so,npcs){


so = NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000,assay = "RNA")
so = ScaleData(so,assay = "RNA")

so = CellCycleScoring(so,s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

so = SCTransform(so)

so <- RunPCA(object = so, features = VariableFeatures(object = so),verbose=F,npcs = 50)

so <- FindNeighbors(so,dims = 1:npcs)
so <- FindClusters(so,dims = 1:npcs, print.output = 0, resolution = 0.8,algorithm = 3)
so <- RunUMAP(so,dims = 1:npcs,n.components = 3)

return(so)
}


runSingleRCell = function(obj,refFile,fineORmain){
obj = DietSeurat(obj, graphs = "umap")
sce = as.SingleCellExperiment(obj,assay = "SCT")
ref = refFile
s = SingleR(test = sce, ref = ref,labels = ref[[fineORmain]])
return(s$pruned.labels)
}


doublets <-function(dfso){

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

runClustIntegrated = function(so,npcs,resolutions,reduction){

so <- RunUMAP(so, reduction = reduction, dims = 1:npcs)
so <- FindNeighbors(so, reduction = reduction, dims = 1:npcs)

for (res in resolutions){
    so <- FindClusters(so,dims = 1:npcs, resolution = res,algorithm = 3)
  }
return(so)

}

singleRClusters = function(obj,refFile,fineORmain){
  avg = AverageExpression(obj,assays = "SCT")
  avg = as.data.frame(avg)
  ref = refFile
  s = SingleR(test = as.matrix(avg),ref = ref,labels = ref[[fineORmain]])

  clustAnnot = s$labels
  names(clustAnnot) = colnames(avg)
  names(clustAnnot) = gsub("SCT.","",names(clustAnnot))

  clustAnnot = clustAnnot[match(obj$seurat_clusters,names(clustAnnot))]
  names(clustAnnot) = colnames(obj)
  obj$clustAnnot = clustAnnot
  return(obj$clustAnnot)
}
