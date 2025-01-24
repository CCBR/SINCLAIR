#' Seurat pre-processing
SEURAT_CLUSTERING <- function(so_in, npcs_in) {
  # Runs Principal Component Analysis, FindNeighbors, clustering with the Smart Local Moving algorithm, and UMAP dimensionality reduction
  so <- RunPCA(
    object = so_in,
    features = VariableFeatures(object = so_in),
    verbose = F,
    npcs = 50
  )
  so <- FindNeighbors(so, dims = 1:npcs_in)
  so <- FindClusters(so, resolution = 0.8, algorithm = 3, verbose = TRUE)
  so <- RunUMAP(so, dims = 1:npcs_in, n.components = 3)
  return(so)
}

CONVERT_TO_HUMAN_GENELIST <- function(gns) {
  egs <- mapIds(org.Hs.eg.db, gns, "ENTREZID", "SYMBOL")
  mapped <- select(Orthology.eg.db, egs, "Mus.musculus", "Homo.sapiens")
  mapped$MUS <- mapIds(org.Mm.eg.db, as.character(mapped$Mus.musculus), "SYMBOL", "ENTREZID")
  return(as.character(unlist(mapped$MUS)))
}


RUN_SINGLEr <- function(obj, refFile, fineORmain) {
  obj <- DietSeurat(obj, graphs = "umap")
  sce <- as.SingleCellExperiment(obj, assay = "SCT")
  ref <- refFile
  s <- SingleR(test = sce, ref = ref, labels = ref[[fineORmain]])
  return(s$pruned.labels)
}

#' run batch corrections
RUN_SINGLEr_AVERAGE <- function(obj, refFile, fineORmain) {
  avg <- AverageExpression(obj, assays = "SCT")
  avg <- as.data.frame(avg)
  ref <- refFile
  s <- SingleR(test = as.matrix(avg), ref = ref, labels = ref[[fineORmain]])

  clustAnnot <- s$labels
  names(clustAnnot) <- colnames(avg)
  names(clustAnnot) <- gsub("SCT.", "", names(clustAnnot))

  annotVect <- clustAnnot[match(obj$seurat_clusters, names(clustAnnot))]
  names(annotVect) <- colnames(obj)
  return(annotVect)
}

#' batch correction function used in multiple rmarkdown notebooks
MAIN_BATCH_CORRECTION <- function(so_in, npcs, species, resolution_list, method_in, reduction_in, v_list) {
  # set assay to RNA to avoid double transform/norm
  DefaultAssay(so_in) <- "RNA"

  # integration method for
  ### SCVI
  ### LIGER
  ### harmony,rpca,cca
  if (method_in == "scVIIntegration") {
    print("--running SVII integration")

    so_transform <- NormalizeData(so_in)
    so_variable <- FindVariableFeatures(so_transform)
    so_scaled <- ScaleData(so_variable)
    so_pca <- RunPCA(so_scaled)

    so_integrate <- IntegrateLayers(
      object = so_pca, method = scVIIntegration,
      new.reduction = "integrated.scvi",
      dims = 1:npcs
    )
  } else if (method_in == "LIGER") {
    print("--running LIGER")

    # preprocess
    so_norm <- NormalizeData(so_in)
    so_norm <- FindVariableFeatures(so_norm)
    so_norm <- ScaleData(so_norm, do.center = FALSE)
    so_norm <- RunOptimizeALS(so_norm, k = 20, lambda = 5)
    so_integrate <- RunQuantileNorm(so_norm)
  } else {
    print("--running SCT")

    # use variables to regress, if provided by user
    if (length(v_list) > 0) {
      so_transform <- SCTransform(so_in, vars.to.regress = v_list)
    } else {
      so_transform <- SCTransform(so_in)
    }

    # runPCA
    so_pca <- RunPCA(so_transform)

    so_integrate <- IntegrateLayers(
      object = so_pca, method = get(method_in),
      normalization.method = "SCT",
      verbose = F, new.reduction = reduction_in
    )
  }

  # run neighbors, clusters
  so <- FindNeighbors(so_integrate, reduction = reduction_in, dims = 1:npcs)
  for (res in resolution_list) {
    so <- FindClusters(so, resolution = res, algorithm = 3)
  }

  # reduction
  so <- RunUMAP(so, reduction = reduction_in, dims = 1:npcs)

  # relabel
  if (species == "hg38" || species == "hg19") {
    so$clustAnnot_HPCA_main <- RUN_SINGLEr_AVERAGE(so, celldex::HumanPrimaryCellAtlasData(), "label.main")
    so$clustAnnot_HPCA <- RUN_SINGLEr_AVERAGE(so, celldex::HumanPrimaryCellAtlasData(), "label.fine")
    so$clustAnnot_BP_encode_main <- RUN_SINGLEr_AVERAGE(so, celldex::BlueprintEncodeData(), "label.main")
    so$clustAnnot_BP_encode <- RUN_SINGLEr_AVERAGE(so, celldex::BlueprintEncodeData(), "label.fine")
    so$clustAnnot_monaco_main <- RUN_SINGLEr_AVERAGE(so, celldex::MonacoImmuneData(), "label.main")
    so$clustAnnot_monaco <- RUN_SINGLEr_AVERAGE(so, celldex::MonacoImmuneData(), "label.fine")
    so$clustAnnot_immu_cell_exp_main <- RUN_SINGLEr_AVERAGE(so, celldex::DatabaseImmuneCellExpressionData(), "label.main")
    so$clustAnnot_immu_cell_exp <- RUN_SINGLEr_AVERAGE(so, celldex::DatabaseImmuneCellExpressionData(), "label.fine")
  } else if (species == "mm10") {
    so$clustAnnot_immgen_main <- RUN_SINGLEr_AVERAGE(so, celldex::ImmGenData(), "label.main")
    so$clustAnnot_immgen <- RUN_SINGLEr_AVERAGE(so, celldex::ImmGenData(), "label.fine")
    so$clustAnnot_mouseRNAseq_main <- RUN_SINGLEr_AVERAGE(so, celldex::MouseRNAseqData(), "label.main")
    so$clustAnnot_mouseRNAseq <- RUN_SINGLEr_AVERAGE(so, celldex::MouseRNAseqData(), "label.fine")
  }
  return(so)
}

#' Integration Report Functions
OBJECT_SELECT <- function(id) {
  obj <- switch(id,
    "merged" = so_merged,
    "integrated" = so_integrated,
    "rpca" = so_rpca,
    "harmony" = so_harmony,
    "scvi" = so_scvi,
    "liger" = so_liger
  )
  return(obj)
}
NAME_SELECT <- function(id) {
  obj <- switch(id,
    "merged" = "Before Batch Correction",
    "integrated" = "Integrated CCA",
    "rpca" = "RPCA",
    "harmony" = "Harmony (Sample)",
    "scvi" = "single-cell Variational Inference",
    "liger" = "Linked Inference of Genomic Experimental Relationships"
  )
  return(obj)
}
