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

MAIN_PROCESS_SO <- function(so_in, species, npcs_in) {
  # assign genes depending on species input
  if (species == "hg38" || species == "hg19") {
    print("--proccesing human data")
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
  } else if (species == "mm10") {
    print("--proccesing mouse data")
    s.genes <- CONVERT_TO_HUMAN_GENELIST(cc.genes$s.genes)
    g2m.genes <- CONVERT_TO_HUMAN_GENELIST(cc.genes$g2m.genes)
  }

  # process
  so_1 <- NormalizeData(so_in,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    assay = "RNA"
  )
  so_2 <- ScaleData(so_1, assay = "RNA")
  so_3 <- CellCycleScoring(so_2,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = TRUE
  )
  so_4 <- SCTransform(so_3)
  so_out <- SEURAT_CLUSTERING(so_4, npcs_in)
  return(so_out)
}

RUN_SINGLEr <- function(obj, refFile, fineORmain) {
  obj <- DietSeurat(obj, graphs = "umap")
  sce <- as.SingleCellExperiment(obj, assay = "SCT")
  ref <- refFile
  s <- SingleR(test = sce, ref = ref, labels = ref[[fineORmain]])
  return(s$pruned.labels)
}

MAIN_SINGLER <- function(so_in, species) {
  if (species == "hg38" || species == "hg19") {
    so_in$HPCA_main <- RUN_SINGLEr(so_in, celldex::HumanPrimaryCellAtlasData(), "label.main")
    so_in$HPCA <- RUN_SINGLEr(so_in, celldex::HumanPrimaryCellAtlasData(), "label.fine")
    so_in$BP_encode_main <- RUN_SINGLEr(so_in, celldex::BlueprintEncodeData(), "label.main")
    so_in$BP_encode <- RUN_SINGLEr(so_in, celldex::BlueprintEncodeData(), "label.fine")
    so_in$monaco_main <- RUN_SINGLEr(so_in, celldex::MonacoImmuneData(), "label.main")
    so_in$monaco <- RUN_SINGLEr(so_in, celldex::MonacoImmuneData(), "label.fine")
    so_in$immu_cell_exp_main <- RUN_SINGLEr(
      so_in, celldex::DatabaseImmuneCellExpressionData(),
      "label.main"
    )
    so_in$immu_cell_exp <- RUN_SINGLEr(
      so_in, celldex::DatabaseImmuneCellExpressionData(),
      "label.fine"
    )
    so_in$annot <- so_in$HPCA_main
  } else if (species == "mm10") {
    so_in$immgen_main <- RUN_SINGLEr(so_in, celldex::ImmGenData(), "label.main")
    so_in$immgen <- RUN_SINGLEr(so_in, celldex::ImmGenData(), "label.fine")
    so_in$mouseRNAseq_main <- RUN_SINGLEr(so_in, celldex::MouseRNAseqData(), "label.main")
    so_in$mouseRNAseq <- RUN_SINGLEr(so_in, celldex::MouseRNAseqData(), "label.fine")
    so_in$annot <- so_in$immgen_main
  }
  return(so_in)
}

MAIN_DOUBLETS <- function(so_in, run_doublet_finder) {
  if (run_doublet_finder == "Y") {
    sweep.res.list_kidney <- paramSweep_v3(so_in, PCs = 1:10, sct = T)
    sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
    bcmvn_kidney <- find.pK(sweep.stats_kidney)

    ## Homotypic Doublet Proportion Estimate
    homotypic.prop <- modelHomotypic(so_in$annot)
    perc <- 0.005 * (length(colnames(so_in)) / 1000)
    nExp_poi <- round(perc * length(colnames(so_in))) # dfso@cell.names
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

    ## Run DoubletFinder with varying classification stringencies
    dfso <- doubletFinder_v3(so_in,
      pN = 0.25, pK = 0.09,
      nExp = nExp_poi,
      reuse.pANN = FALSE, PCs = 1:10, sct = T
    )

    pAAN <- tail(names(dfso@meta.data), 2)[1]
    dfso <- doubletFinder_v3(dfso,
      pN = 0.25, pK = 0.09,
      nExp = nExp_poi.adj,
      reuse.pANN = pAAN, PCs = 1:10, sct = T
    )
    so_in$DF_hi.lo <- dfso[[tail(names(dfso@meta.data), 1)]]
    so_in <- subset(so_in, cells = names(so_in$DF_hi.lo)[so_in$DF_hi.lo == "Singlet"])
  }

  return(so_in)
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
