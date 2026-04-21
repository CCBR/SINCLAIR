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
  mapped$MUS <- mapIds(
    org.Mm.eg.db,
    as.character(mapped$Mus.musculus),
    "SYMBOL",
    "ENTREZID"
  )
  return(as.character(unlist(mapped$MUS)))
}


RUN_SINGLEr <- function(obj, refFile, fineORmain) {
  obj <- DietSeurat(obj, graphs = "umap")
  sce <- as.SingleCellExperiment(obj, assay = "SCT")
  ref <- refFile
  s <- SingleR(test = sce, ref = ref, labels = ref[[fineORmain]])
  return(s$pruned.labels)
}

fetch_celldex_ref <- function(ref_name) {
  ref <- switch(
    ref_name,
    "hpca" = ,
    "HumanPrimaryCellAtlasData" = celldex::fetchReference(
      "hpca",
      version = "2024-02-26",
      realize.assays = TRUE,
      cache = "./"
    ),
    "blueprint_encode" = ,
    "BP_encode" = ,
    "bpencode" = ,
    "BlueprintEncodeData" = celldex::fetchReference(
      "blueprint_encode",
      "2024-02-26",
      realize.assays = TRUE,
      cache = "./"
    ),
    "monaco" = ,
    "MonacoImmuneData" = celldex::fetchReference(
      "monaco_immune",
      "2024-02-26",
      realize.assays = TRUE,
      cache = "./"
    ),
    "immu_cell_exp" = ,
    "DatabaseImmuneCellExpressionData" = ,
    "dice" = celldex::fetchReference(
      "dice",
      "2024-02-26",
      realize.assays = TRUE,
      cache = "./"
    ),
    "immgen" = ,
    "ImmGenData" = celldex::fetchReference(
      "immgen",
      "2024-02-26",
      realize.assays = TRUE,
      cache = "./"
    ),
    "mouseRNAseq" = ,
    "MouseRNAseqData" = celldex::fetchReference(
      "mouse_rnaseq",
      "2024-02-26",
      realize.assays = TRUE,
      cache = "./"
    )
  )
  return(ref)
}

MAIN_SINGLER <- function(so_in, species, cache_path = NULL) {
  if (species == "hg38" || species == "hg19") {
    so_in$HPCA_main <- RUN_SINGLEr(
      so_in,
      fetch_celldex_ref("hpca"),
      "label.main"
    )
    so_in$HPCA <- RUN_SINGLEr(so_in, fetch_celldex_ref("hpca"), "label.fine")
    so_in$BP_encode_main <- RUN_SINGLEr(
      so_in,
      fetch_celldex_ref("BP_encode"),
      "label.main"
    )
    so_in$BP_encode <- RUN_SINGLEr(
      so_in,
      fetch_celldex_ref("BP_encode"),
      "label.fine"
    )
    so_in$monaco_main <- RUN_SINGLEr(
      so_in,
      fetch_celldex_ref("monaco"),
      "label.main"
    )
    so_in$monaco <- RUN_SINGLEr(
      so_in,
      fetch_celldex_ref("monaco"),
      "label.fine"
    )
    so_in$immu_cell_exp_main <- RUN_SINGLEr(
      so_in,
      fetch_celldex_ref("dice"),
      "label.main"
    )
    so_in$immu_cell_exp <- RUN_SINGLEr(
      so_in,
      fetch_celldex_ref("dice"),

      "label.fine"
    )
    so_in$annot <- so_in$HPCA_main
  } else if (species == "mm10") {
    so_in$immgen_main <- RUN_SINGLEr(
      so_in,
      fetch_celldex_ref("immgen"),
      "label.main"
    )
    so_in$immgen <- RUN_SINGLEr(
      so_in,
      fetch_celldex_ref("immgen"),
      "label.fine"
    )
    so_in$mouseRNAseq_main <- RUN_SINGLEr(
      so_in,
      fetch_celldex_ref("mouseRNAseq"),
      "label.main"
    )
    so_in$mouseRNAseq <- RUN_SINGLEr(
      so_in,
      fetch_celldex_ref("mouseRNAseq"),
      "label.fine"
    )
    so_in$annot <- so_in$immgen_main
  }
  return(so_in)
}


MAIN_DOUBLETS <- function(so_in, run_doublet_finder) {
  if (run_doublet_finder == "Y") {
    sweep.res.list_kidney <- paramSweep(so_in, PCs = 1:10, sct = TRUE)
    sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
    bcmvn_kidney <- find.pK(sweep.stats_kidney)

    ## Homotypic Doublet Proportion Estimate
    homotypic.prop <- modelHomotypic(so_in$annot)
    perc <- 0.005 * (length(colnames(so_in)) / 1000)
    nExp_poi <- round(perc * length(colnames(so_in))) # dfso@cell.names
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

    ## Run DoubletFinder with varying classification stringencies
    dfso <- doubletFinder(
      so_in,
      pN = 0.25,
      pK = 0.09,
      nExp = nExp_poi,
      reuse.pANN = NULL, # https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/244
      PCs = 1:10,
      sct = TRUE
    )

    pAAN <- tail(names(dfso@meta.data), 2)[1]
    dfso <- doubletFinder(
      dfso,
      pN = 0.25,
      pK = 0.09,
      nExp = nExp_poi.adj,
      reuse.pANN = pAAN,
      PCs = 1:10,
      sct = T
    )
    so_in$DF_hi.lo <- dfso[[tail(names(dfso@meta.data), 1)]]
    so_in <- subset(
      so_in,
      cells = names(so_in$DF_hi.lo)[so_in$DF_hi.lo == "Singlet"]
    )
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
  names(clustAnnot) <- gsub("^g", "", names(clustAnnot))

  annotVect <- clustAnnot[match(obj$seurat_clusters, names(clustAnnot))]
  names(annotVect) <- colnames(obj)
  return(annotVect)
}

#' batch correction function used in multiple rmarkdown notebooks

MAIN_BATCH_CORRECTION <- function(
  so_in,
  npcs,
  species,
  resolution_list,
  method_in,
  reduction_in,
  v_list = NULL,
  cache_path = NULL
) {
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
      object = so_pca,
      method = scVIIntegration,
      new.reduction = "integrated.scvi",
      dims = 1:npcs
    )
  } else if (method_in == "LIGER") {
    print("--running LIGER")
    #New catch for rliger version, with updated code from
    if (packageVersion("rliger") < "2.0") {
      # preprocess
      so_norm <- Seurat::NormalizeData(so_in)
      so_norm <- Seurat::FindVariableFeatures(so_norm)
      so_norm <- Seurat::ScaleData(so_norm, do.center = FALSE)
      so_norm <- Seurat::RunOptimizeALS(so_norm, k = npcs, lambda = 5)
      so_integrate <- Seurat::RunQuantileNorm(so_norm)
    } else {
      so_norm <- rliger::normalize(so_in)
      so_norm <- rliger::selectGenes(so_norm)
      so_norm <- rliger::scaleNotCenter(so_norm)
      so_norm <- rliger::runINMF(so_norm, k = npcs)
      so_integrate <- rliger::quantileNorm(so_norm)
    }
  } else {
    print("--running SCT")

    # vars.to.regress is NULL by default
    so_transform <- SCTransform(so_in, vars.to.regress = v_list)

    # runPCA
    so_pca <- RunPCA(so_transform)

    so_integrate <- IntegrateLayers(
      object = so_pca,
      method = get(method_in),
      normalization.method = "SCT",
      verbose = F,
      new.reduction = reduction_in
    )
  }

  # run neighbors, clusters
  so <- FindNeighbors(so_integrate, reduction = reduction_in, dims = 1:npcs)
  for (res in resolution_list) {
    so <- FindClusters(so, resolution = res, algorithm = 3)
  }

  # reduction
  so <- RunUMAP(so, reduction = reduction_in, dims = 1:npcs)

  # relabel with cluster-level annotations (uses averaged expression within each cluster)
  if (dir.exists(cache_path)) {
    gypsum::cacheDirectory(cache_path)
  }

  if (species == "hg38" || species == "hg19") {
    so$clustAnnot_HPCA_main <- RUN_SINGLEr_AVERAGE(
      so,
      fetch_celldex_ref("hpca"),
      "label.main"
    )
    so$clustAnnot_HPCA <- RUN_SINGLEr_AVERAGE(
      so,
      fetch_celldex_ref("hpca"),
      "label.fine"
    )
    so$clustAnnot_BP_encode_main <- RUN_SINGLEr_AVERAGE(
      so,
      fetch_celldex_ref("BP_encode"),
      "label.main"
    )
    so$clustAnnot_BP_encode <- RUN_SINGLEr_AVERAGE(
      so,
      fetch_celldex_ref("BP_encode"),
      "label.fine"
    )
    so$clustAnnot_monaco_main <- RUN_SINGLEr_AVERAGE(
      so,
      fetch_celldex_ref("monaco"),
      "label.main"
    )
    so$clustAnnot_monaco <- RUN_SINGLEr_AVERAGE(
      so,
      fetch_celldex_ref("monaco"),
      "label.fine"
    )
    so$clustAnnot_immu_cell_exp_main <- RUN_SINGLEr_AVERAGE(
      so,
      fetch_celldex_ref("dice"),
      "label.main"
    )
    so$clustAnnot_immu_cell_exp <- RUN_SINGLEr_AVERAGE(
      so,
      fetch_celldex_ref("dice"),
      "label.fine"
    )
  } else if (species == "mm10") {
    so$clustAnnot_immgen_main <- RUN_SINGLEr_AVERAGE(
      so,
      fetch_celldex_ref("immgen"),
      "label.main"
    )
    so$clustAnnot_immgen <- RUN_SINGLEr_AVERAGE(
      so,
      fetch_celldex_ref("immgen"),
      "label.fine"
    )
    so$clustAnnot_mouseRNAseq_main <- RUN_SINGLEr_AVERAGE(
      so,
      fetch_celldex_ref("mouseRNAseq"),
      "label.main"
    )
    so$clustAnnot_mouseRNAseq <- RUN_SINGLEr_AVERAGE(
      so,
      fetch_celldex_ref("mouseRNAseq"),
      "label.fine"
    )
  }
  return(so)
}

#' Integration Report Functions
OBJECT_SELECT <- function(id) {
  obj <- switch(
    id,
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
  obj <- switch(
    id,
    "merged" = "Before Batch Correction",
    "integrated" = "Integrated CCA",
    "rpca" = "RPCA",
    "harmony" = "Harmony (Sample)",
    "scvi" = "single-cell Variational Inference",
    "liger" = "Linked Inference of Genomic Experimental Relationships"
  )
  return(obj)
}
