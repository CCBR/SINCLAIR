loc <- "/data/CCBR_Pipeliner/db/PipeDB/scrna4.2Rlibs"
library(Seurat, lib.loc = loc)
library("SingleR", lib.loc = loc)
library(scRNAseq, lib.loc = loc)
library(SingleCellExperiment, lib.loc = loc)
library(dplyr)
library(Matrix)
library(tools)
library(stringr)


args <- commandArgs(trailingOnly = TRUE)

matrix <- as.character(args[1])
outDirSeurat <- as.character(args[2])
ref <- as.character(args[3])
contrasts <- as.character(args[4])
resolutionString <- as.character(strsplit(gsub(",+", ",", as.character(args[5])), split = ",")[[1]])
resolution <- as.numeric(resolutionString)



options(future.globals.maxSize = 96000 * 1024^2)


file.names <- dir(path = matrix, pattern = "rds")

groupFile <- read.delim("groups.tab", header = F, stringsAsFactors = F)
if (grepl("-", contrasts, fixed = TRUE)) {
  groupFile <- groupFile[groupFile$V2 %in% stringr::str_split_fixed(contrasts, pattern = "-", n = Inf)[1, ], ]
  # print("test")
}

splitFiles <- gsub(".rds", "", file.names) # str_split_fixed(file.names,pattern = "[.rd]",n = 2)
splitFiles <- stringr::str_split_fixed(splitFiles, pattern = "__", n = Inf)[, 1]
file.names <- file.names[match(groupFile$V1, splitFiles, nomatch = F)]
print(groupFile$V1)
print(splitFiles)
print(file.names)


readObj <- list()
for (obj in file.names) {
  Name <- strsplit(obj, ".rds")[[1]][1]
  assign(paste0("S_", Name), readRDS(paste0(matrix, "/", obj)))
  readObj <- append(readObj, paste0("S_", Name))
}


combinedObj.list <- list()
i <- 1
for (p in readObj) {
  combinedObj.list[[p]] <- eval(parse(text = readObj[[i]]))
  combinedObj.list[[p]]$Sample <- names(combinedObj.list)[i]
  i <- i + 1
}


reference.list <- combinedObj.list[unlist(readObj)]

print(length(reference.list))
print(reference.list)

reference.list <- lapply(X = reference.list, FUN = SCTransform)
selectFeatures <- SelectIntegrationFeatures(object.list = reference.list, nfeatures = 3000)
reference.list <- PrepSCTIntegration(object.list = reference.list, anchor.features = selectFeatures, verbose = FALSE)
if (sum(sapply(combinedObj.list, ncol)) > 100000) {
  refIndex <- vector()
  for (group in unique(groupFile$V2)) {
    if (length(grep(group, groupFile$V2)) > 2) {
      refIndex <- c(refIndex, sample(grep(group, groupFile$V2), 1))
    }
  }
  if (length(refIndex) == 0) {
    refIndex <- sample(length(reference.list), 1)
  }
  print(refIndex)
  combinedObj.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30, anchor.features = 3000, reference = refIndex)
} else {
  combinedObj.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30, anchor.features = 3000)
}

combinedObj.integrated <- IntegrateData(anchorset = combinedObj.anchors, dims = 1:30)

combinedObj.integrated$Sample <- stringr::str_split_fixed(combinedObj.integrated$Sample, pattern = "__", n = Inf)[, 1]
combinedObj.integrated$groups <- groupFile$V2[match(combinedObj.integrated$Sample, paste0("S_", groupFile$V1), nomatch = F)]

DefaultAssay(object = combinedObj.integrated) <- "integrated"


npcs <- 50 # ifelse(ncol(combinedObj.integratedRNA) > 50,yes = 50, round(ncol(combinedObj.integratedRNA)/1000))


runInt <- function(obj, npcs) {
  if (ncol(obj@assays$integrated@scale.data) != ncol(obj@assays$integrated@data)) {
    obj <- ScaleData(obj)
  }
  obj <- RunPCA(object = obj, npcs = npcs, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:npcs)

  for (res in resolution) {
    obj <- FindClusters(obj, verbose = F, resolution = res, algorithm = 3) #
  }
  colnames(obj@meta.data) <- gsub("integrated_snn_res", "SLM_int_snn_res", colnames(obj@meta.data))
  colnames(obj@meta.data) <- gsub("^SCT_snn_res", "SLM_SCT_snn_res", colnames(obj@meta.data))


  for (res in resolution) {
    obj <- FindClusters(obj, verbose = F, resolution = res, algorithm = 4, method = "igraph") #
  }
  colnames(obj@meta.data) <- gsub("integrated_snn_res", "Leiden_int_snn_res", colnames(obj@meta.data))
  colnames(obj@meta.data) <- gsub("^SCT_snn_res", "Leiden_SCT_snn_res", colnames(obj@meta.data))


  obj <- RunUMAP(
    object = obj, reduction = "pca",
    dims = 1:npcs, n.components = 3
  )

  obj$groups <- groupFile$V2[match(gsub("S_", "", obj$Sample), groupFile$V1, nomatch = F)]

  runSingleR <- function(obj, refFile, fineORmain) {
    avg <- AverageExpression(obj, assays = "SCT")
    avg <- as.data.frame(avg)
    ref <- refFile
    s <- SingleR(test = as.matrix(avg), ref = ref, labels = ref[[fineORmain]])

    clustAnnot <- s$labels
    names(clustAnnot) <- colnames(avg)
    names(clustAnnot) <- gsub("SCT.", "", names(clustAnnot))

    obj$clustAnnot <- clustAnnot[match(obj$seurat_clusters, names(clustAnnot))]
    return(obj$clustAnnot)
  }

  if (ref == "hg38") {
    obj$clustAnnot_HPCA_main <- runSingleR(obj, celldex::HumanPrimaryCellAtlasData(), "label.main")
    obj$clustAnnot_HPCA <- runSingleR(obj, celldex::HumanPrimaryCellAtlasData(), "label.fine")
    obj$clustAnnot_BP_encode_main <- runSingleR(obj, celldex::BlueprintEncodeData(), "label.main")
    obj$clustAnnot_BP_encode <- runSingleR(obj, celldex::BlueprintEncodeData(), "label.fine")
    obj$clustAnnot_monaco_main <- runSingleR(obj, celldex::MonacoImmuneData(), "label.main")
    obj$clustAnnot_monaco <- runSingleR(obj, celldex::MonacoImmuneData(), "label.fine")
    obj$clustAnnot_immu_cell_exp_main <- runSingleR(obj, celldex::DatabaseImmuneCellExpressionData(), "label.main")
    obj$clustAnnot_immu_cell_exp <- runSingleR(obj, celldex::DatabaseImmuneCellExpressionData(), "label.fine")
  }

  if (ref == "mm10") {
    obj$clustAnnot_immgen_main <- runSingleR(obj, celldex::ImmGenData(), "label.main")
    obj$clustAnnot_immgen <- runSingleR(obj, celldex::ImmGenData(), "label.fine")
    obj$clustAnnot_mouseRNAseq_main <- runSingleR(obj, celldex::MouseRNAseqData(), "label.main")
    obj$clustAnnot_mouseRNAseq <- runSingleR(obj, celldex::MouseRNAseqData(), "label.fine")
  }
  return(obj)
}

if (contrasts == "None") {
  saveRDS(combinedObj.integrated, "intermediate_integrated.rds")
}
combinedObj.integrated <- runInt(combinedObj.integrated, npcs)

saveRDS(combinedObj.integrated, outDirSeurat)
