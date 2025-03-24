library(Signac)
library(Seurat)
library(GenomeInfoDb)

args <- commandArgs(trailingOnly = T)

so = readRDS(as.character(args[1]))
counts = Read10X_h5(as.character(args[2]))
fragpath = as.character(args[3])
annotation = readRDS("atachg38Annot.rds")


counts$Peaks = counts$Peaks[,colnames(so)]

# create ATAC assay and add it to the object
so[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)
DefaultAssay(so) <- "ATAC"

so <- NucleosomeSignal(so)
so <- TSSEnrichment(so)





# filter out low quality cells
so <- subset(
  x = so,
  subset = nCount_ATAC < 100000 &
 #   nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
 #   nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)

# call peaks using MACS2
peaks <- CallPeaks(so)

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(so),
  features = peaks,
  cells = colnames(so)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
so[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

DefaultAssay(so) <- "peaks"
so <- FindTopFeatures(so, min.cutoff = 5)
so <- RunTFIDF(so)
so <- RunSVD(so)


so <- FindMultiModalNeighbors(
  object = so,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:30, 2:30),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)
saveRDS(so,paste0("atacSoup/",as.character(args[1])))
