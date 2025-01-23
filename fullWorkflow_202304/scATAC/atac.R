loc = "/data/CCBR_Pipeliner/db/PipeDB/scrna4.2Rlibs"
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86,lib.loc=loc)
library(EnsDb.Mmusculus.v79,lib.loc=loc)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Seurat,lib.loc=loc)
library(Signac,lib.loc=loc)


args <- commandArgs(trailingOnly = T)

h5 = as.character(args[1])
rnaCounts = Read10X_h5(h5)$`Gene Expression`
atacCounts = Read10X_h5(h5)$Peaks

print(head(colnames(atacCounts)))

fragments = as.character(args[2])
ref =  as.character(args[3])

print(ref)

if (ref == "hg38") {
  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotation) <- "UCSC"
#  genome(annotations) <- ref
}

if (ref == "mm10") {
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- ref
}

so = CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

#so <- CreateSeuratObject(counts = atacCounts)
print("xxx")
so["ATAC"] <- CreateChromatinAssay(
  counts = atacCounts,
  sep = c("-", "-"),
  genome = ref,
  sep = c(":", "-"),
  annotation = annotations
)

#Annotation(so[["ATAC"]]) <- annotations
