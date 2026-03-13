install.packages("remotes", "devtools", "BiocManager", "pak")
devtools::install_github("satijalab/seurat", ref = "seurat5")
library(pak)
seuratPkgList <- gsub(".*/", "", pak::pkg_deps("Seurat")$ref)
seuratPkgList <- seuratPkgList[-which(SeuratPkgList == "Seurat")]
install.packages(seuratPkgList)
install.packages("fastDummies", "tinytex")
BiocManager::install(c(
  "AnnotationDbi", "AnnotationFilter", "AnnotationHub", "beachmat",
  "Biobase", "BiocFileCache", "BiocIO", "BiocNeighbors", "BiocParallel", "BiocSingular", "BiocVersion", "Biostrings", "BiocGenerics",
  "biomaRt", "bit", "bit64", "blob", "celldex", "DBI", "dbplyr", "DelayedArray", "DelayedMatrixStats", "dotCall64", "DT",
  "ensembldb", "ExperimentHub", "fields", "filelock", "flexmix", "formatR", "futile.logger", "futile.options",
  "GenomeInfoDb", "GenomeInfoDbData", "GenomicAlignments", "GenomicFeatures", "GenomicRanges", "hdf5r", "hms",
  "interactiveDisplayBase", "IRanges", "KEGGREST", "lambda.r", "maps", "MatrixGenerics", "modeltools", "nnet",
  "Orthology.eg.db", "plogr", "prettyunits", "progress", "ProtGenerics", "RcppHNSW", "RCurl", "restfulr", "Rhtslib",
  "rjson", "R.methodsS3", "R.oo", "Routliers", "Rsamtools", "RSpectra", "RSQLite", "rsvd", "rtracklayer", "R.utils",
  "S4Arrays", "S4Vectors", "ScaledMatrix", "scRNASeq", "SeuratWrappers", "SingleCellExperiment", "SingleR",
  "snow", "spam", "sparseMatrixStats", "SummarizedExperiment", "viridis", "XML", "xml2", "XVector", "zlibbioc"
))
devtools::install_github("chris-mcginnis-ucsf/DoubletFinder")
devtools::install_github("satijalab/seurat-wrappers", ref = "seurat5")
devtools::install_github("satijalab/seurat-data", ref = "seurat5")
devtools::install_github("satijalab/azimuth", ref = "seurat5")
devtools::install_github("stuart-lab/signac", ref = "seurat5")
devtools::install_github("immunogenomics/harmony")
