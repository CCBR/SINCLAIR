.libPaths(c("/data/CCBR_Pipeliner/db/PipeDB/scrna5", .libPaths()))


library(Seurat)
library(stringr)
library(reticulate)
library(SingleR)
library(scRNAseq)
library(SingleCellExperiment)
library(celldex)
library(Orthology.eg.db)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(SeuratWrappers)

source("source.R")
args <- commandArgs(trailingOnly = T)


so <- as.character(args[1])
outDirSeurat <- as.character(args[2])
ref <- as.character(args[3])
pcs <- as.character(args[4])
resolutions <- as.character(strsplit(gsub(",+", ",", as.character(args[5])), split = ",")[[1]])
resolutions <- as.numeric(resolutions)
so <- readRDS(so)


saveRDS(so, outDirSeurat)
