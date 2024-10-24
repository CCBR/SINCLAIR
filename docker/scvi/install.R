#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
sys.setenv(GITHUB_PAT = args[1])

setRepositories(ind = 1:3)
remotes::install_github("satijalab/seurat-wrappers")
remotes::intall_cran("BPCells", repos = c("https://bnprks.r-universe.dev"))
