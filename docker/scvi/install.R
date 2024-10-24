setRepositories(ind = 1:3, addURLs = c("https://satijalab.r-universe.dev", "https://bnprks.r-universe.dev/"))
install.packages(c("BPCells", "glmGamPoi", "Signac"))
remotes::install_github("satijalab/seurat-wrappers")
