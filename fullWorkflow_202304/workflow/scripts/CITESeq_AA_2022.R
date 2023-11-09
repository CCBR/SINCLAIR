library(Seurat)
library(dsb,lib.loc = "/data/khanlab2/abdalla/mijjeCart/hashing")
args <- commandArgs(trailingOnly = TRUE)
rds <- as.character(args[1])
filtMatrix <- as.character(args[2])
rawMatrix <- as.character(args[3])



#read data 
so = readRDS(rds)
filtMatrix = Read10X_h5(filtMatrix)
rawMatrix <- Read10X_h5(rawMatrix)
rawCite = rawMatrix[[2]]
rawCite = rawCite[4:34,]
ctrl = filtMatrix[[2]][32:34,]
cite = filtMatrix[[2]][4:31,]



toSubtract = qlcMatrix::colMax(ctrl)
cite = sweep(cite,2,toSubtract)
cite[cite < 0] <- 0


cite = cite[,colnames(cite) %in% colnames(so)]
rownames(cite) = limma::strsplit2(rownames(cite),"_")[,1]

so[["CLR"]] <- CreateAssayObject(counts = cite)

so <- NormalizeData(so, assay = "CLR", normalization.method = "CLR")
so <- ScaleData(so, assay = "CLR")


SeuratObject = CreateSeuratObject(counts = rawMatrix[[1]])
SeuratObject = subset(SeuratObject, subset = nFeature_RNA < 40)


neg_adt_matrix = as.matrix(rawCite[,!(colnames(rawCite) %in% colnames(cite))])
neg_adt_matrix = neg_adt_matrix[,colnames(neg_adt_matrix) %in% colnames(SeuratObject)]

positive_adt_matrix = as.matrix(rawCite[,colnames(rawCite) %in% colnames(cite)])


normalized_matrix = DSBNormalizeProtein(cell_protein_matrix = positive_adt_matrix,
                                        empty_drop_matrix = neg_adt_matrix, isotype.control.name.vec = tail(rownames(rawCite),3))



normalized_matrix = normalized_matrix[,colnames(normalized_matrix) %in% colnames(so)]
normalized_matrix = normalized_matrix[1:28,]
rownames(normalized_matrix) = limma::strsplit2(rownames(normalized_matrix),"_")[,1]

so[["DSB"]] <- CreateAssayObject(data = normalized_matrix)
so@assays$cite = NULL


saveRDS(so,paste0("final/",rds))
