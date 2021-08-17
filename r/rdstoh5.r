library(Seurat)
library(SeuratDisk)
args <- commandArgs(trailingOnly = TRUE)
rds_file <- args[1]
h5ad_file <- args[2]

filtered_matrix <- readRDS(rds_file)
data.seurat <- CreateSeuratObject(filtered_matrix)
SaveH5Seurat(data.seurat,filename=paste(h5ad_file,".h5Seurat",sep=""))
Convert(paste(h5ad_file,".h5Seurat",sep=""),dest="h5ad")
