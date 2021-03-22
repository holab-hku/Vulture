library(Matrix)
args <- commandArgs(trailingOnly = TRUE)
filtered_matrix_file <- args[1]
genes <- args[2]

filtered_matrix <- readRDS(filtered_matrix_file)


genes_splitted <- strsplit(genes, ",")

for (gene in genes_splitted[[1]]) {
	filtered_matrix_gene <- filtered_matrix[gene, ,drop=FALSE]

	filtered_matrix_gene <- filtered_matrix_gene[, colSums(filtered_matrix_gene != 0) > 0, drop = FALSE]

	cat(gene,colnames(filtered_matrix_gene))
	cat("\n")
	
}