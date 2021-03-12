library(DropletUtils)

sample_name <- ""
args <- commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  output_dir <- "."
} else if (length(args)==1) {
  output_dir <- args[1]
} else if (length(args)==2) {
	output_dir <- args[1]
	sample_name <- args[2]
} else {
	stop("Input a filltered matirx. A maximum of two arguments, the first being the output directory and the second being the sample_name, can be supplied")
}

#filter for lowest total UMI count per cell
lower <- 100

host_gene_id_prefix <- "ENSG"

check_dirs_files_exist <- function(output_dir) {
	
	scvh_map_reads_script_name <- "scvh_map_reads.pl"
	gene_to_accession_and_name_file_name <- "gene_to_accession_and_name.txt"
	
	intermediate_files_dir <- paste0(output_dir,"/intermediate_files")
	if (dir.exists(intermediate_files_dir) == FALSE) {
		stop(paste0("Can't find ", intermediate_files_dir,", please run ",scvh_map_reads_script_name))
	}
	
	gene_to_accession_and_name_path <- paste0(intermediate_files_dir,"/",gene_to_accession_and_name_file_name)
	if (file.exists(gene_to_accession_and_name_path) == FALSE) {
		stop(paste0("Can't find ",gene_to_accession_and_name_path,", please run ",scvh_map_reads_script_name))
	}
	
	STARsolo_dir_name <- "STARsolo_outs"
	STARsolo_path_to_raw_matrix <- "Solo.out/Gene/filtered"
	STARsolo_path_to_raw_matrix_path <- paste0(output_dir,"/",STARsolo_dir_name,"/",STARsolo_path_to_raw_matrix)
	
	if (dir.exists(paste0(output_dir,"/",STARsolo_dir_name)) == FALSE) {
		stop(paste0("Can't find ",output_dir,"/",STARsolo_dir_name,", please run",scvh_map_reads_script_name))
	}
	
	if (dir.exists(STARsolo_path_to_raw_matrix_path) == FALSE) {
		stop(paste0("Can't find ",STARsolo_path_to_raw_matrix_path,", make sure STARsolo output directory structure is unchanged"))
	}
	
	return(STARsolo_path_to_raw_matrix_path)
}

STARsolo_filtered_matrix <- check_dirs_files_exist(output_dir)
gene_to_accession_and_name_path <- paste0(output_dir,"/intermediate_files/gene_to_accession_and_name.txt")
gene_to_accession_and_name <- read.delim(gene_to_accession_and_name_path, sep = "\t", header = T)
#make first column (genes) as row names
rownames(gene_to_accession_and_name) <- gene_to_accession_and_name[,1]
gene_to_accession_and_name[,1] <- NULL

gzip_command <- paste0("gzip -r ",STARsolo_filtered_matrix)
system(gzip_command)

#read 10x data as SingleCellExperiment
filtered_SCE <- read10xCounts(paste0(STARsolo_filtered_matrix), version = "3")

filtered_matrix <- counts(filtered_SCE)
colnames(filtered_matrix) <- filtered_SCE@colData@listData$Barcode

#viral genes
viral_genes_indices <- grep(paste0("^",host_gene_id_prefix), rownames(filtered_matrix), invert = TRUE)

#human genes
human_genes_indices <- grep(paste0("^",host_gene_id_prefix), rownames(filtered_matrix))

#gene names at
#filtered_SCE@rowRanges@elementMetadata@listData$Symbol
rownames(filtered_matrix) <- filtered_SCE@rowRanges@elementMetadata@listData$Symbol

filtered_matrix_viral <- filtered_matrix[viral_genes_indices,]
filtered_matrix_human <- filtered_matrix[human_genes_indices,]

print_viral_gene_UMIs_and_barcode_UMIs <- function(matrix_viral, matrix_human, set) {
	
	viral_UMIs_table <- ""
		
	for (element in c("genes", "barcodes")) {
		if (element == "genes") {
			sum_UMI_viral <- rowSums(matrix_viral)
		} else if (element == "barcodes") {
			sum_UMI_viral <- colSums(matrix_viral)
		}
		
		sum_UMI_viral <- sort(sum_UMI_viral, decreasing = TRUE)
		
		if (max(sum_UMI_viral) == 0) {
			print(paste0("No ",set," viral UMIs"))
			break
		}
		
		viral_UMIs <- sum_UMI_viral[which(sum_UMI_viral > 0)]
		
		if (set == "unfiltered") {
			out_dir <- paste0(output_dir,"/intermediate_files")
		} else if (set == "filtered") {
			out_dir <- output_dir
		} else {
			stop("Unrecognized set ",set)
		}
		
		if (sample_name != "") {
			sample_name_tag <- paste0(sample_name,"_")
		} else {
			sample_name_tag <- ""
		}
		
		if (element == "genes") {
			median_num_human_genes_expressed <- c()
			viral_genes <- names(viral_UMIs)
			for (gene in viral_genes) {
				gene_matrix_viral <- matrix_viral[gene,]
				gene_matrix_viral_expressed_barcodes <- gene_matrix_viral[gene_matrix_viral > 0]
				expressed_barcodes <- names(gene_matrix_viral_expressed_barcodes)
				
				matrix_human_viral_gene_expressed_barcodes <- matrix_human[,expressed_barcodes, drop = FALSE]
				viral_gene_expressed_barcodes_num_human_genes_expressed <- colSums(matrix_human_viral_gene_expressed_barcodes != 0)
				viral_gene_expressed_barcodes_median_num_human_genes_expressed <- median(viral_gene_expressed_barcodes_num_human_genes_expressed)
				
				median_num_human_genes_expressed <- c(median_num_human_genes_expressed, viral_gene_expressed_barcodes_median_num_human_genes_expressed)
			}
			#paste0(set,"_UMI_counts") = as.numeric(viral_UMIs)
			UMI_counts <- paste0(set,"_UMI_counts")
			
			accession_and_names <- gene_to_accession_and_name[names(viral_UMIs), ]
			
			viral_UMIs_table <- data.frame("gene" = names(viral_UMIs), UMI_counts = as.numeric(viral_UMIs), "median_num_human_genes_expressed" = median_num_human_genes_expressed, "accession" = accession_and_names$accession, "reference_name" = accession_and_names$reference_name)
			
			write.table(viral_UMIs_table, paste0(out_dir,"/",sample_name_tag,set,"_matrix_viral_",element,"_info.txt"), quote=F, row.names=F, sep = "\t")
			
		} else if (element == "barcodes") {
			# Fix bug by JY 2021 0309
			viral_UMIs_table[is.na(viral_UMIs_table)]<-"unknown"
			out_matrix_viral <- matrix_viral[viral_UMIs_table$gene,names(viral_UMIs), drop = FALSE]
			out_matrix_viral <- as.data.frame(as.matrix(t(out_matrix_viral)))
			colnames(out_matrix_viral) <- viral_UMIs_table$accession
			out_viral_by_accession <- t(rowsum(t(out_matrix_viral), group = colnames(out_matrix_viral)))
			
			#for reordering columns
			accession_totals <- colSums(out_viral_by_accession)
			descending_accessions <- names(sort(accession_totals, decreasing = TRUE))
			out_viral_by_accession <- out_viral_by_accession[,descending_accessions, drop = FALSE]
			
			write.table(out_viral_by_accession, paste0(out_dir,"/",sample_name_tag,set,"_matrix_viral_",element,"_info.txt"), quote=F, col.names = NA, sep = "\t")
		}
		
	}
}
print_viral_gene_UMIs_and_barcode_UMIs(filtered_matrix_viral, filtered_matrix_human, "filtered")

if (sample_name != "") {
	sample_name_tag <- paste0(sample_name,"_")
} else {
	sample_name_tag <- ""
}
##Save RDS of filtered matrix
saveRDS(filtered_matrix, paste0(output_dir,"/",sample_name_tag,"filtered_matrix.rds"))