if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
install.packages("argparse")
BiocManager::install("DropletUtils",Ncpus=10)
BiocManager::install(c("Biostrings", "ShortRead","doParallel","GenomicAlignments","Gviz","GenomicFeatures","Rsubread"),Ncpus=10)
remotes::install_github("mojaveazure/seurat-disk",upgrade="always",Ncpus=10)
