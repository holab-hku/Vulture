if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DropletUtils")
BiocManager::install(c("Biostrings", "ShortRead","doParallel","GenomicAlignments","Gviz","GenomicFeatures","Rsubread"))
