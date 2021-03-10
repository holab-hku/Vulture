#!/bin/bash
echo $1
cellranger count \
--id=run_covid_host \
--fastqs=/storage/holab/covid_scRNA-seq/$1 \
--sample=SRR11537949 \
--transcriptome=human_host_viruses.viruSITE.with_hg38
