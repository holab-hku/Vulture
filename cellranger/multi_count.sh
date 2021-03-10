#!/bin/bash
cellranger count \
--id=output_SRR11537951 \
--fastqs=/storage/holab/covid_scRNA-seq/SRR11537951 \
--sample=SRR11537951 \
--transcriptome=human_host_viruses.viruSITE.with_hg38 & \
cellranger count \
--id=output_SRR11537950 \
--fastqs=/storage/holab/covid_scRNA-seq/SRR11537950 \
--sample=SRR11537950 \
--transcriptome=human_host_viruses.viruSITE.with_hg38 & \
cellranger count \
--id=output_SRR11537949 \
--fastqs=/storage/holab/covid_scRNA-seq/SRR11537949 \
--sample=SRR11537949 \
--transcriptome=human_host_viruses.viruSITE.with_hg38