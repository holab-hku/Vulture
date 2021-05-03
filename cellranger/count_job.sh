#!/bin/bash
echo $1
cellranger count \
--id=run_covid_vmhost \
--fastqs=/home/jholab/data/COV19/SRR11537949 \
--sample=SRR11537949 \
--localcores = 10 \
--include-introns \
--transcriptome=/home/jholab/data/vmh_genome_dir/references/human_host_viruses_microbes.viruSITE.NCBIprokaryotes.with_hg38 & \
cellranger count \
--id=run_covid_vmhost \
--fastqs=/home/jholab/data/COV19/SRR11537950 \
--sample=SRR11537950 \
--localcores = 10 \
--include-introns \
--transcriptome=/home/jholab/data/vmh_genome_dir/references/human_host_viruses_microbes.viruSITE.NCBIprokaryotes.with_hg38 & \
cellranger count \
--id=run_covid_vmhost \
--fastqs=/home/jholab/data/COV19/SRR11537951 \
--sample=SRR11537951 \
--localcores = 10 \
--include-introns \
--transcriptome=/home/jholab/data/vmh_genome_dir/references/human_host_viruses_microbes.viruSITE.NCBIprokaryotes.with_hg38;