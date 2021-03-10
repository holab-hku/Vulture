#!/bin/bash
echo $1
cellranger mkref \
--genome=human_host_viruses.viruSITE.with_hg38_10 \
--fasta=/home/junyichen/scvh_files/vh_genome_dir/human_host_viruses.viruSITE.with_hg38.fa \
--genes=/home/junyichen/scvh_files/vh_genome_dir/human_host_viruses.viruSITE.with_hg38.removed_amb_viral_exon.gtf \
--nthreads=$1;
mv *out output/mkref.log.out