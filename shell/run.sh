#!/bin/bash
timestamp=$(date +%d%m%Y%H%M%S);
method=STAR;
echo $timestamp;
echo "export PATH=/STAR-2.7.9a/bin/Linux_x86_64_static:\$PATH" >> ~/.bashrc && \
. ~/.bashrc ;
perl ~/scvh_map_reads.pl -t 10 -r 32 -d viruSITE -a $method -o "s3://testscvh/output/SRR11537951$method$timestamp" \
"s3://testscvh/data/genome/" \
"s3://testscvh/data/SRR11537951part.0.R2.fastq.gz" \
"s3://testscvh/data/SRR11537951part.0.R1.fastq.gz" ;

Rscript ~/scvh_filter_matrix.r "s3://testscvh/output/SRR11537951$method$timestamp";
