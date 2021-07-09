#!/bin/bash
timestamp=$(date +%d%m%Y%H%M%S);
method=STAR;
echo $timestamp;
echo "export PATH=/STAR-2.7.9a/bin/Linux_x86_64_static:\$PATH" >> ~/.bashrc && \
. ~/.bashrc ;
OUT_DIR="/storage/output/scvh_out_$method\_$timestamp";
PARAMS=$(cat "~/args.txt")
perl ~/scvh_map_reads.pl -o $OUT_DIR $PARAMS;
Rscript ~/scvh_filter_matrix.r $OUT_DIR;
