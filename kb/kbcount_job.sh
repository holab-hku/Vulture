#!/bin/bash
echo fastq path 1 $1;
echo fastq path 2 $2;
echo kbref path 2 $3;
echo output path $4;
echo threads $5;
echo RAM $6;
kb count \
-x=10XV2 \
-g=$3/transcripts_to_genes.txt \
-i=$3/transcriptome.idx \
-o=$3 \
-t=$5 \
-m=$6 \
$1 $2