#!/bin/bash
echo fasta path $1;
echo gtf path $2;
echo output path $3;
kb ref \
-i=$3/transcriptome.idx \
-g=$3/transcripts_to_genes.txt \
-f1=$3/cdna.fa \
$1 $2