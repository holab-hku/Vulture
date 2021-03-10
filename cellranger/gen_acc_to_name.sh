#!/bin/bash
echo Reference path is:
echo $1
gunzip -k $1/genes/genes.gtf.gz
ln -s $PWD/$1/genes/genes.gtf $1/star/$1.removed_amb_viral_exon.gtf
ln -s $PWD/$1/fasta/genome.fa $1/star/$1.fa
../scvh_get_accesion_and_name.pl $1

# echo Result path is:
# echo $2 
# gunzip -k $2/outs/raw_feature_bc_matrix/*gz

