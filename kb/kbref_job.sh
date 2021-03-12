#!/bin/bash
echo fasta path $1;
echo gtf path $2;
echo output path $3;
kb ref \
-i=$3 \
-g=$3 \
$1 $2