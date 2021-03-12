#!/bin/bash
if [ "$#" -ne 2 ]; then
  echo "$#";
  echo "Usage: $0 <Fasta_directory> <output_directory>" >&2
  exit 1
fi
timestamp=$(date +%m%d%Y%H%M);
fastapath=$1;
outpath=$2;
#rm -r human_host_viruses.viruSITE.with_hg38
nohup /usr/bin/time -v ./kbref_job.sh $fastapath $outpath > output/kbref_$timestamp.out 2>&1