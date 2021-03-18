#!/bin/bash
if [ "$#" -ne 3 ]; then
  echo "$#";
  echo "Usage: $0 <Fasta_directory> <gtf_directory> <output_directory>" >&2
  exit 1
fi
timestamp=$(date +%m%d%Y%H%M);
fastapath=$1;
gtfpath=$2;
outpath=$3;

#rm -r human_host_viruses.viruSITE.with_hg38
nohup /usr/bin/time -v ./kb/kbref_job.sh $fastapath $gtfpath $outpath > ./output/kbref_$timestamp.out 2>&1