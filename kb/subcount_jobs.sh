#!/bin/bash
if [ "$#" -ne 6 ]; then
  echo "$#";
  echo "Usage: $0 <Fastq_R1_directory> <Fastq_R2_directory> <KBref_directory> <KBref_directory> <Output_directory> <threads> <RAM>" >&2
  exit 1
fi
timestamp=$(date +%m%d%Y%H%M);
#rm -r human_host_viruses.viruSITE.with_hg38
nohup /usr/bin/time -v ./kbref_job.sh $1 $2 $3 $4 $5 $6 > ../output/kbcount_$timestamp.out 2>&1