#!/bin/bash
threads=$1
timestamp=$(date +%m%d%Y%H%M)
#rm -r human_host_viruses.viruSITE.with_hg38
/usr/bin/time -v nohup ./mkref_job.sh $threads > output/$timestamp.out 2>&1