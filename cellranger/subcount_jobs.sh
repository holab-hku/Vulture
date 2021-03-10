#!/bin/bash
data=$1
timestamp=$(date +%m%d%Y%H%M)
#rm -r human_host_viruses.viruSITE.with_hg38
/usr/bin/time -v nohup ./count_job.sh $data > output/$timestamp.out 2>&1