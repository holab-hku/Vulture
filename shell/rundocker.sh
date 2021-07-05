#!/bin/bash
docker run -v /home/jholab/drive:/storage scvh /bin/sh "/storage/run.sh"
#nohup python fastq_splitter.py -o /home/junyichen/scvh_splitted_fq -i /home/junyichen/scvh_files/10x_scRNA-seq_fastq/Covid-19/SRR11537950_2.fastq.gz > split2.out 2>&1
