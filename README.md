## <a name="require"></a>Requirements
* Input: 10x Chromium scRNA-seq reads
* STAR >= v2.7.5a
* DropletUtils >= v1.10.2

## <a name="gen_usages"></a>General usage
Map 10x scRNA-seq reads to the viral host reference set using STARsolo
```sh
perl scvh_map_reads.pl -t num_threads -o output_dir vh_genome_dir R2.fastq.gz R1.fastq.gz
```

Filter the mapped UMIs using EmptyDrops to get the viral host filtered UMI counts matrix and also output viral genes and barcodes info files
```sh
Rscript scvh_filter_matrix.r output_dir sample_name
```

*where -t is a user-specified integer indicating number of threads to run STARsolo with, output_dir is a user-specified directory to place the outputs, vh_genome_dir is a pre-generated viral host (human) reference set directory, R2.fastq.gz R1.fastq.gz are input 10x scRNA-seq reads, and sample_name is an optional user-specified tag to be used as a prefix for the output files.*
