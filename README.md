## <a name="require"></a>Requirements
* Input: 10x Chromium scRNA-seq reads
* STAR >= v2.7.8a
* cellranger >= 6.0.0
* Kallisto|bustools >= 0.25.1
* salmon|alevin >= v1.4.0
* DropletUtils >= v1.10.2

## <a name="gen_usages"></a>General usage
Map 10x scRNA-seq reads to the viral host reference set using STARsolo, CellRanger, Kallisto|bustools, and Salmon|Alevin.


```sh
Usage: scvh_map_reads.pl [Options] <vh_genome_dir> <R2> <R1>

Options:                                                                                                                
-o/--output-dir	<string>   the output directory                                                                                        
-t/--threads <int>         number of threads to run alignment with                                                                       
-d/--database <string>     select virus or virus and prokaryotes database, can be 'viruSITE' or 'viruSITE.NCBIprokaryotes'      
-e/--exe <string>          executable command or stand alone executable path of the alignment tool
-s/--soloStrand <string>   STARsolo param: Reverse or Forward used for 10x 5' or 3' protocol, respectively                               
-w/--whitelist <string>    STARsolo param --soloCBwhitelist                                                                            
-r/--ram <int>             Limitation of RAM usage. For STARsolo, param: limitGenomeGenerateRAM unit by GB 
-f/--soloFeatures <string> STARsolo param:  See --soloFeatures in STARsolo manual
-ot/--outSAMtype <string>  STARsolo param:  See --outSAMtype in STARsolo manual                                                          
-a/--alignment <string>    Select alignment methods: 'STAR', 'KB', 'Alevin', or 'CellRanger'                                             
-v/--technology <string>   KB param:  Single-cell technology used (`kb --list` to view)                                                  

```

For STARsolo, Kallisto|bustools, and Salmon|Alevin.

```sh
perl scvh_map_reads.pl -t num_threads -o output_dir vh_genome_dir R2.fastq.gz R1.fastq.gz
```

For CellRanger:

```sh
perl scvh_map_reads.pl -t num_threads -o output_dir vh_genome_dir sample fastqs
```

Filter the mapped UMIs using EmptyDrops to get the viral host filtered UMI counts matrix and also output viral genes and barcodes info files
```sh
Rscript scvh_filter_matrix.r output_dir sample_name
```

*where -t is a user-specified integer indicating number of threads to run with, output_dir is a user-specified directory to place the outputs, vh_genome_dir is a pre-generated viral host (human) reference set directory, R2.fastq.gz R1.fastq.gz are input 10x scRNA-seq reads, and sample_name is an optional user-specified tag to be used as a prefix for the output files.*
