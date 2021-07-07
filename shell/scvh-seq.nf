#!/usr/bin/env nextflow
/*
 * pipeline input parameters
 */
params.reads = "/home/jholab/users/angelayin/input_data/fastq/*_{2,1}.fastq.gz"
params.outdir = "ouput_data"
params.annotation = "/home/jholab/users/angelayin/input_data/ref"
params.codebase = "/home/jholab/users/angelayin/dev/scvh"
params.baseDir = "/home/jholab/users/angelayin/input_data/fastq"
log.info """\
         SCVH - N F   P I P E L I N E
         ===================================
         transcriptome: ${params.annotation}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()


Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .view()
    .set { read_pairs_ch }

/*
 * 1. Mapping 
 */
process Map {
    publishDir "${params.outdir}", mode: "copy", overwrite: "true"
    input:
    tuple val(pair_id), file(reads) from read_pairs_ch
    output:
    set file("intermediate_files/"), file('alignment_outs/') into results_ch
    
    shell
    """
    perl ${params.codebase}/scvh_map_reads.pl -t 10 -r 32 -d "viruSITE" -a "STAR" -o "." \
    "${params.annotation}" \
    "${params.baseDir}/${pair_id}_2.fastq.gz" "${params.baseDir}/${pair_id}_1.fastq.gz"
    """
}

/*
 * 2. Filter
 */
process Filter {
    publishDir "${params.outdir}", mode: "copy"
    input:
    file result from results_ch
    
    shell
    """
    Rscript ${params.codebase}/scvh_filter_matrix.r "${params.outdir}"
    """
}

