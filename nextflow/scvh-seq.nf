#!/usr/bin/env nextflow
/*
 * pipeline input parameters
 */
params.reads = "s3://scvh-input-data/fastq/*_{2,1}.fastq.gz"
params.outdir = "ouput_data"
params.annotation = "s3://scvh-input-data/ref"
params.codebase = "/code"
params.baseDir = "."
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
    cpus 16
    memory '64 GB'

    input:
    tuple val(pair_id), file(reads) from read_pairs_ch
    path ref from params.annotation
    output:
    set file("intermediate_files/"), file('alignment_outs/') into results_ch
    
    shell
    """
    ls -l
    source ~/.bashrc
    ls ref
    perl ${params.codebase}/scvh_map_reads.pl -t 16 -r 32 -d "viruSITE" -a "STAR" -o "." \
    "${ref}" \
    "${params.baseDir}/${pair_id}_2.fastq.gz" "${params.baseDir}/${pair_id}_1.fastq.gz"
    ls
    ls alignment_outs
    """
}

/*
 * 2. Filter
 */
process Filter {
    publishDir "${params.outdir}", mode: "copy"
    input:
    file result from results_ch
    output:
    set file('*.txt'), file('*.rds') into results2_ch
    
    shell
    """
    Rscript ${params.codebase}/scvh_filter_matrix.r "."
    """
}

