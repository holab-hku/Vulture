#!/usr/bin/env nextflow
/*
 * pipeline input parameters
 */
params.reads = "s3://scvhwf/samples/storage/holab/covid_scRNA-seq/*_{2,1}.fastq.gz" 
params.outdir = "s3://scvhwf/prod/output/demo"
params.annotation = "s3://scvhwf/prod/ref"
params.codebase = "/code"
params.baseDir = "."
log.info """\
         SCVH - N F   P I P E L I N E - MULTI
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
    input:
    tuple val(pair_id), file(reads) from read_pairs_ch
    path ref from params.annotation
    output:
    val pair_id into id_ch
    set file("intermediate_files/"), file("alignment_outs/") into results_ch

    shell
    """
    source ~/.bashrc
    df -h
    pwd
    perl ${params.codebase}/scvh_map_reads.pl -f GeneFull -t 24 -r 32 -d "viruSITE" -a "STAR" -o "." -ot "BAM Unsorted" \
    "${ref}" \
    "${params.baseDir}/${pair_id}_2.fastq.gz" "${params.baseDir}/${pair_id}_1.fastq.gz"
    """
}





/*
 * 2. Filter
 */
process Filter {
    publishDir "${params.outdir}/${pair_id}", mode: "copy"
    input:
    file result from results_ch
    val pair_id from id_ch
    output:
    set file('*.txt'), file('*.rds'), file("alignment_outs/Solo.out/") into results2_ch

    shell
    """
    ls
    pwd
    Rscript ${params.codebase}/scvh_filter_matrix.r "."
    """
}