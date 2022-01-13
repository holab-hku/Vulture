#!/usr/bin/env nextflow
/*
 * pipeline input parameters
 */
params.baseDir = "."
params.reads = "s3://sra-pub-sars-cov2/sra-src/SRR11537951/*_R{2,1}.fastq.gz" 
params.ref = "s3://scvhwf/prod/ref/STAR"
params.outdir = "s3://scvhwf/prod/output/demo"
params.codebase = "~"
log.info """\
         S C V H - N F   P I P E L I N E
         ===================================
         transcriptome: ${params.ref}
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

    publishDir "${params.outdir}", mode: "copy"
    cpus 30
    memory '80 GB'

    input:
    tuple val(pair_id), file(reads) from read_pairs_ch
    path ref from params.ref
    output:
    set file("intermediate_files/"), file('alignment_outs/') into results_ch
    
    shell
    """
    source ~/.bashrc
    perl ${params.codebase}/scvh_map_reads.pl \
    --virus_database ${params.virus_database} \
    --threads ${params.threads} --ram ${params.ram} \
    --alignment ${params.alignment} \
    -o "." \
    --barcodes_whitelist ${params.barcodes_whitelist} \
    --soloCBlen ${params.soloCBlen} --soloCBstart ${params.soloCBstart} \
    --soloUMIstart ${params.soloUMIstart} --soloUMIlen ${params.soloUMIlen} \
    --soloStrand ${params.soloStrand} \
    --soloMultiMappers ${params.soloMultiMappers} \
    --soloFeatures ${params.soloFeatures} \
    --outSAMtype ${params.outSAMtype} \
    --soloInputSAMattrBarcodeSeq ${params.soloInputSAMattrBarcodeSeq} \
    --soloInputSAMattrBarcodeQual ${params.soloInputSAMattrBarcodeQual} \ 
    --technology ${params.technology} \
    --pseudoBAM ${params.pseudoBAM} \
    "${ref}" \
    "${params.baseDir}/${pair_id}_R2.fastq.gz" "${params.baseDir}/${pair_id}_R1.fastq.gz";
    rm ./alignment_outs/*.bam
    """
}
/*
 * 2. Filter
 */
process Filter {
    publishDir "${params.outdir}", mode: "copy"
    
    cpus 8
    memory '32 GB'
    input:
    file result from results_ch
    output:
    set file('*.txt'), file('*.rds') into results2_ch
    
    shell
    """
    Rscript ${params.codebase}/scvh_filter_matrix.r "."
    """
}
