#!/usr/bin/env nextflow
/*
 * pipeline input parameters
 */
params.baseDir = "."
params.codebase = "/code"
params.soloCBlen = 16;
params.soloCBstart = 1;
params.soloUMIstart = 17;
params.soloUMIlen = 10;
params.soloStrand = "Forward";
params.threads = 28;
params.ram = 64;
params.alignment = "STAR";
params.technology = "10XV2";
params.virus_database = "viruSITE";
params.pseudoBAM = "";
params.soloMultiMappers  = "EM";
params.soloFeatures = "Gene";
params.outSAMtype = "BAM SortedByCoordinate";
params.soloInputSAMattrBarcodeSeq = "CR UR";
params.soloInputSAMattrBarcodeQual = "-";


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
    ls -l ref
    cd ref 
    ls .
    cd ../
    perl ${params.codebase}/scvh_map_reads.pl -f GeneFull -t 24 -r 32 -d "viruSITE" -a "STAR" -o "." -ot "BAM Unsorted" \
    "${ref}" \
    "${params.baseDir}/${pair_id}_2.fastq.gz" "${params.baseDir}/${pair_id}_1.fastq.gz"

    """
    // """
    // perl ${params.codebase}/scvh_map_reads.pl \
    // --database ${params.virus_database} \
    // --threads ${params.threads} --ram ${params.ram} \
    // --alignment ${params.alignment} \
    // -o "." \
    // --whitelist ${params.barcodes_whitelist} \
    // --soloCBlen ${params.soloCBlen} --soloCBstart ${params.soloCBstart} \
    // --soloUMIstart ${params.soloUMIstart} --soloUMIlen ${params.soloUMIlen} \
    // --soloStrand ${params.soloStrand} \
    // --soloMultiMappers ${params.soloMultiMappers} \
    // --soloFeature ${params.soloFeatures} \
    // --outSAMtype ${params.outSAMtype} \
    // --technology ${params.technology} \
    // --pseudoBAM ${params.pseudoBAM} \
    // "${ref}" \
    // "${params.baseDir}/${pair_id}_R2.fastq.gz" "${params.baseDir}/${pair_id}_R1.fastq.gz";
    // rm ./alignment_outs/*.bam
    // """
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
