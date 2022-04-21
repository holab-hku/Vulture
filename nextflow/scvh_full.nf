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
params.dumpT = 12;
params.soloStrand = "Forward";
params.threads = 28;
params.ram = 64;
params.alignment = "STAR";
params.technology = "10XV2";
params.virus_database = "viruSITE";
params.pseudoBAM = "";
params.soloMultiMappers = "EM";
params.soloFeatures = "Gene";
params.outSAMtype = "BAM SortedByCoordinate";
params.soloInputSAMattrBarcodeSeq = "CR UR";
params.soloInputSAMattrBarcodeQual = "-";
params.barcodes_whitelist = "737K-august-2016.txt"


log.info """\
         S C V H - N F   P I P E L I N E
         ===================================
         transcriptome: ${params.ref}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         database:    : ${params.virus_database}
         threads      : ${params.threads} 
         ram          : ${params.ram} 
         alignment    : ${params.alignment} 
         whitelist    : ${params.barcodes_whitelist} 
         soloCBlen    : ${params.soloCBlen} 
         soloCBstart  : ${params.soloCBstart} 
         soloUMIstart : ${params.soloUMIstart} 
         soloUMIlen   : ${params.soloUMIlen} 
         soloStrand   : ${params.soloStrand} 
         soloMultiMappers: ${params.soloMultiMappers} 
         soloFeature : ${params.soloFeatures} 
         outSAMtype   : ${params.outSAMtype} 
         technology   : ${params.technology} 
         pseudoBAM    : ${params.pseudoBAM} 
         """
         .stripIndent()


Channel
    .fromList( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .view()
    .set { read_pairs_ch }

/*
 * 2. Dump
 */
process Dump {
    cpus 16
    memory '32 GB'
    errorStrategy 'retry'
    maxRetries 3
    input:
    val pair_id from read_pairs_ch
    output:
    set file('*_1.fastq.gz'), file('*_2.fastq.gz') into results3_ch
    val pair_id into id_ch3

    shell
    """
    echo "Start to dump fastq files from .sra files: ${pair_id}";
    fasterq-dump --split-files -e ${params.dumpT} ${pair_id};
    echo "Dump finished";
    pigz -p 16 ${pair_id}_1.fastq ;
    pigz -p 16 ${pair_id}_2.fastq ;
    echo "Compression finished";
    """
}
    
/*
 * 4. Mapping 
 */
process Map {

    publishDir "${params.outdir}", mode: "copy"
    cpus 30
    memory '80 GB'

    input:
    file result from results3_ch
    val pair_id from id_ch3
    path ref from params.ref
    output:
    set file("intermediate_files/"), file('alignment_outs/') into results4_ch
    
    shell
    // """

    // ls -l ref
    // cd ref 
    // ls .
    // cd ../
    // perl ${params.codebase}/scvh_map_reads.pl -f GeneFull -t 24 -r 32 -d "viruSITE" -a "STAR" -o "." -ot "BAM Unsorted" \
    // "${ref}" \
    // "${params.baseDir}/${pair_id}_2.fastq.gz" "${params.baseDir}/${pair_id}_1.fastq.gz"

    // """
    """
    perl ${params.codebase}/scvh_map_reads.pl \
    --output-dir "." \
    --threads ${params.threads} \
    --ram ${params.ram} \
    --database ${params.virus_database} \
    --soloStrand ${params.soloStrand} \
    --whitelist "${ref}/${params.barcodes_whitelist}" \
    --alignment ${params.alignment} \
    --technology ${params.technology} \
    --soloFeature ${params.soloFeatures} \
    --soloMultiMappers ${params.soloMultiMappers} \
    --pseudoBAM ${params.pseudoBAM} \
    --soloCBlen ${params.soloCBlen} \
    --soloCBstart ${params.soloCBstart} \
    --soloUMIstart ${params.soloUMIstart} \
    --soloUMIlen ${params.soloUMIlen} \
    "${ref}" \
    "${params.baseDir}/${pair_id}_2.fastq.gz" \
    "${params.baseDir}/${pair_id}_1.fastq.gz";
    rm ./alignment_outs/*.bam
    """
}
/*
 * 5. Filter
 */
process Filter {
    publishDir "${params.outdir}", mode: "copy"
    
    cpus 8
    memory '32 GB'
    input:
    file result from results4_ch
    output:
    set file('*.txt'), file('*.rds') into results5_ch
    
    shell
    """
    Rscript ${params.codebase}/scvh_filter_matrix.r "."
    perl ${params.codebase}/scvh_analyze_bam.pl "."
    """
}
