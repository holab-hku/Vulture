#!/usr/bin/env nextflow
/*
 * pipeline input parameters
 */
params.baseDir = "."
params.codebase = "/code"
params.dumpT = 12;
params.soloStrand = "Forward";
params.threads = 10;
params.ram = 128;
params.alignment = "STAR";
params.virus_database = "viruSITE";
params.pseudoBAM = "";
params.soloMultiMappers = "EM";
params.soloFeatures = "GeneFull";
params.outSAMtype = "BAM SortedByCoordinate";
params.soloInputSAMattrBarcodeQual = "-";
params.technology = "10XV2";
params.inputformat = "fastq";
params.sampleSubfix1 = "_1";
params.sampleSubfix2 = "_2";
params.barcode = "CR";
params.umi = "UR";
params.soloInputSAMattrBarcodeSeq = "${params.barcode} ${params.umi}";
params.mapqueue = "jy-scvh-queue-r5a8x-1"
params.downloadqueue = "jy-scvh-queue-optimal"


if(params.technology == "10XV2"){
    params.soloCBlen = 16;
    params.soloCBstart = 1;
    params.soloUMIstart = 17;
    params.soloUMIlen = 10;
    params.barcodes_whitelist = "737K-august-2016.txt"
}else if(params.technology == "10XV3"){
    params.soloCBlen = 16;
    params.soloCBstart = 1;
    params.soloUMIstart = 17;
    params.soloUMIlen = 12;
    params.barcodes_whitelist = "3M-february-2018.txt"
}else{
    params.soloCBlen = 16;
    params.soloCBstart = 1;
    params.soloUMIstart = 17;
    params.soloUMIlen = 10;
    params.barcodes_whitelist = None;
}



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
         inputformat    : ${params.inputformat} 
         sampleSubfix1    : ${params.sampleSubfix1} 
         sampleSubfix2    : ${params.sampleSubfix2} 
         """
         .stripIndent()

urls  = Channel.from(params.readurls);
ids = Channel.from(params.reads);
ids.merge(urls)
    .view()
    .set { read_pairs_ch }
/*
 * 1. Dump
 */
process Dump {
    cpus 4
    memory '4 GB'
    queue "${params.downloadqueue}"
    errorStrategy 'finish'

    input:
    tuple val(pair_id), val(url) from read_pairs_ch
    output:
    file('*.bam') into results_ch_dump
    val pair_id into id_ch_dump

    script:
    """
    echo "Start to prefetch fastq files from .sra files: ${pair_id}";
    wget -c ${url} -O ${pair_id}".bam"
    ls -la;
    """
}
    
/*
 * 2. Mapping 
 */
process Map {

    publishDir "${params.outdir}", mode: "copy"
    cpus 16
    memory '192 GB'
    queue "${params.mapqueue}"
    errorStrategy 'finish'

    input:
    file result from results_ch_dump
    val pair_id from id_ch_dump
    path ref from params.ref
    output:
    file("${pair_id}/") into results_ch_map
    val pair_id into id_ch_map
    
    script:
    """
    echo "Make output dir ${pair_id}"
    mkdir ${pair_id}
    ls -la
    echo "Filter BAM ${pair_id}"
    samtools view -@ ${params.threads} -e "[${params.barcode}] && [${params.umi}]" -o "${params.baseDir}/${pair_id}_filtered.bam" "${params.baseDir}/${pair_id}.bam"
    perl ${params.codebase}/scvh_map_reads.pl \
    --output-dir "${pair_id}" \
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
    --soloInputSAMattrBarcodeSeq "${params.soloInputSAMattrBarcodeSeq}" \
    "${ref}" \
    "${params.baseDir}/${pair_id}_filtered.bam"        
    """

}
/*
 * 3. Filter
 */
process Filter {
    publishDir "${params.outdir}", mode: "copy",overwrite: true
    errorStrategy { task.attempt > 2 ? 'ignore' : 'retry' }
    maxRetries 2
    queue 'jy-scvh-queue-r5a4x-1'
    cpus 4
    memory '16 GB'
    
    input:
    file result from results_ch_map
    val pair_id from id_ch_map

    output:
    file("${pair_id}/") into results_ch_fil
    val pair_id into id_ch_fil

    shell
    """
    ls
    echo "Filter output dir ${pair_id}"
    Rscript ${params.codebase}/scvh_filter_matrix.r "${pair_id}"
    """
}
/*
 * 4. Analysis
 */
process Analysis {
    publishDir "${params.outdir}", mode: "copy",overwrite: true
    queue 'jy-scvh-queue-r5a4x-1'
    cpus 16
    memory '64 GB'
    errorStrategy { task.attempt > 2 ? 'ignore' : 'retry' }
    maxRetries 2
    input:
    file result from results_ch_fil
    val pair_id from id_ch_fil

    output:
    file("${pair_id}/") into results_ch4
    
    shell
    """
    if [ -e "${pair_id}/filtered_matrix_viral_barcodes_info.txt" -a -e "${pair_id}/filtered_matrix_viral_genes_info.txt" ]; then
        echo "Analysis bam in the output dir ${pair_id}"
        perl ${params.codebase}/scvh_analyze_bam.pl "${pair_id}" -t 10

    else
        echo "Cannot find virus in the sample"
    fi
    """
}