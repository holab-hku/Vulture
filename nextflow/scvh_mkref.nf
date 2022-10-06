#!/usr/bin/env nextflow
/*
 * pipeline input parameters
 */
params.baseDir = ".";
params.codebase = "/code";
params.threads = 16;
params.virus_database = "viruSITE";
params.output = "result";

log.info """\
         S C V H M K R E F- N F   P I P E L I N E
         ===================================
         transcriptome: ${params.ref}
         human_fa_file:  ${params.humanfa}
         human_gtf_file: ${params.humagtf}
         viruSITE_file:  ${params.viruSITE}
         prokaryotes_file: ${params.prokaryotes}
         """
         .stripIndent()

/*
 * 1. Download and 
 */
process Downloadref {

    publishDir "${params.outdir}", mode: "copy",overwrite: true
    cpus 4
    memory '16 GB'

    input:
    path ref from params.ref
    output:
    file("${params.output}") into results_ch
    
    script:
        """
        ls -la
        mkdir ${params.output};
        ls -la
        perl ${params.codebase}/virusl_et.pl \
        -o "${params.output}" \
        --human_fa "${ref}/${params.humanfa}" \
        --human_gtf "${ref}/${params.humagtf}" \
        --viruSITE "${ref}/${params.viruSITE}" \
        --prokaryotes "${ref}/${params.prokaryotes}" ;
        
        echo "Moving files"
        mv "${params.output}/human_host_viruses_reference_set/with_hg38/human_host_viruses_microbes.viruSITE.NCBIprokaryotes.with_hg38.fa" "${params.output}/human_host_viruses_microbes.viruSITE.NCBIprokaryotes.with_hg38.fa"
        mv "${params.output}/human_host_viruses_reference_set/with_hg38/human_host_viruses_microbes.viruSITE.NCBIprokaryotes.with_hg38.removed_amb_viral_exon.gtf" "${params.output}/human_host_viruses_microbes.viruSITE.NCBIprokaryotes.with_hg38.removed_amb_viral_exon.gtf"

        echo "Downloading the 10x barcode whitelists"
        wget https://raw.githubusercontent.com/10XGenomics/cellranger/master/lib/python/cellranger/barcodes/translation/3M-february-2018.txt.gz
        gunzip -c 3M-february-2018.txt.gz > "${params.output}/3M-february-2018.txt"
        wget https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-august-2016.txt -P "${params.output}"
        ls "${params.output}"
        echo "Finish"
        """
}
/*
 * 2. Filter
 */