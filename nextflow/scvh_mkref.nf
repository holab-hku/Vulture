#!/usr/bin/env nextflow
/*
 * pipeline input parameters
 */
params.baseDir = ".";
params.codebase = "/code";
params.threads = 16;
params.virus_database = "viruSITE.NCBIprokaryotes";
params.output = "newref"

if(params.virus_database == "NCBIprokaryotes"){
    params.genomeprefix = "human_host_microbes"
}else if(params.virus_database == "viruSITE"){
    params.genomeprefix = "human_host_viruses"
}else if(params.virus_database == "viruSITE.NCBIprokaryotes"){
    params.genomeprefix = "human_host_viruses_microbes"
}

log.info """\
         S C V H M K R E F- N F   P I P E L I N E
         ===================================
         transcriptome:     ${params.ref}
         human_fa_file:     ${params.humanfa}
         human_gtf_file:    ${params.humagtf}
         viruSITE_file:     ${params.viruSITE}
         prokaryotes_file:  ${params.prokaryotes}
         virus_database:    ${params.virus_database}
         genomeprefix:      ${params.genomeprefix}
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
    file("${params.output}") into results_ch1
    
    script:
        """
        mkdir ${params.output};
        ls -la
        perl ${params.codebase}/virusl_et.pl -o ${params.output} --human_fa ${ref}/${params.humanfa} --human_gtf ${ref}/${params.humagtf} --viruSITE ${ref}/${params.viruSITE} --prokaryotes ${ref}/${params.prokaryotes} --database ${params.virus_database} --output_prefix ${params.genomeprefix};
        ls -la ${params.output}/human_host_viruses_reference_set/with_hg38
        mv ${params.output}/human_host_viruses_reference_set/with_hg38/human_host_*.with_hg38.fa ${params.output}/
        mv ${params.output}/human_host_viruses_reference_set/with_hg38/human_host_*.with_hg38.removed_amb_viral_exon.gtf ${params.output}/

        wget https://raw.githubusercontent.com/10XGenomics/cellranger/master/lib/python/cellranger/barcodes/translation/3M-february-2018.txt.gz
        gunzip -c 3M-february-2018.txt.gz > ${params.output}/3M-february-2018.txt
        wget https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-august-2016.txt -P ${params.output}
        """
}
/*
 * 2. mkref
 */
 process Mkref {
    container 'junyichen6/vulture:0.0.1'
    publishDir "${params.outdir}", mode: "copy",overwrite: false
    cpus 16
    memory '64 GB'

    input:
    file ref from results_ch1
    output:
    file("${params.output}") into results_ch2
    
    script:
        """
        ls -la
        chmod -R 777 ${ref}
        ls -la ${ref}
        STAR --runThreadN ${params.threads} --runMode genomeGenerate --genomeDir ${ref} --genomeFastaFiles ${ref}/${params.genomeprefix}.${params.virus_database}.with_hg38.fa
        """
}