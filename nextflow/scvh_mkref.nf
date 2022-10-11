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
    params.genomeprefix = "microbes.NCBIprokaryotes"
}else if(params.virus_database == "viruSITE"){
    params.genomeprefix = "viruses.viruSITE"
}else if(params.virus_database == "viruSITE.NCBIprokaryotes"){
    params.genomeprefix = "viruses_microbes.viruSITE.NCBIprokaryotes"
}

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
    file("${params.output}") into results_ch1
    
    script:
        """
        mkdir ${params.output};
        ls -la
        perl ${params.codebase}/virusl_et.pl -o ${params.output} --human_fa ${ref}/${params.humanfa} --human_gtf ${ref}/${params.humagtf} --viruSITE ${ref}/${params.viruSITE} --prokaryotes ${ref}/${params.prokaryotes} --database ${params.virus_database};
        
        mv ${params.output}/human_host_viruses_reference_set/with_hg38/human_host_${params.genomeprefix}.with_hg38.fa ${params.output}/
        mv ${params.output}/human_host_viruses_reference_set/with_hg38/human_host_${params.genomeprefix}.with_hg38.removed_amb_viral_exon.gtf ${params.output}/

        wget https://raw.githubusercontent.com/10XGenomics/cellranger/master/lib/python/cellranger/barcodes/translation/3M-february-2018.txt.gz
        gunzip -c 3M-february-2018.txt.gz > ${params.output}/3M-february-2018.txt
        wget https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-august-2016.txt -P ${params.output}
        """
}
/*
 * 2. mkref
 */
 process Mkref {
    container '795465341760.dkr.ecr.us-east-1.amazonaws.com/dev/scvh:latest'
    publishDir "${params.outdir}", mode: "copy",overwrite: true
    cpus 10
    memory '64 GB'

    input:
    path ref from results_ch1
    output:
    file("${params.output}") into results_ch2
    
    script:
        """
        "STAR --runThreadN 10 --runMode genomeGenerate --genomeDir ${ref} --genomeFastaFiles ${ref}/human_host_${params.genomeprefix}.with_hg38.fa"
        """
}