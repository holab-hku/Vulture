#!/usr/bin/env nextflow
/*
 * pipeline input parameters
 */
params.baseDir = "."
params.codebase = "/code"
params.threads = 16;
params.virus_database = "viruSITE";



log.info """\
         S C V H M K R E F- N F   P I P E L I N E
         ===================================
         transcriptome: ${params.ref}
         """
         .stripIndent()

/*
 * 1. Download and 
 */
process Download {

    publishDir "${params.outdir}", mode: "copy"
    cpus 16
    memory '64 GB'

    input:
    path ref from params.ref

    output:
    set file("/human_host_viruses_reference_set") into results_ch
    
    script:
        """
        ls
        perl ${params.codebase}/virusl_et.pl \
        -o "/" \
        --human_fa "${ref}/${params.humanfa}" \
        --human_gtf "${ref}/${params.humagtf}" \
        --viruSITE "${ref}/${params.viruSITE}" \
        --prokaryotes "${ref}/${params.prokaryotes}" ;
        ls -la
        """
}
/*
 * 2. Filter
 */