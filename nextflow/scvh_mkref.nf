#!/usr/bin/env nextflow
/*
 * pipeline input parameters
 */
params.baseDir = ".";
params.codebase = "/code";
params.threads = 16;
params.virus_database = "viruSITE";
params.output = "human_host_viruses_reference_set";



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
process Download {

    publishDir "./${params.outdir}", mode: "copy",overwrite: true
    cpus 16
    memory '64 GB'

    input:
    path ref from params.ref

    output:
    file("${params.output}") into results_ch
    
    script:
        """
        ls -la
        perl ${params.codebase}/virusl_et.pl \
        -o "./${params.output}" \
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