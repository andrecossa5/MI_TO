// nf-perturbseq-pipeline
 
nextflow.enable.dsl = 2
include { sc_pp } from "./subworkflows/sc_pp/main"
include { maester_pp } from "./subworkflows/maester_pp/main"
include { tenx_pp } from "./subworkflows/tenx_pp/main"

//

// Perturb-seq input sc_fastqs 
ch_perturb = Channel
    .fromPath("${params.sc_pp.indir}/*", type:'dir') 
    .map{ tuple(it.getName(), it) }

//

// MAESTER input sc_fastq
ch_MAESTER = Channel
    .fromPath("${params.maester_pp.indir}/*", type:'dir') 
    .map{ tuple(it.getName(), it) }

// 

// 10x only input sc_fastqs
ch_tenx = Channel
    .fromPath("${params.tenx_pp.indir}/*", type:'dir') 
    .map{ tuple(it.getName(), it) }

//

//----------------------------------------------------------------------------//
// Perturb-seq pipeline entry points and (mock) main workflow
//----------------------------------------------------------------------------//

//

workflow tenx_only {

    tenx_pp(ch_tenx)
    tenx_pp.out.filtered.view()
    tenx_pp.out.bam.view()

}

//

workflow perturbseq_only {

    sc_pp(ch_perturb)
    sc_pp.out.filtered.view()
    sc_pp.out.bam.view()

}

//

workflow tenx_mito {

    tenx_pp(ch_tenx)
    maester_pp(ch_MAESTER, tenx_pp.out.filtered, tenx_pp.out.bam)
    maester_pp.out.outputs.view()
    maester_pp.out.afm.view()

}

//

workflow gbc_mito {

    sc_pp(ch_perturb)
    maester_pp(ch_MAESTER, sc_pp.out.filtered, sc_pp.out.bam)
    maester_pp.out.outputs.view()
    maester_pp.out.afm.view()

}

//

// Mock
workflow  {

    Channel.of(1,2,3,4) | view

}

//----------------------------------------------------------------------------//
