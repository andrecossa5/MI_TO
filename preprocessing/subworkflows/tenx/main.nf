// sc_pp workflow

// Include here
nextflow.enable.dsl = 2
include { MERGE_R1 } from "../maester/modules/merge_R1.nf"
include { MERGE_R2 } from "../maester/modules/merge_R2.nf"
include { SOLO } from "../perturb_sc/modules/Solo.nf"

// 

process publish_tenx {

    publishDir "${params.tenx_outdir}/${sample}/", mode: 'copy'

    input:
    tuple val(sample), val(in_folder)
    path raw
    path filtered
    path stats
    path summary
    path bam

    output:
    path raw
    path filtered
    path stats
    path summary
    path bam

    script:
    """
    echo moving everything to ${params.tenx_outdir}
    """

}

// 


//----------------------------------------------------------------------------//
// Solo subworkflow
//----------------------------------------------------------------------------//

workflow tenx {

    take:
        ch_input

    main:
        MERGE_R1(ch_input)
        MERGE_R2(ch_input)
        SOLO(MERGE_R1.out.R1, MERGE_R2.out.R2)
        publish_tenx(
            ch_input,
            SOLO.out.raw,
            SOLO.out.filtered,
            SOLO.out.stats,
            SOLO.out.summary,
            SOLO.out.bam
        )

    emit:
        filtered = SOLO.out.filtered
        bam = SOLO.out.bam

}

//----------------------------------------------------------------------------//