// measter subworkflow

// Include here
nextflow.enable.dsl = 2
include { MERGE_R1 } from "./modules/merge_R1.nf"
include { MERGE_R2 } from "./modules/merge_R2.nf"
include { ASSEMBLE_FQ } from "./modules/assemble_fastq.nf"
include { STAR } from "./modules/STAR.nf"
include { FILTER_I } from "./modules/filter_bam.nf"
include { FILTER_II } from "./modules/filter_bam.nf"
include { MERGE } from "./modules/merge_bams.nf"
include { INDEX } from "./modules/index_bam.nf"
include { MAEGATK } from "./modules/maegatk.nf"
include { TO_H5AD } from "./modules/to_h5ad.nf"


// 

process publish_maester {

    publishDir "${params.maester_outdir}/${sample}/", mode: 'copy'

    input:
    tuple val(sample), val(in_folder)
    path bam
    path index
    path maegatk_out
    path afm

    output:
    path bam
    path index
    path maegatk_out
    path afm

    script:
    """
    echo moving everything to ${params.maester_outdir}
    """

}

// 

//----------------------------------------------------------------------------//
// maester_pp subworkflow
//----------------------------------------------------------------------------//

workflow maester {
    
    take:
        ch_input
        filtered
        not_enriched_bam  

    main:
        MERGE_R1(ch_input)
        MERGE_R2(ch_input)
        ASSEMBLE_FQ(MERGE_R1.out.R1, MERGE_R2.out.R2)
        STAR(ASSEMBLE_FQ.out.fq)
        FILTER_I(not_enriched_bam)
        FILTER_II(STAR.out.bam)
        MERGE(FILTER_I.out.mitobam, FILTER_II.out.mitobam)
        INDEX(MERGE.out.mitobam)
        MAEGATK(INDEX.out.bam, INDEX.out.index, filtered)
        TO_H5AD(MAEGATK.out.output)
        publish_maester(
            ch_input,
            INDEX.out.bam, 
            INDEX.out.index,
            MAEGATK.out.output, 
            TO_H5AD.out.afm
        )

    emit:
        outputs = MAEGATK.out.output
        afm = TO_H5AD.out.afm

}

//----------------------------------------------------------------------------//