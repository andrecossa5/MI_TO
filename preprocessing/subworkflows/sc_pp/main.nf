// Sc_pp subworkflow

nextflow.enable.dsl = 2

// Include here
include { MERGE_R1 } from "./modules/merge_R1.nf"
include { MERGE_R2 } from "./modules/merge_R2.nf"
include { EXTRACT_R2 } from "./modules/extract_first_33_R2.nf"
include { BOWTIE_INDEX_GBC_PATTERN } from "./modules/create_bowtie_index_pattern.nf"
include { BOWTIE_INDEX_REF } from "./modules/create_bowtie_index_ref.nf"
include { ALIGN_33_R2 } from "./modules/align_first_33_R2.nf"
include { GET_READS_NAMES } from "./modules/get_names_of_all_reads.nf"
include { GET_NAMES_ALIGNED } from "./modules/get_names_aligned.nf"
include { GET_NAMES_NOT_ALIGNED } from "./modules/get_names_not_aligned.nf"
include { PREP_GBC } from "./modules/prep_GBC.nf"
include { PREP_TRANSCRIPT } from "./modules/prep_transcript.nf"
include { SOLO } from "./modules/Solo.nf"
include { FASTA_FROM_REF } from "./modules/fasta_from_ref.nf"
include { GET_GBC_ELEMENTS } from "./modules/filter_and_extract_from_GBC.nf"
include { GBC_TO_FASTA } from "./modules/gbc_to_fasta.nf"
include { ALIGN_GBC } from "./modules/align_GBC_to_ref.nf"
include { CELL_ASSIGNMENT } from "./modules/cell_assignment.nf"

//

process generate_run_summary_sc {
 
    input:
    tuple val(sample), val(in_folder)
    path all_reads
    path reads_transcript
    path reads_aligned
    path GBCs
    path filtered
  
    output:
    path "run_summary.txt", emit: run_summary

    script:
    """
    echo "Summary Step 2, sample ${sample}" > run_summary.txt
    echo "-------------------------------------" >> run_summary.txt
    echo "" >> run_summary.txt
    echo "Overview" >> run_summary.txt
    echo "- Date of analysis:  \$(date)" >> run_summary.txt
    echo "- User:              ${USER}" >> run_summary.txt
    echo "- Working directory: ${PWD}" >> run_summary.txt
    echo "" >> run_summary.txt
    echo "Parameters" >> run_summary.txt
    echo "--indir:                ${params.sc_pp.indir}" >> run_summary.txt
    echo "--outdir:               ${params.sc_pp.outdir}" >> run_summary.txt
    echo "${sample} specific I/O: ${params.sc_pp.indir}/${sample}, ${params.sc_pp.outdir}/${sample}" >> run_summary.txt
    echo "--step_1_out:           ${params.sc_pp.step_1_out}" >> run_summary.txt
    echo "--pattern:              ${params.sc_pp.pattern}" >> run_summary.txt
    echo "--ref:                  ${params.sc_pp.ref}" >> run_summary.txt
    echo "Numbers" >> run_summary.txt
    echo "- Reads in input:                  \$(cat ${all_reads} | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    echo "- Transcriptomic reads:            \$(cat ${reads_transcript} | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    echo "- GBC-containing reads:            \$(cat ${reads_aligned} | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    echo "- Unique GBC in reference:         \$(cat ${params.sc_pp.step_1_out}/${sample}/read_count_by_GBC_corrected.tsv | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    echo "- Unique GBC found in this sample: \$(cat ${GBCs} | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    echo "- Putative cell n (Solo cell-calling): \$(zcat ${filtered}/barcodes.tsv.gz | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    echo "- Total number of transcripts:     \$(zcat ${filtered}/features.tsv.gz | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    """

}

//

process publish_sc_pp {

    publishDir "${params.sc_pp.outdir}/${sample}/", mode: 'copy'

    input:
    tuple val(sample), val(in_folder)
    path CBC_GBC
    path CBC_GBC_plot
    path cells_summary
    path clones_summary
    path bam
    path stats
    path summary
    path filtered
    path raw
    path run_summary

    output:
    path raw
    path CBC_GBC
    path CBC_GBC_plot
    path cells_summary
    path clones_summary
    path bam
    path stats
    path summary
    path filtered
    path run_summary

    script:
    """
    echo "Moving all output files to ${params.sc_pp.outdir}/${sample}/..."
    """

}

//----------------------------------------------------------------------------//
// sc_pp subworkflow
//----------------------------------------------------------------------------//

workflow sc_pp {
    
    take:
        ch_input   

    main:

        // Prep reads
        MERGE_R1(ch_input)
        MERGE_R2(ch_input)
        EXTRACT_R2(MERGE_R2.out.R2)
        BOWTIE_INDEX_GBC_PATTERN()
        ALIGN_33_R2(EXTRACT_R2.out.first_33, BOWTIE_INDEX_GBC_PATTERN.out.index)
        GET_READS_NAMES(MERGE_R1.out.R1)
        GET_NAMES_ALIGNED(ALIGN_33_R2.out.R2_aligned)
        GET_NAMES_NOT_ALIGNED(GET_READS_NAMES.out.names, GET_NAMES_ALIGNED.out.names)
        PREP_GBC(MERGE_R1.out.R1, MERGE_R2.out.R2, GET_NAMES_ALIGNED.out.names) // Control
        PREP_TRANSCRIPT(MERGE_R1.out.R1, MERGE_R2.out.R2, GET_NAMES_NOT_ALIGNED.out.names) // Control
        
        // STARSolo
        SOLO(PREP_TRANSCRIPT.out.R1, PREP_TRANSCRIPT.out.R2)

        // Assign cells to clones
        GET_GBC_ELEMENTS(PREP_GBC.out.R1, PREP_GBC.out.R2, SOLO.out.filtered)
        FASTA_FROM_REF(ch_input)
        BOWTIE_INDEX_REF(FASTA_FROM_REF.out.fasta)
        GBC_TO_FASTA(GET_GBC_ELEMENTS.out.GBCs)
        ALIGN_GBC(BOWTIE_INDEX_REF.out.index, GBC_TO_FASTA.out.fasta)
        CELL_ASSIGNMENT(GET_GBC_ELEMENTS.out.CBCs, GET_GBC_ELEMENTS.out.UMIs, GET_GBC_ELEMENTS.out.GBCs, ALIGN_GBC.out.names)

        // Summary and cleanup
        generate_run_summary_sc(
            ch_input,
            GET_READS_NAMES.out.names,
            GET_NAMES_NOT_ALIGNED.out.names,
            GET_NAMES_ALIGNED.out.names,
            GET_GBC_ELEMENTS.out.GBCs,
            SOLO.out.filtered
        )
        publish_sc_pp(
            ch_input,
            CELL_ASSIGNMENT.out.CBC_GBC_combos,
            CELL_ASSIGNMENT.out.plot,
            CELL_ASSIGNMENT.out.cells_summary,
            CELL_ASSIGNMENT.out.clones_summary,
            SOLO.out.bam,
            SOLO.out.stats,
            SOLO.out.summary,
            SOLO.out.filtered,
            SOLO.out.raw,
            generate_run_summary_sc.out.run_summary
        )

    emit:
        summary = generate_run_summary_sc.out.run_summary
        filtered = SOLO.out.filtered
        bam = SOLO.out.bam

}

//----------------------------------------------------------------------------// 