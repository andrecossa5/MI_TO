// Bulk ppp workflow

// Include here
nextflow.enable.dsl = 2
include { SEARCH_PATTERNS } from "./modules/generate_search_patterns.nf"
include { EXTRACT_READS } from "./modules/extract_reads.nf"
include { FIND_GBC } from "./modules/find_GBC.nf"
include { REMOVE_SHORT } from "./modules/remove_GBC_shorter_than_18bp.nf"
include { WHITELIST } from "./modules/generate_GBC_whitelist.nf"
include { FORMAT_WHITELIST } from "./modules/reformat_GBC_whitelist_for_error_correction.nf"
include { CORRECT } from "./modules/correct_GBC.nf"
include { COUNT } from "./modules/count_reads_for_each_GBC.nf"
include { INFER_PREVALENCES } from "./modules/infer_clone_prevalences.nf"

//

// create summary of run
process generate_run_summary_bulk {

  input:
  tuple val(sample), val(in_folder)
  path reads
  path GBC_not_corrected
  path GBC_not_corrected_18bp
  path read_count_by_GBC_corrected
  path good_GBCs

  output:
  path 'log.txt', emit: summary

  script:
  """
  echo "Summary Step 1, sample ${sample}" > log.txt
  echo "-------------------------------------" >> log.txt
  echo "" >> log.txt
  echo "Overview" >> log.txt
  echo "- Date of analysis:  \$(date)" >> log.txt
  echo "- User:              ${USER}" >> log.txt
  echo "- Working directory: ${PWD}" >> log.txt
  echo "" >> log.txt
  echo "Parameters" >> log.txt
  echo "--indir:        ${params.bulk_pp.indir}" >> log.txt
  echo "--outdir:        ${params.bulk_pp.outdir}" >> log.txt
  echo "--anchor_sequence: ${params.bulk_pp.anchor_sequence}" >> log.txt
  echo "${sample} specific I/O: ${params.bulk_pp.indir}/${sample}, ${params.bulk_pp.outdir}/${sample}" >> log.txt
  echo "" >> log.txt
  echo "Numbers" >> log.txt
  echo "- Reads in input:             \$(zcat ${reads} | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> log.txt
  echo "- Reads with GBC:             \$(zcat ${GBC_not_corrected} | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> log.txt
  echo "- Reads with full length GBC: \$(zcat ${GBC_not_corrected_18bp} | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> log.txt
  echo "- Unique GBCs:" >> log.txt
  echo "  - before error-correction:  \$(zcat ${GBC_not_corrected_18bp} | sort | uniq | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> log.txt
  echo "  - after error-correction (i.e., the used as reference for single-cell): \$(tail -n +2 ${read_count_by_GBC_corrected} | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> log.txt
  echo "- Good GBCs (i.e., at least 10 cells, inferred from spikeins interpolation):   \$(cat ${good_GBCs} | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> log.txt
  """

}

// Prep final output
process cleanup_bulk {

  // Publish
  publishDir "${params.bulk_pp.outdir}/${sample}/", mode: 'copy'

  input:
  tuple val(sample), val(in_folder)
  path GBC_corrected 
  path read_count_by_GBC_corrected
  path clonal_prevalences
  path good_GBCs
  path spikeins_fit
  path run_summary

  output:
  path GBC_corrected
  path read_count_by_GBC_corrected
  path clonal_prevalences
  path good_GBCs
  path spikeins_fit
  path run_summary

  script:
  """
  echo "Moving all necessary files to ${params.bulk_pp.outdir}/${sample}/..."
  """

}


//

//----------------------------------------------------------------------------//
// Bulk preprocessing subworkflow
//----------------------------------------------------------------------------//

workflow bulk_pp {

  take:
      ch_input
      // outdir
      // spikeins_table
      // anchor_sequence

  main:
      SEARCH_PATTERNS()
      EXTRACT_READS(ch_input)
      FIND_GBC(SEARCH_PATTERNS.out.search_patterns, EXTRACT_READS.out.reads)
      REMOVE_SHORT(FIND_GBC.out.GBC)
      WHITELIST(REMOVE_SHORT.out.GBC)
      FORMAT_WHITELIST(WHITELIST.out.whitelist)
      CORRECT(REMOVE_SHORT.out.GBC, FORMAT_WHITELIST.out.formatted_whitelist)
      COUNT(CORRECT.out.GBC, WHITELIST.out.whitelist, REMOVE_SHORT.out.GBC)
      INFER_PREVALENCES(COUNT.out.read_counts)

      // Summary and cleanup
      generate_run_summary_bulk(
        ch_input,
        EXTRACT_READS.out.reads, 
        FIND_GBC.out.GBC, 
        REMOVE_SHORT.out.GBC, 
        COUNT.out.read_counts,
        INFER_PREVALENCES.out.good_GBCs
      )
      cleanup_bulk(
        ch_input, 
        CORRECT.out.GBC, 
        COUNT.out.read_counts, 
        INFER_PREVALENCES.out.stats_table,
        INFER_PREVALENCES.out.good_GBCs,
        INFER_PREVALENCES.out.plot,
        generate_run_summary_bulk.out.summary
      )

  emit:
      summary = generate_run_summary_bulk.out.summary

}

//----------------------------------------------------------------------------//