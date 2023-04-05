// GET_GBC_ELEMENTs module

nextflow.enable.dsl = 2

//

process GET_GBC_ELEMENTS {
  
  input:
  path R1
  path R2
  path filtered

  output:
  path 'filtered_R1.fq.gz', emit: R1_filtered
  path 'filtered_R2.fq.gz', emit: R2_filtered 
  path 'CBCs_by_read.tsv', emit: CBCs
  path 'UMIs_by_read.tsv', emit: UMIs
  path 'GBCs_by_read.tsv', emit: GBCs

  script:
  """
  python ${baseDir}/bin/sc_pp/filter_and_extact_from_GBC_reads.py ${R1} ${R2} ${filtered}
  """

  stub:
  """
  touch filtered_R1.fq.gz
  touch filtered_R2.fq.gz
  touch CBCs_by_read.tsv
  touch UMIs_by_read.tsv
  touch GBCs_by_read.tsv
  """

}

