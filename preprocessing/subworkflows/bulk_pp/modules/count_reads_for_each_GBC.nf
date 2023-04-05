// COUNT module

nextflow.enable.dsl = 2

//

// Counts unique (corrected) GBCs reads 
process COUNT {

  input:
  path GBC_corrected 
  path GBC_whitelist 
  path GBC_not_corrected_18bp 

  output:
  path 'read_count_by_GBC_corrected.tsv', emit: read_counts

  script:
  """
  Rscript \
  ${baseDir}/bin/bulk_pp/count_reads_for_each_GBC.R \
  ${GBC_corrected} \
  ${GBC_whitelist} \
  ${GBC_not_corrected_18bp}
  """

}

