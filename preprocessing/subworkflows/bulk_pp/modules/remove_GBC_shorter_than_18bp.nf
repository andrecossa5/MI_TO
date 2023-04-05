// REMOVE_SHORT module

nextflow.enable.dsl = 2

//

// Process
process REMOVE_SHORT {

  input:
  path GBC_not_corrected

  output:
  path 'GBC_not_corrected_18bp.tsv.gz', emit: GBC

  script:
  """
  Rscript \
  ${baseDir}/bin/bulk_pp/remove_GBC_shorter_than_18bp.R \
  ${GBC_not_corrected}
  """

}