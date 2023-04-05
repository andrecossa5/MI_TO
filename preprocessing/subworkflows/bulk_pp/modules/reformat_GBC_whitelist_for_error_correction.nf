// FORMAT_WHITELIST module

nextflow.enable.dsl = 2

//

process FORMAT_WHITELIST {

  input:
  path GBC_whitelist 

  output:
  path 'GBC_whitelist_proper_format.tsv', emit: formatted_whitelist

  script:
  """
  Rscript \
  ${baseDir}/bin/bulk_pp/reformat_GBC_whitelist_for_error_correction.R \
  ${GBC_whitelist}
  """

}