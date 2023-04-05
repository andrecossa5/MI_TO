// CORRECT module

nextflow.enable.dsl = 2

//

process CORRECT {

  input:
  path GBC_not_corrected_18bp 
  path reformatted_GBC_whitelist 

  output:
  path 'GBC_corrected.tsv.gz', emit: GBC

  script:
  """
  perl \
  ${baseDir}/bin/bulk_pp/correct_barcodes.pl \
  ${reformatted_GBC_whitelist} \
  ${GBC_not_corrected_18bp} \
  | gzip --fast \
  > GBC_corrected.tsv.gz
  """

}