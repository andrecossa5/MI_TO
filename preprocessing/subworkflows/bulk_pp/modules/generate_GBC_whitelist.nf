// WHITELIST module

nextflow.enable.dsl = 2

//

// Process
process WHITELIST {

  input:
  path GBC_not_corrected_18bp

  output:
  path 'GBC_whitelist.tsv', emit: whitelist

  script:
  """
  python3 \
  ${baseDir}/bin/bulk_pp/generate_GBC_whitelist.py \
  --input ${GBC_not_corrected_18bp} \
  --method directional \
  --threshold 1 \
  --output GBC_whitelist.tsv
  """

}

