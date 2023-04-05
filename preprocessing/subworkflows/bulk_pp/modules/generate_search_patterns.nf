// SEARCH_PATTERNS module

nextflow.enable.dsl = 2

//

process SEARCH_PATTERNS {

  output:
  path "search_patterns.tsv", emit: search_patterns

  script:
  """
  Rscript \
  ${baseDir}/bin/bulk_pp/generate_search_patterns.R \
  ${params.bulk_pp.anchor_sequence}
  """

}
