// FIND_GBC module

nextflow.enable.dsl = 2
params.bin_dir = "../bin/"

//

// Process
process FIND_GBC {

  input:
  path search_patterns
  path reads

  output:
  path 'GBC_not_corrected.tsv.gz', emit: GBC

  script:
  """
  zcat ${reads} \
  | egrep -f ${search_patterns} -o \
  | awk '{print substr(\$0, 23, 18);}' \
  | gzip --fast \
  > GBC_not_corrected.tsv.gz
  """

}