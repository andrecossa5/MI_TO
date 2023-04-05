// EXTRACT_READS module

nextflow.enable.dsl = 2

//

process EXTRACT_READS {

  input:
  tuple val(sample), val(in_folder)

  output:
  path 'reads.tsv.gz', emit: reads

  script:
  """
  zcat ${in_folder}/*.fastq.gz | awk 'NR % 4 == 2' | gzip --fast > reads.tsv.gz
  """
  
}