// MERGE_R1 module

nextflow.enable.dsl = 2

//

process MERGE_R1 {

  input:
  tuple val(sample), val(in_folder)

  output:
  path "R1_raw.fastq.gz", emit: R1

  script:
  """
  zcat ${in_folder}/*R1*.fastq.gz | pigz --fast -p ${task.cpus} > R1_raw.fastq.gz
  """

  stub:
  """
  touch R1_raw.fastq.gz
  """

}