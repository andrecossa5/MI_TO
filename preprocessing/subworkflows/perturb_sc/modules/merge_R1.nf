// MERGE_R1 module

nextflow.enable.dsl = 2

//

process MERGE_R1 {

  input:
  tuple val(sample), val(in_folder)

  output:
  path "R1_raw.fq.gz", emit: R1

  script:
  """
  zcat ${in_folder}/*R1*.fastq.gz \
  | awk '{if(NR%4==1){print "@"(NR%1?c+1:++c)} else {print \$0}}' \
  | pigz --fast -p ${task.cpus}  \
  > R1_raw.fastq.gz
  """

  stub:
  """
  touch R1_raw.fastq.gz
  """

}