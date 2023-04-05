// MERGE_R2 module

nextflow.enable.dsl = 2

//

process MERGE_R2 {

  input:
  tuple val(sample), val(in_folder)

  output:
  path "R2_raw.fq.gz", emit: R2

  script:
  """
  # cat ${in_folder}/*R2*.fastq.gz > R2_raw.fq.gz

  zcat ${in_folder}/*R2*.fastq.gz \
  | awk '{if(NR%4==1){print "@"(NR%1?c+1:++c)} else {print \$0}}' \
  | gzip --fast \
  > R2_raw.fq.gz
  """

  stub:
  """
  touch R2_raw.fq.gz
  """
}