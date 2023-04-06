// PREP_GBC module

nextflow.enable.dsl = 2

//

process PREP_GBC {

  input:
  path R1_raw
  path R2_raw
  path reads_aligned

  output:
  path "aligned_R1.fq.gz", emit: R1
  path "aligned_R2.fq.gz", emit: R2

  script:
  """
  seqtk subseq \
  ${R1_raw} \
  ${reads_aligned} \
  | pigz --fast -p ${task.cpus} \
  > aligned_R1.fq.gz

  seqtk subseq \
  ${R2_raw} \
  ${reads_aligned} \
  | pigz --fast -p ${task.cpus} \
  > aligned_R2.fq.gz
  """

  stub:
  """
  touch aligned_R1.fq.gz
  touch aligned_R2.fq.gz
  """

}