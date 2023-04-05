// PREP_GBC module

nextflow.enable.dsl = 2

//

process PREP_TRANSCRIPT {

  input:
  file R1_raw
  file R2_raw
  file reads_transcriptomic

  output:
  path "transcriptomic_R1.fq.gz", emit: R1
  path "transcriptomic_R2.fq.gz", emit: R2

  script:
  """
  seqtk subseq \
  ${R1_raw} \
  ${reads_transcriptomic} \
  | gzip \
  > transcriptomic_R1.fq.gz

  seqtk subseq \
  ${R2_raw} \
  ${reads_transcriptomic} \
  | gzip --fast \
  > transcriptomic_R2.fq.gz
  """

  stub:
  """
  touch transcriptomic_R1.fq.gz
  touch transcriptomic_R2.fq.gz
  """

}