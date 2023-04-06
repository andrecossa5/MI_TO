// GET_NAMES_NOT_ALIGNED module

nextflow.enable.dsl = 2

//

process GET_NAMES_NOT_ALIGNED {

  input:
  path all_reads
  path reads_aligned

  output:
  path "reads_transcriptomic.txt", emit: names

  script:
  """
  LC_ALL=C comm -23 \
  ${all_reads} \
  ${reads_aligned} \
  > reads_transcriptomic.txt
  """

  stub:
  """
  touch reads_transcriptomic.txt
  """

}