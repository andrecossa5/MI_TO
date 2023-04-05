// INDEX module

nextflow.enable.dsl = 2

//

process INDEX {

  input:
  path bam

  output:
  path "${bam}", emit: bam
  path "${bam}.bai", emit: index

  script:
  """
  samtools index -@ 15 ${bam}
  """

  stub:
  """
  touch ${bam}
  touch ${bam}.bai
  """

}