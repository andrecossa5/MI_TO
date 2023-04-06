// FILTER module

nextflow.enable.dsl = 2

//

process FILTER_I {

  input:
  path bam

  output:
  path "mitobam_I.bam", emit: mitobam

  script:
  """
  samtools index -@ ${task.cpus} ${bam}
  samtools view ${bam} -b -@ ${task.cpus} chrM > mitobam_I.bam
  """

  stub:
  """
  touch mitobam_I.bam
  """

}

//

process FILTER_II {

  input:
  path bam

  output:
  path "mitobam_II.bam", emit: mitobam

  script:
  """
  samtools index -@ 8 ${bam}
  samtools view ${bam} -b -@ ${task.cpus} chrM > mitobam_II.bam
  """

  stub:
  """
  touch mitobam_II.bam
  """

}