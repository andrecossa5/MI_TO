// MERGE module

nextflow.enable.dsl = 2

//

process MERGE {

  input:
  path bam_1
  path bam_2

  output:
  path "merged_mitobam.bam", emit: mitobam

  script:
  """
  samtools index -@ 10 ${bam_1}
  samtools index -@ 10 ${bam_2}
  samtools merge -@ 10 ./merged_mitobam.bam ${bam_1} ${bam_2}
  """

  stub:
  """
  touch merged_mitobam.bam
  """

} 