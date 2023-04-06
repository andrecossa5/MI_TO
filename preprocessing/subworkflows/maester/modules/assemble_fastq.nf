// ALIGN_33_R2 module

nextflow.enable.dsl = 2

//

process ASSEMBLE_FQ {

  input:
  path R1
  path R2

  output:
  path "assembled.fastq.gz", emit: fq

  script:
  """
  python ${baseDir}/bin/assemble_trim_fastq.py ${R1} ${R2}
  """

  stub:
  """
  touch assembled_fastq.gz
  """

}