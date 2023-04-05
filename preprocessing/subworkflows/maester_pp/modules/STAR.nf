// STAR module

nextflow.enable.dsl = 2

//

process STAR {

  input:
  path fastq

  output:
  path "boh.bam", emit: bam

  script:
  """
  /STAR-2.7.9a/source/STAR \
    --runThreadN 15 \
    --genomeDir ${params.maester_pp.ref} \
    --readFilesIn ${fastq}
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI nM AS 
  """

  stub:
  """
  touch boh.bam
  """

}