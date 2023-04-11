// STAR module

nextflow.enable.dsl = 2

//

process STAR {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(fastq)

  output:
  tuple val(sample_name), path("Aligned.sortedByCoord.out.bam"), emit: bam

  script:
  """
  /STAR-2.7.9a/source/STAR \
    --runThreadN ${task.cpus} \
    --genomeDir ${params.ref} \
    --readFilesIn ${fastq} \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI nM AS
  """

  stub:
  """
  touch Aligned.sortedByCoord.out.bam
  """

}