// BOWTIE_INDEX_REF module

nextflow.enable.dsl = 2

//

process BOWTIE_INDEX_REF {

  input:
  path GBC_reference_fa

  output:
  path "custom_ref_genome", emit: index

  script:
  """
  mkdir -p custom_ref_genome

  bowtie2-build \
  -f ${GBC_reference_fa} \
  custom_ref_genome/GBC_reference
  """

  stub:
  """
  touch custom_ref_genome
  """

}