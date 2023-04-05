// ALIGN_33_R2 module

nextflow.enable.dsl = 2

//

process ALIGN_33_R2 {

  input:
  path R2_first_33nt
  path bowtie_index_GBC_pattern_folder

  output:
  path "R2_first_33_nt.fq3.gz", emit: R2_aligned

  script:
  """
  bowtie2 \
  -N 1 \
  --norc \
  -x ${bowtie_index_GBC_pattern_folder}/GBC_pattern \
  -f ${R2_first_33nt} \
  | gzip --fast \
  > R2_first_33_nt.fq3.gz
  """

  stub:
  """
  touch R2_first_33_nt.fq3.gz
  """

}