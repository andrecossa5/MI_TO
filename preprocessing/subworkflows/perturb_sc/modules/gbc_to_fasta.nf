// GBC_TO_FASTA module

nextflow.enable.dsl = 2

//

process GBC_TO_FASTA {

  input:
  path GBCs

  output:
  path "GBC_to_align.fa.gz", emit: fasta

  script:
  """
  awk '{ gsub("@",">",\$1); print }' ${GBCs} \
  | tr ' ' '\n' \
  | pigz --fast -p ${task.cpus} \
  > GBC_to_align.fa.gz
  """

  stub:
  """
  touch GBC_to_align.fa.gz
  """
}