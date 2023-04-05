// FASTA_FROM_REF module

nextflow.enable.dsl = 2

//

process FASTA_FROM_REF {
  
  input:
  tuple val(sample), val(in_folder)

  output:
  path "GBC_reference.fa", emit: fasta

  script:
  """
  awk 'FNR > 1 {print ">"NR-1"\\n"\$1}' \
  ${params.sc_pp.step_1_out}/${sample}/read_count_by_GBC_corrected.tsv > GBC_reference.fa
  """

}