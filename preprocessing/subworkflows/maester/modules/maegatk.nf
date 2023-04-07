// MAEGATK module

nextflow.enable.dsl = 2

//

process MAEGATK {

  input:
  path mitobam
  path index
  path filtered

  output:
  path "final", emit: output
 
  script:
  """
  zcat ./filtered/barcodes.tsv.gz > ./barcodes.txt

  python ${baseDir}/bin/maegatk_cli.py \
  ${params.maester_code_dir} \
  ${mitobam} \
  ${task.cpus} \
  ./barcodes.txt \
  ${params.maester_min_reads} 
  """

  stub:
  """
  mkdir final
  """

}