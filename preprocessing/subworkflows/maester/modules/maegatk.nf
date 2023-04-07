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

  maegatk bcall -i ${mitobam} \
  -g rCRS \
  -c ${task.cpus} \
  -b ./barcodes.txt \
  -mr 3 -c 4 -sr
  """

  stub:
  """
  mkdir final
  """

}


// python ${baseDir}/bin/maegatk_cli.py \
// ${params.maester_code_dir} \
// ${mitobam} \
// ${task.cpus} \
// ./barcodes.txt \
// ${params.maester_min_reads} 