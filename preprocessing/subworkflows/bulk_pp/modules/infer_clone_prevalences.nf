// INFER_PREVALENCES module

nextflow.enable.dsl = 2

//

// Infere GBC clones frequencies from their read counts using ground truth spikins data
process INFER_PREVALENCES {

  input:
  path read_count_by_GBC_corrected

  output:
  path 'clonal_prevalences.csv', emit: stats_table
  path 'good_GBCs_bulk.txt', emit: good_GBCs
  path 'spikeins_fit.png', emit: plot

  script:
  """
  python3 \
  ${baseDir}/bin/bulk_pp/infer_clone_prevalences.py \
  ${params.bulk_pp.spikeins_table} \
  ${read_count_by_GBC_corrected}
  """

}
