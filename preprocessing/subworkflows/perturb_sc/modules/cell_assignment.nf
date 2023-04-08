// CELL_ASSIGNMENT module

nextflow.enable.dsl = 2

//

process CELL_ASSIGNMENT {

  input:
  // tuple val(sample), val(in_folder)
  path CBCs
  path UMIs 
  path GBCs
  file reads_aligned_to_ref

  output:
  path "CBC_GBC_combos.tsv.gz", emit: CBC_GBC_combos
  path "clones_summary_table.csv", emit: clones_summary
  path "cells_summary_table.csv", emit: cells_summary
  path "CBC_GBC_combo_status.png", emit: plot

  script:
  """
  python ${baseDir}/bin/cell_assignment.py \
  ${CBCs} ${UMIs} ${GBCs} ${reads_aligned_to_ref}
  """

  stub:
  """
  touch CBC_GBC_combos.tsv.gz
  touch clones_summary_table.csv
  touch cells_summary_table.csv
  touch CBC_GBC_combo_status.png
  """

}
