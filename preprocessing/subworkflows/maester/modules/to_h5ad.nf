// FILTER module

nextflow.enable.dsl = 2

//

process TO_H5AD {

  input:
  path output

  output:
  path "AFM.h5ad", emit: afm

  script:
  """
  Rscript ${baseDir}/bin/to_h5ad.r ${output}/final AFM
  """

  stub:
  """
  touch AFM.h5ad
  """

}