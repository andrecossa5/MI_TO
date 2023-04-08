// FILTER module

nextflow.enable.dsl = 2

//

process TO_H5AD {

  input:
  path final

  output:
  path "AFM.h5ad", emit: afm

  script:
  """
  python ${baseDir}/bin/to_h5ad.py ${final}
  """

  stub:
  """
  touch AFM.h5ad
  """

}