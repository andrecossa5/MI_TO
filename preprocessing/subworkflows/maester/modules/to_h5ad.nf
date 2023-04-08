// TO_H5AD module

nextflow.enable.dsl = 2

//

process TO_H5AD {

  input:
  path output

  output:
  path "AFM.h5ad", emit: afm

  script:
  """
  python ${baseDir}/bin/to_h5ad.py ${output}
  """
  
  stub:
  """
  touch AFM.h5ad
  """

}