// FIX_TAGS module

nextflow.enable.dsl = 2

//

process FIX_TAGS {

  input:
  path mitobam_no_UB_CB

  output:
  path "mitobam_fixed.bam", emit: mitobam

  script:
  """
  python ${baseDir}/bin/fix_tags.py ${mitobam_no_UB_CB}
  """

  stub:
  """
  touch mitobam_fixed.bam
  """

}