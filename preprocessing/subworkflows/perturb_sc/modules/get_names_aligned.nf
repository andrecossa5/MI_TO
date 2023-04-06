// GET_NAMES_ALIGNED module

nextflow.enable.dsl = 2

//

process GET_NAMES_ALIGNED {

  input:
  path R2_first_33nt_aligned

  output:
  path "reads_aligned.txt", emit: names

  script:
  """
  zcat ${R2_first_33nt_aligned} \
  | grep -P '\tconstruct' \
  | cut -f1 \
  > reads_aligned_unsorted.txt

  LC_ALL=C sort --parallel=4 reads_aligned_unsorted.txt > reads_aligned.txt

  rm reads_aligned_unsorted.txt
  """

  stub:
  """
  touch reads_aligned.txt
  """

}
