// EXTRACT_R2 module

nextflow.enable.dsl = 2

//

process EXTRACT_R2 {

  input:
  path R2_raw

  output:
  path "R2_first_33_nt.fa.gz", emit: first_33

  script:
  """
  zcat ${R2_raw} \
  | sed -n '1~4s/^@/>/p;2~4p' \
  | awk '{if(NR%2){print \$0}else{print substr(\$0, 1, 33);}}' \
  | gzip --fast \
  > R2_first_33_nt.fa.gz
  """

}