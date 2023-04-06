// SOLO module

nextflow.enable.dsl = 2

//

process SOLO {

  input:
  path transcript_R1
  path transcript_R2

  output:
  path 'raw', emit: raw 
  path 'filtered', emit: filtered 
  path 'Features.stats', emit: stats 
  path 'Summary.csv', emit: summary 
  path 'Aligned.sortedByCoord.out.bam', emit: bam 

  script:
  """
  STAR \
    --runThreadN ${task.cpus} \
    --genomeDir ${params.ref} \
    --readFilesIn ${transcript_R1} ${transcript_R1} \
    --readFilesCommand zcat \
    --outTmpDir tmp \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI nM AS CR UR CB UB \
    --soloType CB_UMI_Simple \
    --soloBarcodeReadLength 28 \
    --soloCBstart 1 \
    --soloCBlen 16 \
    --soloUMIstart 17 \
    --soloUMIlen 12 \
    --soloCBwhitelist ${params.whitelist} \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --soloUMIdedup 1MM_CR \
    --soloCellFilter EmptyDrops_CR 
  mv Solo.out/Gene/* ./ &&
  rm -rf Solo.out &&
  gzip raw/*.tsv
  gzip filtered/*.tsv
  """

  stub:
  """
  touch raw
  touch filtered
  touch Features.stats 
  touch Summary.csv
  touch Aligned.sortedByCoord.out.bam
  """

}