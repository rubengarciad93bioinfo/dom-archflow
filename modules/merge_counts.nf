process MERGE_COUNTS {
  publishDir "${params.outdir}", mode: 'copy'
  input:
    path counts_files 
  output:
    path "summary.txt"
  script:
  """
  echo "Summary dom-archflow" > summary.txt
  cat ${counts_files} >> summary.txt
  """
}
