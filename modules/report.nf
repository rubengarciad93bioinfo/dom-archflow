process REPORT {
  publishDir "${params.outdir}", mode: 'copy'
  input:
    tuple val(org), path(counts)
  output:
    path "summary.txt"
  script:
  """
  echo "Summary dom-archflow" > summary.txt
  printf "%s: " "${org}" >> summary.txt
  cat ${counts} >> summary.txt
  """
}

