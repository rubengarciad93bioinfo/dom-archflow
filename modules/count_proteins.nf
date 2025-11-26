process COUNT {
  tag { org }
  input:
    tuple val(org), path(faa)
  output:
    tuple val(org), path("counts_${org}.txt")
  script:
  """
  awk -v ORG="${org}" 'BEGIN{c=0} /^>/{c++} END{printf "%s: proteins=%d\\n", ORG, c}' \
    "${faa}" > "counts_${org}.txt"
  """
}
