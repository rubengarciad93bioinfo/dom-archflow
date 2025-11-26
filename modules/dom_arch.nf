// modules/dom_arch.nf
nextflow.enable.dsl=2

process DOM_ARCH {
  tag { org }
  cpus 1
  publishDir "${params.outdir}/architectures", mode: 'copy'

  input:
    // org, pfam_out, tax_id, taxonomy, completeness
    tuple val(org), path(pfam_out), val(tax_id), val(taxonomy), val(completeness)

  output:
    tuple val(org), path("${org}_architectures.csv")

  script:
  """
  python3 ${projectDir}/scripts/build_architectures.py \\
    --pfam ${pfam_out} \\
    --proteome-id ${org} \\
    --out ${org}_architectures.csv \\
    --version ${params.dataset_version} \\
    --overlap ${params.overlap_threshold} \\
    --tax-id ${tax_id} \\
    --taxonomy "${taxonomy}" \\
    --completeness ${completeness}
  """
}
