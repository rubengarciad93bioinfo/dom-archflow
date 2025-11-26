// modules/phylo_matrices.nf
nextflow.enable.dsl=2

process PHYLO_MATRICES {
  tag "phylo_matrices"
  cpus 1
  publishDir "${params.outdir}/phylo", mode: 'copy'

  input:
    // un único CSV: all_architectures o tu subset
    path all_arch

  output:
    // solo nos interesa sacar los .phylip como canal
    path "*.phylip"

  script:
  """
  export PATH=/opt/conda/bin:\$PATH

  python \\
    ${projectDir}/scripts/build_phylo_matrices.py \\
    --all-architectures ${all_arch}
  """
}
