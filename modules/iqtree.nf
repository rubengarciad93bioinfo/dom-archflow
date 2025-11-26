// modules/iqtree.nf
nextflow.enable.dsl=2

process IQTREE_RUN {
  tag { phylip.baseName }

  // Número de hilos que IQ-TREE va a usar
  cpus params.iqtree_cpus ?: 4

  // Muy importante: un árbol cada vez
  maxForks 1

  publishDir "${params.outdir}/trees", mode: 'copy'

  input:
    path phylip

  output:
    path "${phylip.baseName}.*"

  script:
  """
  export PATH=/opt/conda/bin:\$PATH

  iqtree2 \\
    -s ${phylip} \\
    -nm ${params.iqtree_nm ?: 5000} \\
    -bb ${params.iqtree_bb ?: 1000} \\
    -nt ${task.cpus} \\
    -m MFP \\
    -bnni \\
    -pre ${phylip.baseName}
  """
}
