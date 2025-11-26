// modules/pfam_scan.nf

process PFAM_SCAN {
  tag { org }

  /*
   * cpus = nº de cores que le pasamos a hpc_pfam_search (-cpu)
   * maxForks = cuántos PFAM_SCAN (organismos) a la vez
   */
  cpus params.pfam_cpus
  maxForks params.pfam_maxForks

  // etiquetamos la carpeta de salida con el perfil de pfam
  publishDir "${params.outdir}/pfam_${params.pfam_label}", mode: 'copy'

  input:
    tuple val(org), path(faa)

  output:
    tuple val(org), path("pfam_${org}.out")

  script:
  """
  # Aseguramos conda en PATH
  export PATH=/opt/conda/bin:\$PATH
  echo "[$org] task.cpus=${task.cpus}" >&2
  # Llamamos a hpc_pfam_search igual que en el benchmark
  /opt/conda/share/pfam_scan-1.6-5/pfam_scan.pl \\
    -fasta ${faa} \\
    -dir /db \\
    -cut_ga \\
    -cpu ${task.cpus} \\
    -outfile pfam_${org}.out
  """
}
