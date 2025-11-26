process PFAM_SCAN {
  cpus params.pfam_cpus
  maxForks params.pfam_maxForks
  tag { org }
 // publishDir "${params.outdir}/pfam", mode: 'copy'
  publishDir "${params.outdir}/pfam_${params.pfam_label}", mode: 'copy'

  input:
    tuple val(org), path(faa)

  output:
    tuple val(org), path("pfam_${org}.out")

  script:
  """
  # Aseguramos hmmscan en el PATH
  export PATH=/opt/conda/bin:\$PATH

  # pfam_scan con su ruta y librerías
  perl -I /opt/conda/share/pfam_scan-1.6-5 \\
    /opt/conda/share/pfam_scan-1.6-5/pfam_scan.pl \\
      -fasta ${faa} \\
      -dir /db \\
      -cut_ga \\
      -outfile pfam_${org}.out \\
      -cpu ${task.cpus}
  """
}