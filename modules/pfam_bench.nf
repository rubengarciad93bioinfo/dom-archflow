// modules/pfam_bench.nf

/*
 * 1) pfam_scan "gordo": usa muchos cores por tarea,
 *    pero dejamos que pfam_scan decida los hilos (sin -cpu).
 */
process PFAM_SCAN_GORDO {
  tag { org }
  cpus params.pfam_cpus_all
  maxForks params.pfam_maxForks_all
  publishDir "${params.outdir}/pfam/benchmark/gordo", mode: 'copy'

  input:
    tuple val(org), path(faa)

  output:
    tuple val(org), path("pfam_${org}.gordo.out")

  script:
  """
  export PATH=/opt/conda/bin:\$PATH

  perl -I /opt/conda/share/pfam_scan-1.6-5 \\
    /opt/conda/share/pfam_scan-1.6-5/pfam_scan.pl \\
      -fasta ${faa} \\
      -dir /db \\
      -cut_ga \\
      -outfile pfam_${org}.gordo.out
  """
}

/*
 * 2) pfam_scan "many small": pocas CPUs por tarea, varias en paralelo.
 */
process PFAM_SCAN_MANY_SMALL {
  tag { org }
  cpus params.pfam_cpus_small
  maxForks params.pfam_maxForks_many
  publishDir "${params.outdir}/pfam/benchmark/many_small", mode: 'copy'

  input:
    tuple val(org), path(faa)

  output:
    tuple val(org), path("pfam_${org}.many_small.out")

  script:
  """
  export PATH=/opt/conda/bin:\$PATH

  perl -I /opt/conda/share/pfam_scan-1.6-5 \\
    /opt/conda/share/pfam_scan-1.6-5/pfam_scan.pl \\
      -fasta ${faa} \\
      -dir /db \\
      -cut_ga \\
      -outfile pfam_${org}.many_small.out \\
      -cpu ${task.cpus}
  """
}

/*
 * 3a) hpc_pfam_search "normal": más conservador
 */
process PFAM_SCAN_HPC {
  tag { org }
  cpus params.pfam_cpus_hpc_std
  maxForks params.pfam_maxForks_hpc_std
  publishDir "${params.outdir}/pfam/benchmark/hpc_std", mode: 'copy'

  input:
    tuple val(org), path(faa)

  output:
    tuple val(org), path("pfam_${org}.hpc.out")

  script:
  """
  export PATH=/opt/conda/bin:\$PATH

  /opt/hpc_pfam_search/hpc_pfam_search.pl \\
    -fasta ${faa} \\
    -dir /db \\
    -cut_ga \\
    -outfile pfam_${org}.hpc.out \\
    -cpu ${task.cpus}
  """
}

/*
 * 3b) hpc_pfam_search "a saco": menos CPUs por tarea pero más tareas a la vez
 */
process PFAM_SCAN_HPC_MAX {
  tag { org }
  cpus params.pfam_cpus_hpc_max
  maxForks params.pfam_maxForks_hpc_max
  publishDir "${params.outdir}/pfam/benchmark/hpc_max", mode: 'copy'

  input:
    tuple val(org), path(faa)

  output:
    tuple val(org), path("pfam_${org}.hpc_max.out")

  script:
  """
  export PATH=/opt/conda/bin:\$PATH

  /opt/hpc_pfam_search/hpc_pfam_search.pl \\
    -fasta ${faa} \\
    -dir /db \\
    -cut_ga \\
    -outfile pfam_${org}.hpc_max.out \\
    -cpu ${task.cpus}
  """
}

