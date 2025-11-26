process PFAM_SCAN {
  tag { org }

  /*
   * cpus = nº de jobs GNU parallel en paralelo (cada uno hpc_pfam_search -cpu 1)
   * maxForks = cuántos PFAM_SCAN (organismos) a la vez
   */
  cpus params.pfam_cpus
  maxForks params.pfam_maxForks

  publishDir "${params.outdir}/pfam_${params.pfam_label}", mode: 'copy'

  input:
    tuple val(org), path(faa)

  output:
    tuple val(org), path("pfam_${org}.out")

  script:
  """
  # Aseguramos conda en PATH
  export PATH=/opt/conda/bin:\$PATH

  # Enlace cómodo
  ln -sf ${faa} in.fasta

  # Decidir tamaño de chunk
  CHUNKSIZE=\$(decide_chunk_size.sh in.fasta "${task.cpus}" 100)
  echo "[$org] Usando CHUNKSIZE=\$CHUNKSIZE con ${task.cpus} jobs paralelos" >&2

  parallel \\
    --halt now,fail=1 \\
    --joblog log.txt \\
    -j "${task.cpus}" \\
    -N \${CHUNKSIZE} \\
    --recstart '>' \\
    --line-buffer \\
    --cat 'tmp=\$(mktemp chunk_XXXX.out); \\
           /opt/hpc_pfam_search/hpc_pfam_search.pl -fasta "{}" -dir /db -cut_ga -cpu 1 -outfile "\$tmp"; \\
           cat "\$tmp"; \\
           rm "\$tmp"' \\
    < in.fasta \\
    > pfam_${org}.out
  """
}
