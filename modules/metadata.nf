process METADATA {
  tag { org }
  cpus 1
  maxForks 1   // para no spamear NCBI con varias llamadas a la vez
  publishDir "${params.outdir}/metadata", mode: 'copy'

  input:
    tuple val(org), path(faa), path(gff)

  output:
    tuple val(org), path("metadata_${org}.csv")

  script:
  """
  export PATH=/opt/conda/bin:\$PATH

  # 👇 Cambia este numerito cuando quieras invalidar la caché de METADATA
  echo "METADATA_SCRIPT_VERSION v2" >&2

  python ${projectDir}/scripts/build_metadata.py \\
    --proteome-id ${org} \\
    --faa ${faa} \\
    --gff ${gff} \\
    --out metadata_${org}.csv \\
    --debug
  """
}
