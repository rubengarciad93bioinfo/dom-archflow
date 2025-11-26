process MERGE_METADATA {
  cpus 1
  publishDir "${params.outdir}/metadata", mode: 'copy'

  input:
    path(meta_files)

  output:
    path("metadata_all.csv")

  script:
  """
  # 👇 También forzamos invalidar la caché de MERGE_METADATA
  echo "MERGE_METADATA_SCRIPT_VERSION v2" >&2

  # cabecera del primero
  head -n 1 ${meta_files[0]} > metadata_all.csv

  # resto sin cabecera
  for f in ${meta_files}; do
    tail -n +2 "\$f" >> metadata_all.csv
  done
  """
}
