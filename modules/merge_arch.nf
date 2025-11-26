process MERGE_ARCH {
  cpus 1
  publishDir "${params.outdir}/architectures", mode: 'copy'

  input:
    path csvs

  output:
    path "all_architectures.csv"

  script:
  """
  # Tomamos la cabecera del primer fichero
  head -n 1 ${csvs[0]} > all_architectures.csv

  # Añadimos el contenido (sin cabecera) de todos los CSV
  for f in ${csvs}; do
    tail -n +2 "\$f" >> all_architectures.csv
  done
  """
}
