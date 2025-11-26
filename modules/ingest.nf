// modules/ingest.nf
nextflow.enable.dsl=2

process INGEST {
  tag { org }

  input:
    tuple val(org), val(base_dir)

  output:
    tuple val(org), path("proteins.faa"), path("genes.gff")

  script:
  """
  ln -sf ${base_dir}/${org}/proteins.faa proteins.faa
  ln -sf ${base_dir}/${org}/genes.gff  genes.gff
  """
}
