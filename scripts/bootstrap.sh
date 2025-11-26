#!/usr/bin/env bash
set -euo pipefail

# FASTA toy
mkdir -p assets/toys/Ecoli
cat > assets/toys/Ecoli/proteins.faa << 'FASTA'
>protA
MKKLLTLALLALATAMA
>protB
MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDR
>protC
MTEITAAMVKELRESTGAGMMDCKNALSETQ
FASTA

# nextflow.config
cat > nextflow.config << 'CFG'
params {
  input_dir = "assets/toys"
  outdir    = "results"
  organisms = ["Ecoli"]
}
profiles {
  local {
    process.executor = 'local'
    docker.enabled = false
  }
}
CFG

# main.nf
cat > main.nf << 'NF'
nextflow.enable.dsl=2

include { INGEST } from './modules/ingest.nf'
include { COUNT  } from './modules/count_proteins.nf'
include { REPORT } from './modules/report.nf'

workflow {
  Channel.fromList(params.organisms).set { orgs }
  ing = INGEST(orgs, params.input_dir)   // encuentra FASTA
  cnt = COUNT(ing)                        // cuenta cabeceras '>'
  REPORT(cnt, params.outdir)              // escribe summary
}
NF

# modules/ingest.nf
cat > modules/ingest.nf << 'NF'
process INGEST {
  tag { org }
  input:
    val org
    val base_dir
  output:
    tuple val(org), path("proteins.faa")
  script:
  '''
  ln -sf ${base_dir}/${org}/proteins.faa proteins.faa
  '''
}
NF

# modules/count_proteins.nf
cat > modules/count_proteins.nf << 'NF'
process COUNT {
  tag { org }
  input:
    tuple val(org), path(faa)
  output:
    tuple val(org), path("counts_${org}.txt")
  script:
  '''
  grep -c '^>' ${faa} | awk '{print "proteins=" $1}' > counts_${org}.txt
  '''
}
NF

# modules/report.nf
cat > modules/report.nf << 'NF'
process REPORT {
  publishDir "${params.outdir}", mode: 'copy'
  input:
    tuple val(org), path(counts)
  output:
    path "summary.txt"
  script:
  '''
  echo "Summary dom-archflow" > summary.txt
  echo "${org}: $(cat ${counts})" >> summary.txt
  '''
}
NF

