// benchmark.nf
nextflow.enable.dsl=2

include { INGEST } from './modules/ingest.nf'
include { PFAM_SCAN_GORDO
        ; PFAM_SCAN_MANY_SMALL
        ; PFAM_SCAN_HPC
        ; PFAM_SCAN_HPC_MAX } from './modules/pfam_bench.nf'

workflow {

  // 1) Canal de entrada: (org, base_dir)
  Channel.fromList(params.organisms)
         .map { org -> tuple(org, params.input_dir) }
         .set { ingest_in }

  // 2) INGEST -> (org, faa, gff)
  ing = INGEST(ingest_in)

  // 3) Para el benchmark sólo necesitamos (org, faa)
  fafas = ing.map { org, faa, gff -> tuple(org, faa) }

  // 4) Elegimos QUÉ benchmark correr según params.bench_mode
  if (params.bench_mode == 'gordo') {
    log.info "Benchmark PFAM: modo GORDO (pfam_scan.pl, tarea gorda por proteoma)"
    PFAM_SCAN_GORDO(fafas)
  }
  else if (params.bench_mode == 'many') {
    log.info "Benchmark PFAM: modo MANY_SMALL (pfam_scan.pl, muchas tareas pequeñas en paralelo)"
    PFAM_SCAN_MANY_SMALL(fafas)
  }
  else if (params.bench_mode == 'hpc') {
    log.info "Benchmark PFAM: modo HPC estándar (hpc_pfam_search, moderado)"
    PFAM_SCAN_HPC(fafas)
  }
  else if (params.bench_mode == 'hpc_max') {
    log.info "Benchmark PFAM: modo HPC_MAX (hpc_pfam_search, más paralelización)"
    PFAM_SCAN_HPC_MAX(fafas)
  }
  else {
    throw new IllegalArgumentException("bench_mode desconocido: ${params.bench_mode}")
  }
}
