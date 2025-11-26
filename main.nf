nextflow.enable.dsl=2

include { INGEST        } from './modules/ingest.nf'
include { COUNT         } from './modules/count_proteins.nf'
include { MERGE_COUNTS  } from './modules/merge_counts.nf'
include { PFAM_SCAN     } from './modules/pfam_scan_3_hpc.nf'
include { DOM_ARCH      } from './modules/dom_arch.nf'
include { METADATA      } from './modules/metadata.nf'
include { MERGE_METADATA} from './modules/merge_metadata.nf'
include { MERGE_ARCH    } from './modules/merge_arch.nf'
// 🚨 nuevos módulos:
include { FEATURES      } from './modules/features.nf'
include { DOMAIN_COOCC  } from './modules/domains_coocc.nf'
include { VIS_QC        } from './modules/vis_qc.nf'
include { PHYLO_MATRICES } from './modules/phylo_matrices.nf'
include { IQTREE_RUN    } from './modules/iqtree.nf'
params.subset_arch  = params.subset_arch ?: null
params.iqtree_cpus  = params.iqtree_cpus ?: 4

workflow {
  Channel.fromList(params.organisms)
         .map { org -> tuple(org, params.input_dir) }
         .set { ingest_in }

  // 1) INGEST: (org, faa, gff)
  ing = INGEST(ingest_in)

  // 2) COUNT + MERGE_COUNTS (usa faa)
  cnt   = COUNT( ing.map { org, faa, gff -> tuple(org, faa) } )
  files = cnt.map { it[1] }
  counts_file = MERGE_COUNTS( files.collect() )

  // 3) METADATA: un CSV por organismo
  meta_per_org = METADATA(ing)    // (org, metadata_org.csv)
  meta_files   = meta_per_org.map { it[1] }

  // 4) MERGE_METADATA -> metadata_all.csv
  meta_all_file = MERGE_METADATA( meta_files.collect() )

  // 5) Leer metadata_all.csv como filas
  meta_all = meta_all_file
    .splitCsv(header: true)
    .map { row ->
      tuple(row.proteome_id, row.tax_id, row.taxonomy, row.completeness)
    }

  // 6) PFAM_SCAN: sólo faa
  pfam = PFAM_SCAN( ing.map { org, faa, gff -> tuple(org, faa) } )
  // pfam: (org, pfam_out)

  // 7) join pfam + metadata por proteome_id = org
  pfam_with_meta = pfam.join(meta_all).map { org, pfam_out, tax_id, taxonomy, completeness ->
    tuple(org, pfam_out, tax_id, taxonomy, completeness)
  }

  dom_arch_out = DOM_ARCH(pfam_with_meta)   // ya lo tienes

  arch_csvs = dom_arch_out.map { it[1] }    // ajusta al esquema real
  all_arch_file = MERGE_ARCH( arch_csvs.collect() )
  
  // 10) FEATURES (usa all_architectures.csv + summary.txt)
  features_in = all_arch_file.combine(counts_file)

  def (features_per_ch, features_global_ch) = FEATURES(features_in)

  // 11) GRAFO DE CO-OCURRENCIA DE DOMINIOS
  def (coocc_edges, coocc_gexf) = DOMAIN_COOCC(all_arch_file)

    // 12) VIS_QC recibe 3 canales de entrada (workflow)
  VIS_QC(features_per_ch, features_global_ch, coocc_edges)

  // 13) Elegir qué arquitecturas usar para las matrices filogenéticas:
  //     - si params.subset_arch está definido -> usar ese CSV
  //     - si no, usar el all_arch_file del propio pipeline
  def arch_for_phylo_ch

  if( params.subset_arch ) {
    log.info "Usando subset para phylo: ${params.subset_arch}"
    arch_for_phylo_ch = Channel.of( file(params.subset_arch) )
  }
  else {
    log.info "Usando all_arch_file completo para phylo"
    arch_for_phylo_ch = all_arch_file
  }

  // 14) Matrices filogenéticas (dominios, pares, dom+pares, arqs completas)
  //     -> devuelve un canal con TODOS los .phylip creados
  phylo_phylip_ch = PHYLO_MATRICES( arch_for_phylo_ch )

  // 15) IQ-TREE PARA TODOS LOS PHYLIP
  IQTREE_RUN( phylo_phylip_ch.flatten() )
}