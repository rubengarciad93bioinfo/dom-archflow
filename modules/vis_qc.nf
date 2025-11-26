process VIS_QC{
  tag "vis_qc"
  cpus 1
  publishDir "${params.outdir}/qc", mode: 'copy'

  input:
    path per_feat      // features_per_proteome.csv
    path glob_feat     // features_global.csv
    path edge_file     // domain_cooccurrence_edges.csv

  output:
    path "qc_report.html"
    path "qc_frac_domains_per_proteome.png"
    path "qc_entropy_scatter.png"
    path "qc_domain_type_composition.png"
    path "qc_coocc_weight_hist.png"
    path "qc_top_pairs.png"
    path "qc_coocc_graph.png"   // ← solo una vez

  script:
  """
  export PATH=/opt/conda/bin:\$PATH

  python ${projectDir}/scripts/vis_qc.py \\
    --features-per-proteome ${per_feat} \\
    --features-global ${glob_feat} \\
    --edges ${edge_file}
  """
}
