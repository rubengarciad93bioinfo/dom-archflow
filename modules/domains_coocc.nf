// modules/domain_coocc.nf
nextflow.enable.dsl=2

process DOMAIN_COOCC {
  cpus 1
  publishDir "${params.outdir}/graphs", mode: 'copy'

  input:
    path architectures  // all_architectures.csv

  output:
    path "domain_cooccurrence_edges.csv"
    path "domain_cooccurrence.gexf"

  script:
  """
  python << 'EOF'
  import pandas as pd
  from collections import Counter
  import networkx as nx

  arch = pd.read_csv("${architectures}")

  # mapeo proteína -> proteoma
  prot2prot = dict(zip(arch["seqId"], arch["proteome_id"]))

  pair_protein_counts = Counter()
  pair_proteome_sets  = {}

  for _, row in arch.iterrows():
      seq = row["seqId"]
      prot = prot2prot[seq]
      doms = str(row["doms"]).split("-") if pd.notna(row["doms"]) else []
      if len(doms) < 2:
          continue
      pairs = list(zip(doms, doms[1:]))
      for p in pairs:
          pair = tuple(p)  # ya está en orden de izquierda a derecha
          pair_protein_counts[pair] += 1
          pair_proteome_sets.setdefault(pair, set()).add(prot)

  rows = []
  for pair, c_prot in pair_protein_counts.items():
      c_protomes = len(pair_proteome_sets[pair])
      rows.append(dict(
          dom1 = pair[0],
          dom2 = pair[1],
          weight_proteins = c_prot,
          weight_proteomes = c_protomes,
      ))

  df = pd.DataFrame(rows)
  df.to_csv("domain_cooccurrence_edges.csv", index=False)

  # grafo para Gephi/Cytoscape
  G = nx.Graph()
  for _, r in df.iterrows():
      G.add_edge(
          r["dom1"], r["dom2"],
          weight_proteins = int(r["weight_proteins"]),
          weight_proteomes = int(r["weight_proteomes"])
      )

  nx.write_gexf(G, "domain_cooccurrence.gexf")
  EOF
  """
}
