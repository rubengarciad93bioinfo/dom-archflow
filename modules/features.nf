// modules/features.nf
nextflow.enable.dsl=2

process FEATURES {
  cpus 1
  publishDir "${params.outdir}/features", mode: 'copy'

  input:
    tuple path(architectures), path(counts)   // counts = summary.txt

  output:
    path "features_per_proteome.csv"
    path "features_global.csv"

  script:
  """
  python << 'EOF'
  import pandas as pd
  import math
  from collections import Counter

  # all_architectures.csv
  arch = pd.read_csv("${architectures}")

  # summary.txt -> DataFrame con proteome_id, n_proteins_total
  proteomes = []
  with open("${counts}", "r") as fh:
      for line in fh:
          line = line.strip()
          if not line:
              continue
          if line.startswith("Summary"):
              continue
          if ":" not in line or "proteins=" not in line:
              continue
          left, right = line.split(":", 1)
          pid = left.strip()
          # right ~ ' proteins=4042'
          try:
              n_str = right.split("proteins=")[1]
              n_str = n_str.strip()
              # por si acaso hubiera más texto, nos quedamos con dígitos
              n_digits = "".join(ch for ch in n_str if ch.isdigit())
              n_total = int(n_digits)
          except Exception:
              continue
          proteomes.append(dict(proteome_id=pid, n_proteins_total=n_total))

  counts_df = pd.DataFrame(proteomes)

  def shannon(counts_iter):
      total = sum(counts_iter)
      if total == 0:
          return 0.0
      H = 0.0
      for c in counts_iter:
          if c <= 0:
              continue
          p = c / total
          H -= p * math.log2(p)
      return H

  rows = []

  for pid, sub in arch.groupby("proteome_id"):
      # nº proteínas con dominios
      n_with = sub["seqId"].nunique()

      # nº proteínas total: del summary; si falta, usamos n_with
      row_c = counts_df[counts_df["proteome_id"] == pid]
      if len(row_c) == 1:
          n_total = int(row_c.iloc[0]["n_proteins_total"])
      else:
          n_total = n_with

      frac_with = n_with / n_total if n_total > 0 else 0.0

      # dominios y arquitecturas
      dom_lists = sub["doms"].dropna().astype(str).str.split("-")
      all_doms = [d for lst in dom_lists for d in lst]
      dom_counts = Counter(all_doms)
      n_dom_dist = len(dom_counts)
      n_arch_dist = sub["doms"].nunique()

      # pares adyacentes
      adj_pairs = []
      for lst in dom_lists:
          if len(lst) < 2:
              continue
          adj_pairs.extend(zip(lst, lst[1:]))
      pair_counts = Counter(adj_pairs)
      n_adj_total = sum(pair_counts.values())
      n_adj_dist  = len(pair_counts)

      # tipos
      types_lists = sub["types"].dropna().astype(str).str.split("-")
      all_types = [t for lst in types_lists for t in lst]
      tcounts = Counter(all_types)
      total_doms = len(all_doms)

      def tstat(t):
          c = tcounts.get(t, 0)
          return c, (c / total_doms if total_doms > 0 else 0.0)

      n_dom_dom, frac_dom_dom = tstat("Domain")
      n_dom_fam, frac_dom_fam = tstat("Family")
      n_dom_rep, frac_dom_rep = tstat("Repeat")
      n_dom_mot, frac_dom_mot = tstat("Motif")

      # entropías
      ent_domains = shannon(dom_counts.values())
      arch_counts = Counter(sub["doms"].dropna().astype(str))
      ent_arch = shannon(arch_counts.values())

      tax_id   = sub["tax_id"].iloc[0]
      taxonomy = sub["taxonomy"].iloc[0]

      rows.append(dict(
          proteome_id = pid,
          tax_id = tax_id,
          taxonomy = taxonomy,
          n_proteins_total = n_total,
          n_proteins_with_domains = n_with,
          frac_proteins_with_domains = frac_with,
          n_domains_distinct = n_dom_dist,
          n_architectures_distinct = n_arch_dist,
          n_adj_pairs_total = n_adj_total,
          n_adj_pairs_distinct = n_adj_dist,
          n_domains_type_Domain = n_dom_dom,
          frac_domains_type_Domain = frac_dom_dom,
          n_domains_type_Family = n_dom_fam,
          frac_domains_type_Family = frac_dom_fam,
          n_domains_type_Repeat = n_dom_rep,
          frac_domains_type_Repeat = frac_dom_rep,
          n_domains_type_Motif = n_dom_mot,
          frac_domains_type_Motif = frac_dom_mot,
          entropy_domains = ent_domains,
          entropy_architectures = ent_arch,
      ))

  per_df = pd.DataFrame(rows)
  per_df.to_csv("features_per_proteome.csv", index=False)

  # ---- Global ----
  # nº proteínas total global: suma de todos los proteomas en el summary
  if not counts_df.empty:
      n_total_global = int(counts_df["n_proteins_total"].sum())
  else:
      n_total_global = int(arch["seqId"].nunique())

  n_with_global = int(arch["seqId"].nunique())
  frac_with_global = n_with_global / n_total_global if n_total_global > 0 else 0.0

  dom_lists_all = arch["doms"].dropna().astype(str).str.split("-")
  all_doms_all = [d for lst in dom_lists_all for d in lst]
  dom_counts_all = Counter(all_doms_all)
  n_dom_dist_all = len(dom_counts_all)
  n_arch_dist_all = arch["doms"].nunique()

  adj_pairs_all = []
  for lst in dom_lists_all:
      if len(lst) < 2:
          continue
      adj_pairs_all.extend(zip(lst, lst[1:]))
  pair_counts_all = Counter(adj_pairs_all)
  n_adj_total_all = sum(pair_counts_all.values())
  n_adj_dist_all  = len(pair_counts_all)

  types_lists_all = arch["types"].dropna().astype(str).str.split("-")
  all_types_all = [t for lst in types_lists_all for t in lst]
  tcounts_all = Counter(all_types_all)
  total_doms_all = len(all_doms_all)

  def tstat_all(t):
      c = tcounts_all.get(t, 0)
      return c, (c / total_doms_all if total_doms_all > 0 else 0.0)

  n_dom_dom_all, frac_dom_dom_all = tstat_all("Domain")
  n_dom_fam_all, frac_dom_fam_all = tstat_all("Family")
  n_dom_rep_all, frac_dom_rep_all = tstat_all("Repeat")
  n_dom_mot_all, frac_dom_mot_all = tstat_all("Motif")

  ent_domains_all = shannon(dom_counts_all.values())
  arch_counts_all = Counter(arch["doms"].dropna().astype(str))
  ent_arch_all = shannon(arch_counts_all.values())

  global_row = dict(
      n_proteins_total = n_total_global,
      n_proteins_with_domains = n_with_global,
      frac_proteins_with_domains = frac_with_global,
      n_domains_distinct = n_dom_dist_all,
      n_architectures_distinct = n_arch_dist_all,
      n_adj_pairs_total = n_adj_total_all,
      n_adj_pairs_distinct = n_adj_dist_all,
      n_domains_type_Domain = n_dom_dom_all,
      frac_domains_type_Domain = frac_dom_dom_all,
      n_domains_type_Family = n_dom_fam_all,
      frac_domains_type_Family = frac_dom_fam_all,
      n_domains_type_Repeat = n_dom_rep_all,
      frac_domains_type_Repeat = frac_dom_rep_all,
      n_domains_type_Motif = n_dom_mot_all,
      frac_domains_type_Motif = frac_dom_mot_all,
      entropy_domains = ent_domains_all,
      entropy_architectures = ent_arch_all,
  )

  gdf = pd.DataFrame([global_row])
  gdf.to_csv("features_global.csv", index=False)

  EOF
  """
}
