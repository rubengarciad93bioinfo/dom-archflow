#!/usr/bin/env python3
import argparse
from pathlib import Path
from collections import defaultdict
import pandas as pd
import numpy as np


def build_binary_matrix(feature_map):
    """
    feature_map: dict[proteome_id] -> set(features)
    Devuelve DataFrame (rows=proteome_id, cols=features) con 0/1
    (solo las features, sin metadata).
    """
    proteomes = sorted(feature_map.keys())
    all_features = sorted({f for s in feature_map.values() for f in s})
    if not all_features:
        # matriz vacía: devolvemos df con solo índice
        return pd.DataFrame(index=proteomes)

    mat = np.zeros((len(proteomes), len(all_features)), dtype=int)
    feat_index = {f: i for i, f in enumerate(all_features)}

    for r, pid in enumerate(proteomes):
        for f in feature_map[pid]:
            j = feat_index[f]
            mat[r, j] = 1

    df = pd.DataFrame(mat, index=proteomes, columns=all_features)
    return df


def filter_rare_features(mat, min_present=2, max_present=None):
    """
    mat: DataFrame (rows=proteome_id, cols=features) con 0/1
    min_present: número mínimo de taxa en los que debe aparecer una feature.
    max_present: opcional, máximo de taxa donde puede estar (para quitar casi-invariantes).

    Devuelve una matriz con solo las columnas que pasan los filtros.
    """
    if mat.empty:
        return mat

    col_sums = mat.sum(axis=0)

    mask = col_sums >= min_present
    if max_present is not None:
        mask &= (col_sums <= max_present)

    kept = mask[mask].index.tolist()
    print(f"[PHYLO] {mat.shape[1]} columnas -> {len(kept)} tras filtrar (min_present={min_present}, max_present={max_present})")

    return mat[kept]


def write_phylip(df, out_path, id_col="proteome_id"):
    """
    Escribe una matriz binaria en formato PHYLIP secuencial.
    Usa id_col como nombre del taxón, y todas las columnas 0/1 como caracteres.
    """
    non_feature_cols = {id_col, "tax_id", "taxonomy", "completeness"}
    feature_cols = [c for c in df.columns if c not in non_feature_cols]

    n_taxa = len(df)
    n_chars = len(feature_cols)

    print(f"[PHYLO] Escribiendo PHYLIP {out_path} con {n_taxa} taxa y {n_chars} caracteres")

    with open(out_path, "w") as fh:
        fh.write(f"{n_taxa} {n_chars}\n")
        for _, row in df.iterrows():
            pid = str(row[id_col])
            seq = "".join(str(int(row[c])) for c in feature_cols) if n_chars > 0 else ""
            fh.write(f"{pid}\t{seq}\n")


def main():
    ap = argparse.ArgumentParser(
        description="Construye matrices filogenéticas (dominios, pares, dominios+pares, arquitecturas completas) a partir de all_architectures.csv"
    )
    ap.add_argument("--all-architectures", required=True,
                    help="Fichero all_architectures.csv (o subset_architectures.csv)")
    ap.add_argument("--min-present", type=int, default=2,
                    help="Nº mínimo de taxa en los que debe estar presente una feature (por defecto 2)")
    args = ap.parse_args()

    all_arch_path = Path(args.all_architectures)
    print(f"[PHYLO] Leyendo arquitecturas desde: {all_arch_path}")
    df = pd.read_csv(all_arch_path)

    required_cols = {"proteome_id", "doms"}
    if not required_cols.issubset(df.columns):
        missing = required_cols - set(df.columns)
        raise SystemExit(f"Faltan columnas en {all_arch_path}: {missing}")

    # --- metadata por proteoma ---
    meta_cols = [c for c in ["proteome_id", "tax_id", "taxonomy", "completeness"] if c in df.columns]
    meta = (
        df[meta_cols]
        .drop_duplicates(subset=["proteome_id"])
        .set_index("proteome_id")
        .sort_index()
    )
    print(f"[PHYLO] Nº de proteomas (taxa): {meta.shape[0]}")

    # Mapas: proteome_id -> set(features)
    dom_features = defaultdict(set)
    pair_features = defaultdict(set)
    arch_features = defaultdict(set)

    for _, row in df.iterrows():
        pid = row["proteome_id"]
        doms_str = str(row["doms"])
        if doms_str == "nan" or not doms_str:
            continue
        dom_list = doms_str.split("-")

        # dominios individuales
        for d in dom_list:
            dom_features[pid].add(d)

        # pares adyacentes
        if len(dom_list) >= 2:
            for i in range(len(dom_list) - 1):
                pair = f"{dom_list[i]}-{dom_list[i+1]}"
                pair_features[pid].add(pair)

        # arquitectura completa (cadena tal cual)
        arch_features[pid].add(doms_str)

    # union dominios + pares
    dompair_features = defaultdict(set)
    all_pids = set(dom_features.keys()) | set(pair_features.keys()) | set(arch_features.keys())
    for pid in all_pids:
        doms = dom_features.get(pid, set())
        pairs = pair_features.get(pid, set())
        dompair_features[pid] = doms | pairs

    # --- construir matrices binarias ---
    print("[PHYLO] Construyendo matriz de dominios...")
    mat_domains = build_binary_matrix(dom_features)
    print("[PHYLO] Construyendo matriz de pares...")
    mat_pairs = build_binary_matrix(pair_features)
    print("[PHYLO] Construyendo matriz dominios+pares...")
    mat_dompairs = build_binary_matrix(dompair_features)
    print("[PHYLO] Construyendo matriz de arquitecturas completas...")
    mat_arch = build_binary_matrix(arch_features)

    # --- filtrar columnas raras: presentes en menos de min_present taxa ---
    min_present = args.min_present
    mat_domains = filter_rare_features(mat_domains, min_present=min_present)
    mat_pairs = filter_rare_features(mat_pairs, min_present=min_present)
    mat_dompairs = filter_rare_features(mat_dompairs, min_present=min_present)
    mat_arch = filter_rare_features(mat_arch, min_present=min_present)

    # añadir metadata delante (si existe)
    def attach_meta(mat, name):
        if mat.empty:
            print(f"[PHYLO] MATRIZ VACÍA después de filtrar para {name}, se devuelve solo metadata")
            return meta.reset_index()
        mat = mat.copy()
        mat.index.name = "proteome_id"
        mat = mat.reset_index()
        if not meta.empty:
            mat = mat.merge(meta.reset_index(), on="proteome_id", how="left")
            # reordenar: proteome_id, tax_id, taxonomy, completeness, features...
            front = [c for c in ["proteome_id", "tax_id", "taxonomy", "completeness"] if c in mat.columns]
            feat_cols = [c for c in mat.columns if c not in front]
            mat = mat[front + feat_cols]
        print(f"[PHYLO] {name}: {mat.shape[0]} taxa x {mat.shape[1]-len(front)} caracteres (después de filtrar)")
        return mat

    mat_domains = attach_meta(mat_domains, "domains_only")
    mat_pairs = attach_meta(mat_pairs, "pairs_only")
    mat_dompairs = attach_meta(mat_dompairs, "domains_pairs")
    mat_arch = attach_meta(mat_arch, "full_arch")

    # --- escribir CSV ---
    mat_domains.to_csv("phylo_domains_only.csv", index=False)
    mat_pairs.to_csv("phylo_pairs_only.csv", index=False)
    mat_dompairs.to_csv("phylo_domains_pairs.csv", index=False)
    mat_arch.to_csv("phylo_full_arch.csv", index=False)

    # --- escribir PHYLIP ---
    write_phylip(mat_domains, "phylo_domains_only.phylip")
    write_phylip(mat_pairs, "phylo_pairs_only.phylip")
    write_phylip(mat_dompairs, "phylo_domains_pairs.phylip")
    write_phylip(mat_arch, "phylo_full_arch.phylip")


if __name__ == "__main__":
    main()
