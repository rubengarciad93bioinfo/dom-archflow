#!/usr/bin/env python3
import argparse
from pathlib import Path

import math
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # backend no interactivo
import matplotlib.pyplot as plt
import networkx as nx   # NEW


def plot_frac_domains_per_proteome(df: pd.DataFrame, out_path: Path):
    required = {"proteome_id", "frac_proteins_with_domains"}
    if not required.issubset(df.columns):
        print(f"[VIS_QC] Saltando frac_domains: faltan columnas {required - set(df.columns)}")
        return

    df_plot = df.copy()
    df_plot = df_plot.sort_values("proteome_id")

    x = df_plot["proteome_id"].astype(str).tolist()
    y = df_plot["frac_proteins_with_domains"].astype(float).tolist()

    plt.figure(figsize=(8, 4))
    plt.bar(x, y)
    plt.ylabel("Frac. proteínas con dominios")
    plt.ylim(0, 1.05)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"[VIS_QC] Guardado {out_path}")


def plot_entropy_scatter(df: pd.DataFrame, out_path: Path):
    required = {"proteome_id", "entropy_domains", "entropy_architectures"}
    if not required.issubset(df.columns):
        print(f"[VIS_QC] Saltando entropy_scatter: faltan columnas {required - set(df.columns)}")
        return

    df_plot = df.copy()
    df_plot = df_plot.sort_values("proteome_id")

    x = df_plot["entropy_domains"].astype(float)
    y = df_plot["entropy_architectures"].astype(float)

    sizes = None
    if "n_architectures_distinct" in df_plot.columns:
        a = df_plot["n_architectures_distinct"].astype(float)
        max_a = a.max() if len(a) else 0.0
        if max_a > 0:
            sizes = 20 + 80 * (a / max_a)
        else:
            sizes = 40
    else:
        sizes = 40

    plt.figure(figsize=(6, 5))
    plt.scatter(x, y, s=sizes)
    for _, row in df_plot.iterrows():
        plt.text(
            row["entropy_domains"],
            row["entropy_architectures"],
            str(row["proteome_id"]),
            fontsize=7,
            ha="center",
            va="center",
        )
    plt.xlabel("Entropía dominios")
    plt.ylabel("Entropía arquitecturas")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"[VIS_QC] Guardado {out_path}")


def plot_domain_type_composition(df: pd.DataFrame, out_path: Path):
    base = "frac_domains_type_"
    cols = {
        "Domain": base + "Domain",
        "Family": base + "Family",
        "Repeat": base + "Repeat",
        "Motif":  base + "Motif",
    }
    if not set(cols.values()).issubset(df.columns):
        print(f"[VIS_QC] Saltando domain_type_composition: faltan columnas de tipo dominio")
        return

    df_plot = df.copy()
    df_plot = df_plot.sort_values("proteome_id")
    proteomes = df_plot["proteome_id"].astype(str).tolist()

    import numpy as np

    vals = {t: df_plot[c].fillna(0.0).astype(float).values for t, c in cols.items()}
    x = np.arange(len(proteomes))
    bottom = np.zeros(len(proteomes))

    plt.figure(figsize=(8, 4))
    for tipo in ["Domain", "Family", "Repeat", "Motif"]:
        v = vals[tipo]
        plt.bar(x, v, bottom=bottom, label=tipo)
        bottom = bottom + v

    plt.xticks(x, proteomes, rotation=45, ha="right")
    plt.ylabel("Proporción de dominios por tipo")
    plt.ylim(0, 1.05)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"[VIS_QC] Guardado {out_path}")


def plot_coocc_hist_and_top(edges: pd.DataFrame, hist_out: Path, top_out: Path, top_n: int = 20):
    if edges.empty:
        print("[VIS_QC] Edges vacío, saltando co-ocurrencias")
        return

    weight_col = None
    for c in ["weight_proteins", "weight_architectures", "weight"]:
        if c in edges.columns:
            weight_col = c
            break
    if weight_col is None:
        num_cols = [c for c in edges.columns if pd.api.types.is_numeric_dtype(edges[c])]
        if not num_cols:
            print("[VIS_QC] Ninguna columna numérica encontrada en edges; no se puede hacer histograma")
            return
        weight_col = num_cols[0]

    w = edges[weight_col].dropna().astype(float)
    if w.empty:
        print("[VIS_QC] Columna de pesos vacía")
        return

    # Histograma
    plt.figure(figsize=(6, 4))
    plt.hist(w, bins=30)
    plt.xlabel(f"Peso ({weight_col})")
    plt.ylabel("Frecuencia")
    plt.tight_layout()
    plt.savefig(hist_out, dpi=150)
    plt.close()
    print(f"[VIS_QC] Guardado {hist_out}")

    # Top pairs
    dom1_col = None
    dom2_col = None
    for c in edges.columns:
        if c.lower() in ("dom1", "domain1", "pfam1"):
            dom1_col = c
        if c.lower() in ("dom2", "domain2", "pfam2"):
            dom2_col = c
    if dom1_col is None or dom2_col is None:
        print("[VIS_QC] No se encontraron columnas dom1/dom2; saltando top_pairs")
        return

    top_df = edges.sort_values(weight_col, ascending=False).head(top_n).copy()
    top_df["pair"] = top_df[dom1_col].astype(str) + "-" + top_df[dom2_col].astype(str)

    plt.figure(figsize=(8, 4))
    plt.barh(top_df["pair"], top_df[weight_col].astype(float))
    plt.gca().invert_yaxis()
    plt.xlabel(weight_col)
    plt.tight_layout()
    plt.savefig(top_out, dpi=150)
    plt.close()
    print(f"[VIS_QC] Guardado {top_out}")


# NEW: dibujar un pequeño grafo de co-ocurrencias
def plot_coocc_graph(edges: pd.DataFrame, out_path: Path, top_n: int = 50):
    """
    Construye un grafo a partir de domain_cooccurrence_edges.csv y dibuja
    un subgrafo con las top N aristas por peso.
    """
    if edges.empty:
        print("[VIS_QC] Edges vacío, no se dibuja grafo")
        return

    # Columnas de dominios
    dom1_col = None
    dom2_col = None
    for c in edges.columns:
        lc = c.lower()
        if lc in ("dom1", "domain1", "pfam1"):
            dom1_col = c
        if lc in ("dom2", "domain2", "pfam2"):
            dom2_col = c
    if dom1_col is None or dom2_col is None:
        print("[VIS_QC] No se encontraron columnas dom1/dom2; no se dibuja grafo")
        return

    # Columna de peso
    weight_col = None
    for c in ["weight_proteins", "weight_architectures", "weight"]:
        if c in edges.columns:
            weight_col = c
            break
    if weight_col is None:
        num_cols = [c for c in edges.columns if pd.api.types.is_numeric_dtype(edges[c])]
        if not num_cols:
            print("[VIS_QC] Ninguna columna numérica encontrada en edges; no se dibuja grafo")
            return
        weight_col = num_cols[0]

    df = edges[[dom1_col, dom2_col, weight_col]].dropna().copy()
    df[weight_col] = df[weight_col].astype(float)

    # Top N aristas más pesadas
    df = df.sort_values(weight_col, ascending=False).head(top_n)

    G = nx.Graph()
    for _, row in df.iterrows():
        u = str(row[dom1_col])
        v = str(row[dom2_col])
        w = float(row[weight_col])
        G.add_edge(u, v, weight=w)

    if G.number_of_edges() == 0:
        print("[VIS_QC] Grafo sin aristas después de filtrar; no se dibuja")
        return

    # Usar la componente conexa más grande
    largest_cc = max(nx.connected_components(G), key=len)
    H = G.subgraph(largest_cc).copy()

    degrees = dict(H.degree())
    max_deg = max(degrees.values()) if degrees else 1

    # Layout
    pos = nx.spring_layout(H, k=0.5, iterations=100)

    # tamaños de nodo ~ grado
    node_sizes = [200 + 800 * (degrees[n] / max_deg) for n in H.nodes()]
    edge_widths = [0.5 + 3.0 * (d["weight"] / df[weight_col].max()) for _, _, d in H.edges(data=True)]

    plt.figure(figsize=(8, 6))
    nx.draw_networkx_edges(H, pos, width=edge_widths, alpha=0.6)
    nx.draw_networkx_nodes(H, pos, node_size=node_sizes)
    # etiquetas solo para nodos con grado >= 2 para no hacer una sopa de letras
    labels = {n: n for n, deg in degrees.items() if deg >= 2}
    nx.draw_networkx_labels(H, pos, labels=labels, font_size=7)

    plt.axis("off")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"[VIS_QC] Guardado {out_path}")


def main():
    ap = argparse.ArgumentParser(description="Visualizaciones QC para dom-archflow")
    ap.add_argument("--features-per-proteome", required=True)
    ap.add_argument("--features-global", required=True)
    ap.add_argument("--edges", required=True)
    args = ap.parse_args()

    fpp = Path(args.features_per_proteome)
    fg  = Path(args.features_global)
    ed  = Path(args.edges)

    print(f"[VIS_QC] Leyendo {fpp}")
    per_df = pd.read_csv(fpp)

    print(f"[VIS_QC] Leyendo {fg}")
    try:
        global_df = pd.read_csv(fg)
    except Exception as e:
        print(f"[VIS_QC] Aviso: no se pudo leer features_global ({e}), continúo sin usarlo")
        global_df = None

    print(f"[VIS_QC] Leyendo {ed}")
    try:
        edges_df = pd.read_csv(ed)
    except Exception as e:
        print(f"[VIS_QC] Aviso: no se pudo leer edges ({e}), edgelist vacía")
        edges_df = pd.DataFrame()

    plot_frac_domains_per_proteome(per_df, Path("qc_frac_domains_per_proteome.png"))
    plot_entropy_scatter(per_df, Path("qc_entropy_scatter.png"))
    plot_domain_type_composition(per_df, Path("qc_domain_type_composition.png"))
    plot_coocc_hist_and_top(
        edges_df,
        Path("qc_coocc_weight_hist.png"),
        Path("qc_top_pairs.png"),
    )
    plot_coocc_graph(edges_df, Path("qc_coocc_graph.png"))

    # === Índice HTML muy simple ===
    html_path = Path("qc_report.html")
    imgs = [
        "qc_frac_domains_per_proteome.png",
        "qc_entropy_scatter.png",
        "qc_domain_type_composition.png",
        "qc_coocc_weight_hist.png",
        "qc_top_pairs.png",
        "qc_coocc_graph.png",
    ]
    with html_path.open("w") as f:
        f.write("<html><body>\n")
        f.write("<h1>dom-archflow QC report</h1>\n")
        for img in imgs:
            if Path(img).exists():
                f.write(f"<h2>{img}</h2>\n")
                f.write(f"<img src='{img}' style='max-width:900px;'><br/>\n")
        f.write("</body></html>\n")
    print("[VIS_QC] Guardado qc_report.html")


if __name__ == "__main__":
    main()
