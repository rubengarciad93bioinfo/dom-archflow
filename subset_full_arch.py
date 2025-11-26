#!/usr/bin/env python3
"""
Script único para Spyder:

Objetivo:
- Leer combined_with_uabimicrobiia.csv (sin cabecera, 9 columnas).
- Construir metadata por proteoma.
- Mantener TODOS los Uabimicrobiia.
- Añadir el mínimo número de proteomas extra tal que, para cada clase,
  haya al menos MIN_PER_CLASS proteomas (si existen suficientes en el dataset).
- Guardar:
    - subset_min5perclass_metadata.csv
    - subset_min5perclass_architectures.csv
"""

from pathlib import Path
import pandas as pd
import numpy as np

# ===========================
# CONFIGURACIÓN
# ===========================

INFILE = "/Users/rubensiok/uabi_genomes/trees/bacteria/combined_with_uabimicrobiia.csv"

MIN_PER_CLASS = 11        # mínimo deseado por clase
RANDOM_STATE  = 42       # para muestreo reproducible

OUT_PREFIX = f"subset_min{MIN_PER_CLASS}perclass"
META_OUT  = f"{OUT_PREFIX}_metadata.csv"
ARCH_OUT  = f"{OUT_PREFIX}_architectures.csv"


# ===========================
# FUNCIONES
# ===========================

def parse_taxonomy(taxonomy: str):
    """
    Recibe algo tipo:
      d__Bacteria;p__Planctomycetota;c__Uabimicrobiia;o__...
    Devuelve (domain, phylum, class_name).
    """
    domain = "Unknown"
    phylum = "Unknown"
    clazz  = "Unknown"

    if pd.isna(taxonomy):
        return domain, phylum, clazz

    parts = [p.strip() for p in str(taxonomy).split(";") if p.strip()]

    for p in parts:
        if p.startswith("d__"):
            domain = p[3:] or "Unknown"
        elif p.startswith("p__"):
            phylum = p[3:] or "Unknown"
        elif p.startswith("c__"):
            clazz = p[3:] or "Unknown"

    return domain, phylum, clazz


# def load_arch_and_metadata(infile):
#     """
#     Lee combined_with_uabimicrobiia.csv SIN cabecera, con columnas:

#         0: seqId
#         1: proteome_id
#         2: tax_id
#         3: taxonomy
#         4: doms
#         5: names
#         6: clans
#         7: types
#         8: completeness

#     Devuelve:
#         arch_df: tabla completa de arquitecturas
#         meta_df: una fila por proteome_id, con domain/phylum/class/etc.
#     """
#     infile = Path(infile)

#     colnames = [
#         "seqId",
#         "proteome_id",
#         "tax_id",
#         "taxonomy",
#         "doms",
#         "names",
#         "clans",
#         "types",
#         "completeness",
#     ]

#     arch_df = pd.read_csv(infile, header=None, names=colnames)

#     records = []
#     for pid, sub in arch_df.groupby("proteome_id"):
#         tax_id       = str(sub["tax_id"].iloc[0])
#         taxonomy     = sub["taxonomy"].iloc[0]
#         completeness = sub["completeness"].iloc[0]

#         domain, phylum, clazz = parse_taxonomy(taxonomy)

#         records.append(
#             dict(
#                 proteome_id   = pid,
#                 tax_id        = tax_id,
#                 taxonomy      = taxonomy,
#                 completeness  = completeness,
#                 domain        = domain,
#                 phylum        = phylum,
#                 class_name    = clazz,
#             )
#         )

#     meta_df = pd.DataFrame(records)

#     # Flag Uabimicrobiia (por clase o en la cadena taxonomy)
#     meta_df["is_uabi"] = meta_df["class_name"].str.contains(
#         "Uabimicrobiia", case=False, na=False
#     ) | meta_df["taxonomy"].astype(str).str.contains(
#         "Uabimicrobiia", case=False, na=False
#     )

#     return arch_df, meta_df

def load_arch_and_metadata(infile):
    """
    Lee combined_with_uabimicrobiia.csv SIN cabecera, con columnas:

        0: seqId
        1: proteome_id
        2: tax_id
        3: taxonomy
        4: doms
        5: names
        6: clans
        7: types
        8: completeness

    Devuelve:
        arch_df: tabla completa de arquitecturas
        meta_df: una fila por proteome_id, con domain/phylum/class/etc.
    """
    infile = Path(infile)

    colnames = [
        "seqId",
        "proteome_id",
        "tax_id",
        "taxonomy",
        "doms",
        "names",
        "clans",
        "types",
        "completeness",
    ]

    arch_df = pd.read_csv(infile, header=None, names=colnames)

    records = []
    for pid, sub in arch_df.groupby("proteome_id"):
        tax_id       = str(sub["tax_id"].iloc[0])
        taxonomy     = sub["taxonomy"].iloc[0]
        completeness = sub["completeness"].iloc[0]

        domain, phylum, clazz = parse_taxonomy(taxonomy)

        records.append(
            dict(
                proteome_id   = pid,
                tax_id        = tax_id,
                taxonomy      = taxonomy,
                completeness  = completeness,
                domain        = domain,
                phylum        = phylum,
                class_name    = clazz,
            )
        )

    meta_df = pd.DataFrame(records)

    # 💥 NUEVO: quedarnos solo con Bacteria
    meta_df = meta_df[meta_df["domain"] == "Bacteria"].copy()

    # 💥 y filtrar arquitecturas al mismo conjunto de proteomas
    arch_df = arch_df[arch_df["proteome_id"].isin(meta_df["proteome_id"])].copy()

    # Flag Uabimicrobiia (por clase o en la cadena taxonomy)
    meta_df["is_uabi"] = meta_df["class_name"].str.contains(
        "Uabimicrobiia", case=False, na=False
    ) | meta_df["taxonomy"].astype(str).str.contains(
        "Uabimicrobiia", case=False, na=False
    )

    return arch_df, meta_df



def select_min_for_min_per_class(meta_df, min_per_class=5, random_state=42):
    """
    Devuelve:

        meta_subset, selected_ids, class_counts_subset, classes_below_min

    Reglas:
      - Todos los proteomas con is_uabi=True se seleccionan siempre.
      - Para cada (domain, phylum, class_name):
            * Sea total_ids el conjunto de proteomas en esa clase.
            * Sea already_selected los ya en selected_ids (por Uabi o porque otra clase los cogió, etc.).
            * Si len(total_ids) <= min_per_class:
                  -> se añaden TODOS (para no perder info, aunque no se llegue a min_per_class).
            * Si len(total_ids)  > min_per_class:
                  -> si len(already_selected) >= min_per_class: no se añade nada.
                  -> si no, se añaden (min_per_class - len(already_selected))
                     seleccionados al azar entre los que faltan.

    Así obtienes el conjunto "mínimo" que satisface:
      - Uabi completo
      - Al menos min_per_class por clase siempre que el dataset lo permita.
    """
    rng = np.random.default_rng(random_state)
    meta = meta_df.copy()

    # 1) Todos los Uabi siempre dentro
    uabi_df = meta[meta["is_uabi"]]
    selected_ids = set(uabi_df["proteome_id"])
    print(f"[INFO] Proteomas totales: {len(meta_df)}")
    print(f"[INFO] Proteomas UABI:    {len(selected_ids)}")

    # 2) Recorremos clases
    for (dom, phyl, cls), sub in meta.groupby(["domain", "phylum", "class_name"]):
        total_ids = set(sub["proteome_id"])
        already   = total_ids & selected_ids

        n_total = len(total_ids)
        n_sel   = len(already)

        if n_total == 0:
            continue

        print(f"[CLASE] {dom} / {phyl} / {cls}  total={n_total}, selected={n_sel}")

        # Si hay <= min_per_class en todo el dataset, cogemos todos
        if n_total <= min_per_class:
            # Añadimos todo ese grupo (no podemos llegar a min_per_class, pero al menos no perdemos nadie)
            selected_ids.update(total_ids)
            print(f"   -> Clase con <= {min_per_class} proteomas; añado TODOS ({n_total})")
            continue

        # Si ya tenemos >= min_per_class en esa clase, listo
        if n_sel >= min_per_class:
            print("   -> Ya hay suficientes seleccionados para esta clase; no añado más.")
            continue

        # Si hay más de min_per_class disponibles y tenemos menos, completamos hasta min_per_class
        need = min_per_class - n_sel
        candidates = list(total_ids - selected_ids)

        if need <= 0:
            continue

        if len(candidates) <= need:
            # Añadimos todos los candidatos (es el mínimo para llegar o acercarse a min_per_class)
            to_add = candidates
        else:
            # Elegimos need al azar entre los candidatos
            to_add = list(rng.choice(candidates, size=need, replace=False))

        selected_ids.update(to_add)
        print(f"   -> Añadidos {len(to_add)} proteomas para llegar a {min_per_class}.")

    # 3) Construimos subset final
    meta_subset = meta_df[meta_df["proteome_id"].isin(selected_ids)].copy()

    # Estadísticas por clase en el subset
    class_counts_subset = (
        meta_subset.groupby(["domain", "phylum", "class_name"])
        .size()
        .reset_index(name="n_proteomes")
    )

    # Clases que siguen por debajo del mínimo (porque el dataset no tenía suficientes)
    classes_below_min = class_counts_subset[class_counts_subset["n_proteomes"] < min_per_class]

    return meta_subset, selected_ids, class_counts_subset, classes_below_min


def show_summary(meta_subset, class_counts_subset, classes_below_min, min_per_class=5):
    """
    Saca un resumen rápido en stdout.
    """
    print("\n===== RESUMEN SUBSET =====")
    print(f"Nº de proteomas en subset: {len(meta_subset)}")

    print("\n[Por dominio]")
    print(meta_subset["domain"].value_counts())

    print("\n[Top (domain, phylum, class_name) por nº de proteomas]")
    print(
        class_counts_subset
        .sort_values("n_proteomes", ascending=False)
        .head(30)
        .to_string(index=False)
    )

    if len(classes_below_min) > 0:
        print(f"\n[Clases con < {min_per_class} proteomas en el subset (no había suficientes en el dataset)]")
        print(classes_below_min.to_string(index=False))
    else:
        print(f"\nTodas las clases tienen al menos {min_per_class} proteomas en el subset.")

    print("==========================\n")

# %%


# ===========================
# EJECUCIÓN (para Spyder)
# ===========================

print(f"[INFO] Leyendo archivo: {INFILE}")
arch_df, meta_df = load_arch_and_metadata(INFILE)

print("[INFO] Ejemplo de metadata completa:")
print(meta_df.head())

meta_subset, selected_ids, class_counts_subset, classes_below_min = select_min_for_min_per_class(
    meta_df,
    min_per_class = MIN_PER_CLASS,
    random_state  = RANDOM_STATE,
)

show_summary(meta_subset, class_counts_subset, classes_below_min, min_per_class=MIN_PER_CLASS)

# Filtrar arquitecturas al subset
subset_arch = arch_df[arch_df["proteome_id"].isin(selected_ids)].copy()

print(f"[INFO] Guardando metadata subset en:      {META_OUT}")
meta_subset.to_csv(META_OUT, index=False)

print(f"[INFO] Guardando arquitecturas subset en: {ARCH_OUT}")
subset_arch.to_csv(ARCH_OUT, index=False)

print("[INFO] Hecho.")
