#!/usr/bin/env python3
from pathlib import Path
import argparse
from typing import Optional


def calcular_solapamiento(inicio1: int, fin1: int, inicio2: int, fin2: int) -> float:
    solapamiento = max(0, min(fin1, fin2) - max(inicio1, inicio2))
    longitud_max = max(fin1 - inicio1, fin2 - inicio2)
    return (solapamiento / longitud_max) * 100 if longitud_max > 0 else 0.0


def ajustar_solapamientos(data, umbral_solapamiento: Optional[int] = None):
    data.sort(key=lambda x: (x["seqId"], x["inicio"]))
    resultados_filtrados = []
    arquitectura_name, arquitectura_dominios, arquitectura_clanes, arquitectura_tipos = {}, {}, {}, {}

    for item in data:
        seqId = item["seqId"]
        item["clan"] = item["clan"] if item["clan"] not in ("No_clan", "Sin_clan") else item["hmm"]

        if seqId not in arquitectura_dominios:
            resultados_filtrados.append(item)
            arquitectura_dominios[seqId] = [item["hmm"]]
            arquitectura_clanes[seqId] = [item["clan"]]
            arquitectura_tipos[seqId] = [item["tipo"]]
            arquitectura_name[seqId] = [item["hmm_name"]]
        else:
            ultimo_item = resultados_filtrados[-1]
            solapamiento = calcular_solapamiento(
                ultimo_item["inicio"], ultimo_item["fin"], item["inicio"], item["fin"]
            )
            if umbral_solapamiento is None or solapamiento <= umbral_solapamiento:
                resultados_filtrados.append(item)
                arquitectura_dominios[seqId].append(item["hmm"])
                arquitectura_clanes[seqId].append(item["clan"])
                arquitectura_tipos[seqId].append(item["tipo"])
                arquitectura_name[seqId].append(item["hmm_name"])
            elif item["bit_score"] > ultimo_item["bit_score"]:
                resultados_filtrados[-1] = item
                arquitectura_dominios[seqId][-1] = item["hmm"]
                arquitectura_clanes[seqId][-1] = item["clan"]
                arquitectura_tipos[seqId][-1] = item["tipo"]
                arquitectura_name[seqId][-1] = item["hmm_name"]

    return resultados_filtrados, arquitectura_dominios, arquitectura_clanes, arquitectura_tipos, arquitectura_name


def analizar_dominios_un_archivo(
    pfam_file: Path,
    output_file: Path,
    proteome_id: str,
    version: str,
    umbral_solapamiento: int,
    tax_id: str,
    taxonomy: str,
    completeness: str,
):
    pfam_file = pfam_file.resolve()
    output_file = output_file.resolve()

    dominios = []
    with pfam_file.open("rt") as data:
        for linea in data:
            if linea.startswith("#") or linea.strip() == "":
                continue
            dominios.append(linea)

    data2 = []
    for linea in dominios:
        elementos = linea.strip().split()
        try:
            evalue = float(elementos[12])
            if evalue < 0.001:
                data2.append(
                    {
                        "seqId": elementos[0],
                        "inicio": int(elementos[1]),
                        "fin": int(elementos[2]),
                        "bit_score": float(elementos[11]),
                        "i_evalue": evalue,
                        "hmm": elementos[5].split(".")[0],
                        "hmm_name": elementos[6],
                        "tipo": elementos[7],
                        "clan": elementos[14] if len(elementos) > 14 else "Sin_clan",
                    }
                )
        except Exception:
            continue

    if not data2:
        print(f"⚠️ No se encontraron dominios en {pfam_file}")
        return

    _, arq_dominios, arq_clanes, arq_tipos, arq_name = ajustar_solapamientos(
        data2, umbral_solapamiento
    )

    with output_file.open("w") as f:
        f.write(
            "seqId,proteome_id,tax_id,taxonomy,doms,names,clans,types,completeness,version\n"
        )
        for seqId, doms in arq_dominios.items():
            doms_str = "-".join(doms)
            tipos_str = "-".join(arq_tipos[seqId])
            clanes_str = "-".join(arq_clanes[seqId])
            names_str = "-".join(arq_name[seqId])

            f.write(
                f"{seqId},{proteome_id},{tax_id},{taxonomy},"
                f"{doms_str},{names_str},{clanes_str},{tipos_str},"
                f"{completeness},{version}\n"
            )


def main():
    ap = argparse.ArgumentParser(
        description="Construye arquitecturas de dominios a partir de salida PfamScan"
    )
    ap.add_argument("--pfam", required=True)
    ap.add_argument("--proteome-id", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--version", default="v1")
    ap.add_argument("--overlap", type=int, default=30)

    # AHORA: metadata viene ya calculada
    ap.add_argument("--tax-id", default="NA")
    ap.add_argument("--taxonomy", default="Undetermined_taxon")
    ap.add_argument("--completeness", default="NA")

    args = ap.parse_args()

    analizar_dominios_un_archivo(
        pfam_file=Path(args.pfam),
        output_file=Path(args.out),
        proteome_id=args.proteome_id,
        version=args.version,
        umbral_solapamiento=args.overlap,
        tax_id=args.tax_id,
        taxonomy=args.taxonomy,
        completeness=args.completeness,
    )


if __name__ == "__main__":
    main()
