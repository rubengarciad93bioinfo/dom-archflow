#!/usr/bin/env python3
from pathlib import Path
import argparse
import re
import sys
import json
from typing import Optional, Tuple

try:
    import requests
except ImportError:
    requests = None


def dprint(msg: str, debug: bool):
    """Imprime mensajes de debug a stderr si debug=True."""
    if debug:
        print(f"[DEBUG] {msg}", file=sys.stderr)


def extract_assembly_from_gff(gff: Path, debug: bool = False) -> str:
    """
    Extrae un accession tipo GCF_... del GFF.

    Cubre:
    - NCBI:  ##genome-build-accession NCBI_Assembly:GCF_000009045.1
    - SGD:   #!refseq-version GCF_000146045.2
    """
    if not gff.exists():
        dprint(f"GFF no existe: {gff}", debug)
        return "NA"

    assembly = "NA"
    with gff.open("rt") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("##genome-build-accession") or line.startswith("#!genome-build-accession"):
                parts = line.split()
                if len(parts) >= 2 and "GCF_" in parts[1]:
                    assembly = parts[1].split(":")[-1]
                    break
            if line.startswith("#!refseq-version"):
                parts = line.split()
                if len(parts) >= 2 and parts[1].startswith("GCF_"):
                    assembly = parts[1]
                    break

    dprint(f"Assembly extraído de {gff}: {assembly}", debug)
    return assembly


def extract_taxid_from_gff(gff: Path, debug: bool = False) -> str:
    """
    Intenta extraer TaxID de una línea ##species ... id=XXXX (formato NCBI).
    """
    if not gff.exists():
        dprint(f"GFF no existe: {gff}", debug)
        return "NA"

    taxid_re = re.compile(r"id=(\d+)")
    tax_id = "NA"

    with gff.open("rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if line.startswith("##species"):
                m = taxid_re.search(line)
                if m:
                    tax_id = m.group(1)
                    break

    dprint(f"TaxID extraído de {gff}: {tax_id}", debug)
    return tax_id


def extract_organism_from_faa(faa: Path, debug: bool = False) -> str:
    """
    Coge el primer header FASTA y, si puede, lo que haya entre [ ].
    Si no, devuelve la descripción completa (menos '>').
    """
    if not faa.exists():
        dprint(f"FAA no existe: {faa}", debug)
        return "NA"

    organism = "NA"
    with faa.open("rt") as fh:
        for line in fh:
            if line.startswith(">"):
                header = line[1:].strip()
                m = re.search(r"\[([^\]]+)\]", header)
                if m:
                    organism = m.group(1)
                else:
                    parts = header.split(None, 1)
                    organism = parts[1] if len(parts) > 1 else parts[0]
                break

    dprint(f"Organism extraído de {faa}: {organism}", debug)
    return organism


def fetch_from_taxid(
    tax_id: str,
    api_key: Optional[str] = None,
    debug: bool = False,
) -> Tuple[str, str]:
    """
    Usa NCBI E-utilities EFetch (db=taxonomy, XML) para obtener:
      - Lineage completo (<Lineage>)
      - ScientificName como backup.

    Devuelve (tax_id, taxonomy_string).
    """
    if not requests:
        dprint("requests no está disponible, no se puede llamar a NCBI Taxonomy", debug)
        return tax_id, "NA"
    if tax_id in ("NA", "", None):
        dprint("TaxID es NA, no se llama a EFetch", debug)
        return "NA", "NA"

    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "taxonomy",
        "id": tax_id,
        "retmode": "xml",
        "tool": "dom-archflow",
        "email": "tu_email@ejemplo.org",
    }
    if api_key:
        params["api_key"] = api_key

    dprint(f"Llamando a EFetch Taxonomy con tax_id={tax_id}", debug)

    try:
        resp = requests.get(base, params=params, timeout=10)
        dprint(f"[EFETCH] URL: {resp.url}", debug)
        dprint(f"[EFETCH] Status: {resp.status_code}", debug)
        resp.raise_for_status()

        text_snip = resp.text[:400].replace("\n", " ")
        dprint(f"[EFETCH] Respuesta (snippet 400 chars): {text_snip}", debug)

        import xml.etree.ElementTree as ET
        root = ET.fromstring(resp.text)
        taxon = root.find(".//Taxon")
        if taxon is None:
            dprint("[EFETCH] No se encontró nodo <Taxon>", debug)
            return tax_id, "NA"

        sci_name = taxon.findtext("ScientificName") or "NA"
        lineage = taxon.findtext("Lineage") or ""
        lineage = lineage.strip()

        dprint(f"[EFETCH] ScientificName: {sci_name}", debug)
        dprint(f"[EFETCH] Lineage: {lineage}", debug)

        taxonomy = lineage if lineage else sci_name
        return str(tax_id), taxonomy

    except Exception as e:
        dprint(f"[EFETCH] Error: {e}", debug)
        return tax_id, "NA"


def fetch_from_assembly(
    assembly: str,
    api_key: Optional[str] = None,
    debug: bool = False,
) -> Tuple[str, str]:
    """
    Usa NCBI Datasets Genome para obtener tax_id y nombre científico
    a partir de un assembly accession (GCF_...).

    Devuelve (tax_id, sci_name).
    """
    if not requests:
        dprint("requests no está disponible, no se puede llamar a NCBI Datasets", debug)
        return "NA", "NA"
    if assembly in ("NA", "", None):
        dprint("Assembly es NA, no se llama a Datasets", debug)
        return "NA", "NA"

    base = "https://api.ncbi.nlm.nih.gov/datasets/v2"
    url = f"{base}/genome/accession/{assembly}"
    headers = {"accept": "application/json"}
    if api_key:
        headers["X-Api-Key"] = api_key

    dprint(f"Llamando a Datasets Genome con assembly={assembly}", debug)

    try:
        resp = requests.get(url, headers=headers, timeout=10)
        dprint(f"[DATASETS] URL: {url}", debug)
        dprint(f"[DATASETS] Status: {resp.status_code}", debug)
        resp.raise_for_status()

        try:
            data = resp.json()
            dprint(f"[DATASETS] JSON keys top-level: {list(data.keys())}", debug)
        except Exception as e:
            dprint(f"[DATASETS] Error parseando JSON: {e}", debug)
            return "NA", "NA"

        reports = data.get("reports") or data.get("genome_reports") or []
        if isinstance(reports, list) and reports:
            rep = reports[0]
        elif isinstance(reports, dict):
            rep = reports
        else:
            rep = data

        tax_node = rep.get("taxonomy", {})
        dprint(f"[DATASETS] taxonomy node: {json.dumps(tax_node)[:400]}", debug)

        tax_id = tax_node.get("tax_id") or rep.get("tax_id") or "NA"
        tax_id = str(tax_id)

        sci_name = (
            tax_node.get("organism_name")
            or rep.get("organism", {}).get("organism_name")
            or rep.get("organism_name")
            or "NA"
        )

        dprint(f"[DATASETS] tax_id: {tax_id}", debug)
        dprint(f"[DATASETS] organism_name: {sci_name}", debug)

        return tax_id, sci_name

    except Exception as e:
        dprint(f"[DATASETS] Error: {e}", debug)
        return "NA", "NA"


def main():
    ap = argparse.ArgumentParser(
        description="Construye metadata para un proteoma (assembly, tax_id, taxonomy, completeness)"
    )
    ap.add_argument("--proteome-id", required=True)
    ap.add_argument("--faa", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ncbi-api-key", default=None)
    ap.add_argument(
        "--offline-only",
        action="store_true",
        help="No usar NCBI; solo lo que se pueda inferir de GFF/FAA",
    )
    ap.add_argument(
        "--debug",
        action="store_true",
        help="Imprimir información detallada de NCBI y parsing",
    )
    args = ap.parse_args()

    faa = Path(args.faa)
    gff = Path(args.gff)

    dprint(f"=== build_metadata.py DEBUG ===", args.debug)
    dprint(f"proteome_id = {args.proteome_id}", args.debug)
    dprint(f"FAA = {faa}", args.debug)
    dprint(f"GFF = {gff}", args.debug)

    assembly = extract_assembly_from_gff(gff, debug=args.debug)
    tax_id = extract_taxid_from_gff(gff, debug=args.debug)
    organism = extract_organism_from_faa(faa, debug=args.debug)

    completeness = "NA"

    taxonomy: Optional[str] = None
    sci_from_asm: Optional[str] = None

    # 0) Si no es offline y tenemos assembly pero NO tax_id, intenta Datasets
    if (
        not args.offline_only
        and requests is not None
        and tax_id in ("NA", "", None)
        and assembly not in ("NA", "", None)
    ):
        dprint("No hay tax_id, pero sí assembly -> llamar a Datasets primero", args.debug)
        tid2, sci2 = fetch_from_assembly(assembly, api_key=args.ncbi_api_key, debug=args.debug)
        if tid2 not in ("NA", "", None):
            tax_id = tid2
        if sci2 not in ("NA", "", None):
            sci_from_asm = sci2

    # 1) Intentar SIEMPRE por TaxID vía EFetch (Taxonomy)
    if (
        not args.offline_only
        and requests is not None
        and tax_id not in ("NA", "", None)
    ):
        dprint(f"Llamando a EFetch con tax_id final={tax_id}", args.debug)
        tid, tx = fetch_from_taxid(tax_id, api_key=args.ncbi_api_key, debug=args.debug)
        tax_id = tid
        if tx not in ("NA", "", None):
            taxonomy = tx

    # 2) Si Taxonomy no ha dado nada, usar sci_name de Datasets si lo tenemos
    if taxonomy is None and sci_from_asm not in (None, "", "NA"):
        dprint("Usando sci_name de Datasets como taxonomy", args.debug)
        taxonomy = sci_from_asm

    # 3) Si aún no hay taxonomy, usar organism del FASTA
    if taxonomy is None and organism not in ("NA", "", None):
        dprint("Usando organism del FASTA como taxonomy", args.debug)
        taxonomy = organism

    # 4) Si todavía no hay taxonomy, pero sí tax_id o assembly, poner marcadores
    if taxonomy is None:
        if tax_id not in ("NA", "", None):
            taxonomy = f"NCBI_taxid:{tax_id}"
        elif assembly not in ("NA", "", None):
            taxonomy = f"NCBI_Assembly:{assembly}"
        else:
            taxonomy = "Undetermined_taxon"

    safe_taxonomy = (taxonomy or "").replace('"', "")

    dprint(f"== RESUMEN FINAL ==", args.debug)
    dprint(f"assembly   = {assembly}", args.debug)
    dprint(f"tax_id     = {tax_id}", args.debug)
    dprint(f"taxonomy   = {safe_taxonomy}", args.debug)
    dprint(f"organism   = {organism}", args.debug)

    out_path = Path(args.out)
    with out_path.open("w") as f:
        f.write("proteome_id,proteins_faa,genes_gff,assembly,tax_id,taxonomy,completeness\n")
        f.write(
            f"{args.proteome_id},{faa.name},{gff.name},"
            f"{assembly},{tax_id},{safe_taxonomy},{completeness}\n"
        )


if __name__ == "__main__":
    main()
