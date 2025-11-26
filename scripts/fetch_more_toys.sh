#!/usr/bin/env bash
set -euo pipefail

have() { command -v "$1" >/dev/null 2>&1; }

if ! have curl; then
  echo "ERROR: se requiere 'curl' instalado." >&2
  exit 1
fi

mkdir -p assets/toys

########################################
# 1) Cyanobacteria: Synechocystis sp. PCC 6803
#    RefSeq: GCF_000009725.1_ASM972v1
########################################
CYANO_DIR="assets/toys/Synechocystis"
mkdir -p "$CYANO_DIR"

CYANO_BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/725/GCF_000009725.1_ASM972v1"

echo "[INFO] Descargando Synechocystis sp. PCC 6803..."
curl -L "$CYANO_BASE/GCF_000009725.1_ASM972v1_protein.faa.gz" \
  -o "$CYANO_DIR/proteins.faa.gz"
curl -L "$CYANO_BASE/GCF_000009725.1_ASM972v1_genomic.gff.gz" \
  -o "$CYANO_DIR/genes.gff.gz"

gunzip -f "$CYANO_DIR/proteins.faa.gz"
gunzip -f "$CYANO_DIR/genes.gff.gz"


########################################
# 2) Planctomycete: Planctopirus limnophila DSM 3776
#    👉 TÚ rellenas el accession exacto y el nombre del assembly.
#    Pasos:
#      - Ve a NCBI Assembly, busca "Planctopirus limnophila DSM 3776 RefSeq"
#      - Copia el accession tipo GCF_xxx... y el "Assembly name", ej. ASM2592v1
#      - Construye la ruta 'all' como en el ejemplo de Synechocystis.
########################################
PLANC_DIR="assets/toys/Planctopirus"
mkdir -p "$PLANC_DIR"

# TODO: RELLENAR ESTAS DOS VARIABLES CON LO QUE SAQUES DE NCBI
PLANC_ACC="GCF_000092105.1"          # ej. GCF_000025925.1
PLANC_ASM="ASM9210v1"                 # ej. ASM2592v1

# La ruta 'all' se construye partiendo del accession:
#   GCF_000025925.1 -> GCF/000/025/925/GCF_000025925.1_ASM2592v1
# Rellena el bloque de 3x3 dígitos acorde al accession.
PLANC_BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/092/105/${PLANC_ACC}_${PLANC_ASM}"

echo "[INFO] Descargando Planctopirus limnophila (Planctomycetota)..."
curl -L "${PLANC_BASE}/${PLANC_ACC}_${PLANC_ASM}_protein.faa.gz" \
  -o "$PLANC_DIR/proteins.faa.gz"
curl -L "${PLANC_BASE}/${PLANC_ACC}_${PLANC_ASM}_genomic.gff.gz" \
  -o "$PLANC_DIR/genes.gff.gz"

gunzip -f "$PLANC_DIR/proteins.faa.gz"
gunzip -f "$PLANC_DIR/genes.gff.gz"


########################################
# 3) Asgardarchaeota: p.ej. Candidatus Lokiarchaeota archaeon GC14_75
#    👉 Igual que arriba: rellena accession + assembly name.
########################################
ASGARD_DIR="assets/toys/Lokiarchaeota"
mkdir -p "$ASGARD_DIR"

# TODO: RELLENAR ESTAS DOS VARIABLES
ASGARD_ACC="GCA_000986845.1"   # accession RefSeq para ese genoma
ASGARD_ASM="ASM98684v1"         # assembly name

ASGARD_BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/986/845/${ASGARD_ACC}_${ASGARD_ASM}"

echo "[INFO] Descargando Asgardarchaeota (Lokiarchaeota)..."
curl -L "${ASGARD_BASE}/${ASGARD_ACC}_${ASGARD_ASM}_protein.faa.gz" \
  -o "$ASGARD_DIR/proteins.faa.gz"
curl -L "${ASGARD_BASE}/${ASGARD_ACC}_${ASGARD_ASM}_genomic.gff.gz" \
  -o "$ASGARD_DIR/genes.gff.gz"

gunzip -f "$ASGARD_DIR/proteins.faa.gz"
gunzip -f "$ASGARD_DIR/genes.gff.gz"


########################################
# 4) Eucariota "grande": p.ej. Caenorhabditis elegans
#    (O el que prefieras, pero C. elegans es bastante estándar)
#    👉 De nuevo, rellenas accession + assembly.
#    OJO: esto es mucho más grande que E. coli / B. subtilis.
########################################
CELEG_DIR="assets/toys/Celegans"
mkdir -p "$CELEG_DIR"

# TODO: RELLENAR ESTAS DOS VARIABLES
CELEG_ACC="GCF_000002985.6"   # este es probable, pero verifícalo en NCBI
CELEG_ASM="WBcel235"          # assembly name típico de C. elegans

# Para WBcel235, la ruta en 'all' es algo tipo:
# GCF/000/002/985/GCF_000002985.6_WBcel235
CELEG_BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/${CELEG_ACC}_${CELEG_ASM}"

echo "[INFO] Descargando Caenorhabditis elegans..."
curl -L "${CELEG_BASE}/${CELEG_ACC}_${CELEG_ASM}_protein.faa.gz" \
  -o "$CELEG_DIR/proteins.faa.gz"
curl -L "${CELEG_BASE}/${CELEG_ACC}_${CELEG_ASM}_genomic.gff.gz" \
  -o "$CELEG_DIR/genes.gff.gz"

gunzip -f "$CELEG_DIR/proteins.faa.gz"
gunzip -f "$CELEG_DIR/genes.gff.gz"


########################################
# 5) Validación rápida (igual que en tu script original)
########################################
for org in Ecoli Bsubtilis Mjannaschii Scerevisiae \
           Synechocystis Planctopirus Lokiarchaeota Celegans; do
  f="assets/toys/${org}/proteins.faa"
  if [ -f "$f" ]; then
    grep -q '^>' "$f" || { echo "ERROR: $f no es FASTA válido"; exit 1; }
  fi
done

for g in assets/toys/*/genes.gff; do
  [ -f "$g" ] || continue
  grep -q -E '^\w+\t\w+\tgene\t' "$g" || {
    echo "ADVERTENCIA: $g no parece GFF con features 'gene' estándar"; 
  }
done

echo "[OK] Descargas y validaciones (hasta donde están rellenas) completas."
