#!/usr/bin/env bash
set -euo pipefail

# Requisitos mínimos: curl y gzip (o zcat/gunzip)
have() { command -v "$1" >/dev/null 2>&1; }

if ! have curl; then
  echo "ERROR: se requiere 'curl' instalado." >&2
  exit 1
fi

mkdir -p assets/toys

# ---------- E. coli K-12 MG1655 (RefSeq) ----------
# GCF_000005845.2 ASM584v2
E_DIR="assets/toys/Ecoli"
mkdir -p "$E_DIR"
E_BASE="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Escherichia_coli/reference/GCF_000005845.2_ASM584v2"

curl -L "$E_BASE/GCF_000005845.2_ASM584v2_protein.faa.gz" -o "$E_DIR/proteins.faa.gz"
curl -L "$E_BASE/GCF_000005845.2_ASM584v2_genomic.gff.gz"  -o "$E_DIR/genes.gff.gz"
gunzip -f "$E_DIR/proteins.faa.gz"
gunzip -f "$E_DIR/genes.gff.gz"

# ---------- Bacillus subtilis 168 (RefSeq) ----------
# GCF_000009045.1 ASM904v1
B_DIR="assets/toys/Bsubtilis"
mkdir -p "$B_DIR"
B_BASE="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Bacillus_subtilis/reference/GCF_000009045.1_ASM904v1"

curl -L "$B_BASE/GCF_000009045.1_ASM904v1_protein.faa.gz" -o "$B_DIR/proteins.faa.gz"
curl -L "$B_BASE/GCF_000009045.1_ASM904v1_genomic.gff.gz"  -o "$B_DIR/genes.gff.gz"
gunzip -f "$B_DIR/proteins.faa.gz"
gunzip -f "$B_DIR/genes.gff.gz"

# ---------- Methanocaldococcus jannaschii (RefSeq) ----------
# GCF_000091665.1 ASM9166v1
M_DIR="assets/toys/Mjannaschii"
mkdir -p "$M_DIR"
M_BASE="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/Methanocaldococcus_jannaschii/reference/GCF_000091665.1_ASM9166v1"

curl -L "$M_BASE/GCF_000091665.1_ASM9166v1_protein.faa.gz" -o "$M_DIR/proteins.faa.gz"
curl -L "$M_BASE/GCF_000091665.1_ASM9166v1_genomic.gff.gz"  -o "$M_DIR/genes.gff.gz"
gunzip -f "$M_DIR/proteins.faa.gz"
gunzip -f "$M_DIR/genes.gff.gz"

# ---------- Saccharomyces cerevisiae S288C (SGD) ----------

S_DIR="assets/toys/Scerevisiae"
mkdir -p "$S_DIR"

# OJO: esta ruta es un ejemplo de patrón, tienes que mirar en NCBI
# cuál es exactamente el subdirectorio (nombre completo del asm)
S_BASE="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Saccharomyces_cerevisiae/reference/GCF_000146045.2_R64"

curl -L "$S_BASE/GCF_000146045.2_R64_protein.faa.gz" -o "$S_DIR/proteins.faa.gz"
curl -L "$S_BASE/GCF_000146045.2_R64_genomic.gff.gz"  -o "$S_DIR/genes.gff.gz"

gunzip -f "$S_DIR/proteins.faa.gz"
gunzip -f "$S_DIR/genes.gff.gz"

# ---------- Validación mínima (tiene que haber '>' en los FASTA) ----------
for org in Ecoli Bsubtilis Mjannaschii Scerevisiae; do
  f="assets/toys/${org}/proteins.faa"
  grep -q '^>' "$f" || { echo "ERROR: $f no es FASTA válido"; exit 1; }
done
# opcional: validar GFF si existe
for g in assets/toys/*/genes.gff; do
  [ -f "$g" ] || continue
  grep -q -E '^\w+\t\w+\tgene\t' "$g" || { echo "ERROR: $g no parece GFF gene features"; exit 1; }
done
echo "[OK] Descargas y validaciones completas."

