FROM debian:stable-slim

RUN apt-get update && apt-get install -y \
    hmmer \
    perl \
    wget \
    tar \
    gzip \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Instalar PfamScan
RUN wget https://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/pfam_scan.tar.gz \
 && tar -xzf pfam_scan.tar.gz \
 && mv pfam_scan /opt/pfam \
 && chmod +x /opt/pfam/pfam_scan.pl

ENV PATH="/opt/pfam:${PATH}"

WORKDIR /work
