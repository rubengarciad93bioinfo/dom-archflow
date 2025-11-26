FROM mambaorg/micromamba:1.5.6

# Usar el shell especial del base image (activa micromamba en cada RUN)
SHELL ["/usr/local/bin/_dockerfile_shell.sh"]
ARG MAMBA_DOCKERFILE_ACTIVATE=1

###############################
# 1) Perl + Moose + BioPerl
###############################
USER root
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
      perl \
      libmoose-perl \
      bioperl \
    && rm -rf /var/lib/apt/lists/*

###############################
# 2) Directorios de trabajo
###############################
RUN mkdir -p /db /work && \
    chown -R $MAMBA_USER:$MAMBA_USER /db /work

USER $MAMBA_USER
WORKDIR /work

###############################
# 3) HMMER + pfam_scan (conda)
###############################
RUN micromamba install -y -n base -c conda-forge -c bioconda \
      hmmer=3.4 \
      pfam_scan \
      requests \
      pandas \
      networkx \
      matplotlib \
      iqtree=2  \
    && micromamba clean -a -y


##############################
# 4) Dependencias extra del sistema
###############################
USER root
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
      git \
      python3 \
      python3-pip \
      parallel \
    && rm -rf /var/lib/apt/lists/*

# Clonar hpc_pfam_search
USER root
RUN git clone https://github.com/apcamargo/hpc_pfam_search.git /opt/hpc_pfam_search && \
    chmod +x /opt/hpc_pfam_search/hpc_pfam_search.pl && \
    ln -s /opt/hpc_pfam_search/hpc_pfam_search.pl /usr/local/bin/hpc_pfam_search
USER $MAMBA_USER

# Dentro de tu Dockerfile de pfam_scan_arm
# Clonar hpc_pfam_search
# RUN git clone https://github.com/apcamargo/hpc_pfam_search.git /opt/hpc_pfam_search

# (opcional) instalarlo con pip si lo necesitas como paquete Python
# RUN pip install /opt/hpc_pfam_search

###############################
# 5) Variables de entorno
###############################
ENV PERL5LIB="/opt/conda/share/pfam_scan-1.6-5:$PERL5LIB"
ENV PATH="/opt/conda/bin:${PATH}"

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["bash"]
