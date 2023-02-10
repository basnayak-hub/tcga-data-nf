## base container
FROM mambaorg/micromamba:0.15.3 as base_container
USER root
RUN apt-get update \
    && apt-get install --yes rename tini procps curl \
    && apt-get autoremove \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

USER micromamba
ENTRYPOINT ["tini", "--"]
CMD ["/bin/bash"]

## main container
ARG CONDA_FILE=containers/env.base.python.yml
FROM base_container
# adding opencontainer labels to link registry to github repository
LABEL org.opencontainers.image.title="tcga-data-nf"
LABEL org.opencontainers.image.description="Workflow to download and prepare TCGA data"
LABEL org.opencontainers.image.url="https://github.com/QuackenbushLab/tcga-data-nf"
LABEL org.opencontainers.image.documentation="https://github.com/QuackenbushLab/tcga-data-nf"
LABEL org.opencontainers.image.source="https://github.com/QuackenbushLab/tcga-data-nf"
LABEL org.opencontainers.image.vendor="QuackenbushLab"
LABEL org.opencontainers.image.authors="Viola Fanfani"
LABEL org.opencontainers.image.revision="v0.0.9"

COPY --chown=micromamba:micromamba ${CONDA_FILE} /tmp
RUN micromamba install -y -n base -f /tmp/`basename ${CONDA_FILE}` && \
    micromamba clean --all --yes


#FROM rocker/tidyverse:latest
# Install R packages


# Docker inheritance
FROM bioconductor/bioconductor_docker:devel

# copy install scripts
COPY containers/*.sh /opt 

RUN mkdir -p /opt/install_scripts && \ 
    mv /opt/install*.sh /opt/install_scripts && \
    chmod +x /opt/install_scripts/install*.sh

# install tidyverse
RUN /opt/install_scripts/install_tidyverse.sh

###
RUN apt-get update
## Install cran packages
RUN R -e "install.packages(c('devtools','optparse','ggplot2','cowplot'),\
                            dependencies=TRUE, \
                           repos='http://cran.rstudio.com/')"


## Install our own package

#COPY /Users/violafanfani/Documents/uni-harvard/repos/NetSciDataCompanion /opt/NetSciDataCompanion 
#RUN R -e "devtools::install_local('/opt/NetSciDataCompanion', force=TRUE)"
## Install bioconductior packages                           
RUN R -e 'BiocManager::install(ask = F)' && R -e 'BiocManager::install(c("recount3", \
    "edgeR", "GenomicAlignments", "Biostrings", "SummarizedExperiment", "GenomicFiles","recount",\
    "TCGAbiolinks", "GenomicDataCommons","limma",ask = F))'


ADD  NetSciDataCompanion/ /opt/NetSciDataCompanion/
RUN R -e "devtools::install_github('pmandros/TCGAPurityFiltering')"
RUN R -e "devtools::install_local('/opt/NetSciDataCompanion/', force=TRUE, dependencies=T)"
