FROM violafanfani/tcga-data-base-r:0.0.5

USER root

# Combine apt-get update and install dependencies
RUN apt-get update && apt-get install -y \
    iputils-ping \
    ca-certificates \
    python3 \
    python3-pip \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev

RUN echo "alias python=python3" >> ~/.bashrc
RUN echo "alias python=python3" >> ~/.profile

# Upgrade pip and install Python dependencies
RUN python3 -m pip install --upgrade pip
COPY containers/requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt
RUN pip install --no-cache-dir git+https://github.com/netZoo/netZooPy.git

# Install R dependencies
RUN R -e "devtools::install_github('pmandros/TCGAPurityFiltering')"
RUN R -e "devtools::install_github('immunogenomics/presto')"
RUN R -e "devtools::install_github('aet21/EpiSCORE')"
RUN R -e "devtools::install_github('https://github.com/QuackenbushLab/NetworkDataCompanion')"
RUN R -r "BiocManager::install('GENIE3')" 

# Final clean-up if necessary
# RUN apt-get clean && rm -rf /var/lib/apt/lists/*