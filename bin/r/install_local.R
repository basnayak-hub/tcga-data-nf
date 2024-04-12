library(devtools)
devtools::install_github('pmandros/TCGAPurityFiltering')
devtools::install_github('immunogenomics/presto')
devtools::install_github('aet21/EpiSCORE')

ADD  NetSciDataCompanion/ /opt/NetSciDataCompanion/
RUN R -e "devtools::install_local('/opt/NetSciDataCompanion/', force=TRUE, dependencies=T)"