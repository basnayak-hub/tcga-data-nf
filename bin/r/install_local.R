library(devtools)


print('Installing')


devtools::install_github('pmandros/TCGAPurityFiltering')
devtools::install_github('immunogenomics/presto')
devtools::install_github('aet21/EpiSCORE')
print(getwd())
devtools::install_local('/home/ubuntu/viola/tcga-data-nf/NetSciDataCompanion/', force=TRUE, dependencies=T)