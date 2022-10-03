library(TCGAbiolinks)
library(tidyverse)
args<-commandArgs(TRUE)

# args[1] = project for TCGA (e.g., TCGA-BRCA)
# args[2] = data category, should be "Clinical"
# args[3] = data type, should be "Clinical Supplement"
# args[4] = data format, should be "BCR Biotab"
# args[5] = directory, where to download the data

print("I am your assistant for today. I will be 
      downloading the clinical data for you.")

project = args[1]
data_category = args[2]
data_type = args[3]
data_format = args[4]
directory = args[5]

###some asserts
stopifnot(grepl("TCGA", args[1]))

print('Printing parameters...')
print(args)

print('Downloading the data...')
query <- GDCquery(project = project, 
                  data.category = data_category,
                  data.type = data_type, 
                  data.format = data_format)
GDCdownload(query)
clinical.BCRtab.all <- GDCprepare(query)
print('Data most probably have been downloaded...')

print('Names of clinical files: ')
print(names(clinical.BCRtab.all))

print('Saving the data...')
for(clinical_file in names(clinical.BCRtab.all)){
  # remove escapes
  data <- data.frame(lapply(clinical.BCRtab.all[[clinical_file]], function(x) {gsub("'", " ", x)}))
  
  write_csv(data,
            paste(paste(directory,clinical_file,sep = .Platform$file.sep),'.csv', sep=''))
}
print('Data most probably have been saved...')

print('Goodbye and I hope you had fun.')

