library("recount3")
library("NetSciDataCompanion")

args<-commandArgs(TRUE)
print(paste("Downloading data from recount3 for:",args[2],'project:',args[1]))
print(paste("Saving to:",args[6]))

# args[1] = project for gtex
# args[2] = project for tcga
# args[3] = organism
# args[4] = gencode version
# args[5] = type

project = args[1]
project_home = args[2]
organism = args[3]
gencodev = args[4]
recount_type = args[5]
sample_list = args[6]
output_fn = args[7]

print(args)
print('Downloading...')
## Download data and create the RSE object
recount_data = recount3::create_rse_manual(
  project = project,
  project_home = project_home, #"data_sources/gtex",
  organism = organism,
  annotation = gencodev,
  type = recount_type
)

print('Data have been downloaded...')
print('Clinical Data has shape (samples X clinical variables):')
print(dim(recount_data@colData))
print('Expression Data has shape (genes X samples):')
print(dim(recount_data))

print(nchar(sample_list))
print(sample_list)

# Adding support to filter the sample list
if (nchar(sample_list)>3){
  print('Filtering only samples in the sample list...')
  ong = NetSciDataCompanion::CreateNetSciDataCompanionObject()
  submitters = read.table(sample_list, header = FALSE, sep = "", dec = ".")
  recount_data = recount_data[,ong$extractSampleAndTypeAndVial(recount_data$tcga.tcga_barcode) %in% submitters$V1]

  print('Data have been filtered by sample list...')
  print('Clinical Data has shape (samples X clinical variables):')
  print(dim(recount_data@colData))
  print('Expression Data has shape (genes X samples):')
  print(dim(recount_data))

}

print('Saving...')
saveRDS(recount_data, output_fn)
print('DONE!')