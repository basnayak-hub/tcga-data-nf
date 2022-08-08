library("recount3")

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
output_fn = args[6]

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
print('Saving...')
saveRDS(recount_data, output_fn)
print('DONE!')