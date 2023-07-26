library(NetSciDataCompanion)

args = commandArgs(trailingOnly=TRUE)

baseDir = args[1]
project = args[2]
outputFn = args[3]
methPath = args[4]
tissueType = args[5]
to_npn = as.logical(args[6])
to_mval = as.logical(args[7])

source(paste0(baseDir,"/bin/r/methylationCleaningFunctions.R"))

# read methylation data from file
#print(paste(c(baseDir,"/",resultsDir,methPath),collapse=""))
meth_raw = read.csv(methPath,row.names=1)
print(names(meth_raw)[1:5])

# there are duplicates to handle

my_friend = NetSciDataCompanion::CreateNetSciDataCompanionObject()

# tissueDF = my_friend$getTissueType(meth_raw$TCGA_barcode[1])
# for(i in 2:nrow(meth_raw))
# {
#   tissueDF = rbind.data.frame(tissueDF,my_friend$getTissueType(meth_raw$TCGA_barcode[i]))
# }                   

# for tumor, use tissueType = "Primary Solid Tumor"

#primaryTumorBarcodes = tissueDF %>% 
#  dplyr::filter(description == tissueType) %>%
#  dplyr::select(TCGA_barcode)

dupes = meth_raw$TCGA_barcode[duplicated(meth_raw$TCGA_barcode)]

# filter duplicates based on name
meth_tumor = meth_raw %>% 
  dplyr::filter(!(TCGA_barcode %in% dupes)) 
  
row.names(meth_tumor) = meth_tumor$TCGA_barcode
print('Methylation data dimension:')
print(dim(meth_tumor))

print(to_npn)
print(to_mval)
print(class(to_npn))
if (to_mval){print('yes')}
meth_tumor_clean = meth_tumor %>% dplyr::relocate(TCGA_barcode) %>%
  cleanMethylationData(npn=to_npn,mval=to_mval)

meth_tumor_clean$TCGAbarcode = row.names(meth_tumor_clean)
write.csv(meth_tumor_clean,outputFn)
