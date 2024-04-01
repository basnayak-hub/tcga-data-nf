library(NetSciDataCompanion)

library("optparse")

################################################
############ CLEANING FUNCTIONS ################
################################################

library(huge)
library(tidyverse)

meanImpute = function(x)
{
  x[is.na(x)] <- mean(x,na.rm=T)
  return(x)
}

# this function removes genes that weren't measured anywhere and should be done agnostic to subtype
removeMissing = function(meth_df, thres = 0.2)
{
  propMiss = apply(meth_df,2,function(x){sum(is.na(x))/length(x)})
  skipIndices = which(propMiss > thres)
  print(paste("Number of genes with > 0.2 missing:",length(skipIndices)))
  print("Genes omitted:")
  print(names(meth_df)[skipIndices])
  meth_dfClean = meth_df[,-skipIndices]
  return(meth_dfClean)
}

betaToM = function(beta)
{
  return(log2(beta/(1-beta))) # reference: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
}

# this function does imputation and transformation and should be applied within subtype
cleanMethylationData = function(meth_df, npn=T, mval=F) # meth_df is a data frame of beta means, 
# rows=samples, first column=sample ids,  cols=genes
{
  meth_dfWhole = removeMissing(meth_df,thres = 0.2)
  # skip anything w/more than 0.2 missing
  
  meth_dfComplete = meth_dfWhole
  for(i in 2:ncol(meth_dfComplete))
  {
    meth_dfComplete[,i] = meanImpute(meth_dfWhole[,i])
  }
  
  # sanity check the imputation
  summary(apply(meth_dfComplete,2,function(x){sum(is.na(x))}))
  summary(apply(meth_dfComplete[,-1],2,sd))
  
  if(mval & !npn)
  {
    transf_data = data.frame(apply(meth_dfComplete[,-1],2,betaToM))
    row.names(transf_data)= meth_dfComplete[,1]
    return(transf_data)
  }
  
  if(npn)
  {
    # do the nonparanormal transformation
    transf_data = data.frame(huge.npn(meth_dfComplete[,-1]))
    row.names(transf_data)= meth_dfComplete[,1]
    return(transf_data)
  }
  
  if(mval & npn)
  {
    # first apply M value transformation
    transf_data_m = apply(meth_dfComplete[,-1],2,betaToM)
    # next apply npn
    transf_data_m_npn = data.frame(huge.npn(transf_data_m))
    row.names(transf_data_m_npn)= meth_dfComplete[,1]
    return(transf_data_m_npn)
  }
  
  row.names(meth_dfComplete) = meth_dfComplete[,1]
  
  return(meth_dfComplete[,-1])
  
}

################################################
################################################
################################################

option_list = list(
  make_option(c("-p", "--project"), type="character", default=NULL, 
              help="project name (e.g. tcga_luad)", metavar="project_name"),
  make_option(c("-m", "--methylation"), type="character", default=NULL, 
              help="methylation data filename", metavar="methylation"),
 make_option(c("-o", "--output_table"), type="character", default='gene_level_methylation.txt', 
              help="output filename for table", metavar="output_txt"),
 make_option(c("--tissue_type"), type="character", default=" ", 
              help="tissue type", metavar="character"),
 make_option(c("--to_npn"), type="character", default="FALSE", 
              help="to_npn", metavar="character"),
 make_option(c("--to_mval"), type="character", default="TRUE", 
              help="to_mval", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

project = opt$project #args[1]
methpath = opt$methylation #args[2]
output_fn = opt$output_table #args[3]
probe_map_fn = opt$probemap #args[4]
# These: if none, keep everything
tissue_type =  opt$tissue_type #args[5]
to_npn = as.logical(opt$to_npn) #args[6]
to_mval = as.logical(opt$to_mval) #args[7]

meth_raw = read.csv(methpath,row.names=1)

# there are duplicates to handle
my_friend = NetSciDataCompanion::CreateNetSciDataCompanionObject()

# tissueDF = my_friend$getTissueType(meth_raw$TCGA_barcode[1])
# for(i in 2:nrow(meth_raw))
# {
#   tissueDF = rbind.data.frame(tissueDF,my_friend$getTissueType(meth_raw$TCGA_barcode[i]))
# }
# for tumor, use tissueType = "Primary Solid Tumor"

# primaryTumorBarcodes = tissueDF %>% 
#  dplyr::filter(description == tissueType) %>%
#  dplyr::select(TCGA_barcode)

dupes = meth_raw$TCGA_barcode[duplicated(meth_raw$TCGA_barcode)]

# filter duplicates based on name
meth_tumor = meth_raw %>% 
  dplyr::filter(!(TCGA_barcode %in% dupes)) 
  
row.names(meth_tumor) = meth_tumor$TCGA_barcode
print('Methylation data dimension:')
print(dim(meth_tumor))

# filter based on tissue type
if (to_mval){print('yes')
}
meth_tumor_clean = meth_tumor %>% dplyr::relocate(TCGA_barcode) %>%
  cleanMethylationData(npn=to_npn,mval=to_mval)

meth_tumor_clean$TCGAbarcode = row.names(meth_tumor_clean)
write.csv(meth_tumor_clean,output_fn)
