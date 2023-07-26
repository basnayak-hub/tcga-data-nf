#library(aws.s3)
#library(config)
library(data.table)
library(NetSciDataCompanion)
library(stringr)

#Sys.setenv("AWS_PROFILE" = "MFA")
args = commandArgs(trailingOnly=TRUE)

#baseDir = args[1]
project = args[1]
methpath = args[2]
#resultsDir = args[4]
outputFn = args[3]
probeMapFn = args[4]
# These: if none, keep everything
tfListFn = args[5]

#project = gsub("-","_",str_to_lower(gsub("\\[|\\]", "", project)))
print(paste("Project:",project))
print(paste("Methylation file:",methpath))

#R config = config::get(config="default", file=paste0(baseDir,"/conf/r_config/config_get_gene_map.yml"))

# meth_df = fread(paste(c(baseDir,"/ebs_work/results/",project,"_methylations.txt"),collapse=""),data.table=F)
#filepath = paste(c(baseDir,"/",resultsDir,"/",methpath),collapse="")
#print(paste("filepath",filepath))

meth_df = fread(methpath,data.table=F)
print(head(meth_df)[1:10])
probe_list = meth_df$probeID

my_friend = NetSciDataCompanion::CreateNetSciDataCompanionObject(project_name=project)

my_map = data.frame(fread(probeMapFn,sep=",",header=T),row.names=1)
# get only probes annotated to transcription factors at this point
# download.file('http://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt', 
#              destfile = "TF_names_v_1.01.txt")

# Get all gene names
#names(my_map)[2]="geneName"
print(names(my_map))
# filter methylation beta to only probes that are in the tf list
allGenes = unique(na.omit(my_map$geneNames))

if (tfListFn!=''){

# Read TF list
tfList = read.table(tfListFn)[,1]#read.table(paste(baseDir,"ext/TF_names_v_1.01.txt",sep="/"))[,1]
print(paste(c('Reading TF list', tfListFn)), collapse = " ")
print(head(tfList))

# To be removed
# if (levineExclusionsFn!=''){}
# add exclusion for Levine here
#levineExclusions = read.csv(levineExclusionsFn)[,1]#read.csv(paste(baseDir,"ext/levine_gene_exclusions.csv",sep="/"))[,1]
#print(paste('Reading Levine exclusions', levineExclusionsFn))
#print(head(levineExclusions))
#tfListFiltered = tfList[!(tfList %in% levineExclusions)]
#print("things in tfList that are not in tfListFiltered:")
#print(setdiff(tfList, tfListFiltered))

tfListFiltered = intersect(tfList,allGenes)
print(paste(c('There are', length(tfListFiltered), 'genes that are going to be used'), collapse = ' '))

geneMap = my_friend$probeToMeanPromoterMethylation(methylation_betas = meth_df, 
                                   probe_gene_map = my_map, 
                                   genesOfInterest = tfListFiltered)
} else {
  # Use all genes
  geneMap = my_friend$probeToMeanPromoterMethylation(methylation_betas = meth_df, 
                                   probe_gene_map = my_map, genesOfInterest = allGenes)
}


myTCGAbarcode = my_friend$mapUUIDtoTCGA(row.names(geneMap))
# sanity check
print("gene map row names")
print(row.names(geneMap)[1:5])
print("Sanity check on order")
sum(order(row.names(myTCGAbarcode)) != order(row.names(geneMap)))

# add TCGA and write to file 
geneMap %>% as.data.frame() %>%
  mutate(TCGA_barcode = myTCGAbarcode$submitter_id) %>%
  write.csv(outputFn,row.names=T, quote=F)

