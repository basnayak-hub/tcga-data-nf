library(data.table)
library(stringr)
library("recount3")
library(NetworkDataCompanion)
library("optparse")
library(dplyr)

option_list = list(
  make_option(c("-p", "--project"), type="character", default=NULL, 
              help="project name (e.g. tcga_luad)", metavar="project_name"),
  make_option(c("-m", "--methylation"), type="character", default=NULL, 
              help="methylation data filename", metavar="methylation"),
 make_option(c("-o", "--output_table"), type="character", default='gene_level_methylation.txt', 
              help="output filename for table", metavar="output_txt"),
make_option(c("--probemap"), type="character", default = 'assets/450k_promoter_probe_map_TSS200_TSS0_one_probe_to_many_genes.csv',
              help="probe map file", metavar="probemap"),
 make_option(c("--tf_list"), type="character", default=" ", 
              help="TF list filename. Pass a text file with TF names to filter the output", metavar="character"),
  make_option("--filter_duplicates_missingness", action = "store_true", default = FALSE,
              help = "Filter duplicate or missing values [default: %default]")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

project = opt$project #args[1]
methpath = opt$methylation #args[2]
output_fn = opt$output_table #args[3]
probe_map_fn = opt$probemap #args[4]
# These: if none, keep everything
tf_list_fn =  opt$tf_list #args[5]
filter_duplicates_missingness = opt$filter_duplicates_missingness #args[7]

#project = gsub("-","_",str_to_lower(gsub("\\[|\\]", "", project)))
print(paste("Project:",project))
print(paste("Methylation file:",methpath))

# Read methylation data
meth_df = fread(methpath,data.table=F)

# Probe list
probe_list = meth_df$probeID
# Get NetsciDataCompanion object
my_friend = NetworkDataCompanion::CreateNetworkDataCompanionObject(project_name=project)

# Read probe map
my_map = data.frame(fread(probe_map_fn,sep=",",header=T),row.names=1)
# get only probes annotated to transcription factors at this point
# download.file('http://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt', 
#              destfile = "TF_names_v_1.01.txt")

# Get all gene names
#names(my_map)[2]="geneName"
#print(head(my_map))
#print('all genes')
all_genes = unique(na.omit(my_map$geneNames))
#print(head(all_genes))

if(filter_duplicates_missingness)
{
  keep_uuids = my_friend$filterDuplicatesMethylationMissingness(meth_df)
  print(paste(c("Keeping",length(keep_uuids),"samples after filtering duplicates")))
  meth_df <- meth_df[, c("probeID", keep_uuids)]
}

# filter methylation beta to only probes that are in the tf list
# Manage list of TF
if (tf_list_fn!=' '){

  # Read TF list
  print(paste(c('Reading TF list', tf_list_fn)), collapse = " ")
  tf_list = read.table(tf_list_fn)[,1] #read.table(paste(baseDir,"ext/TF_names_v_1.01.txt",sep="/"))[,1]
  print(head(tf_list))

  tf_list_filtered = intersect(tf_list,all_genes)
  print(paste(c('There are', length(tf_list_filtered), 'genes that are going to be used'), collapse = ' '))

  gene_map = my_friend$probeToMeanPromoterMethylation(methylation_betas = meth_df, 
                                    probe_gene_map = my_map, 
                                    genesOfInterest = tf_list_filtered)
} else {
  # Use all genes
  gene_map = my_friend$probeToMeanPromoterMethylation(methylation_betas = meth_df, 
                                   probe_gene_map = my_map, genesOfInterest = NULL)
}


# Get TCGA barcodes
barcodes = my_friend$mapUUIDtoTCGA(row.names(gene_map))
# sanity check
print("gene map row names")
print(row.names(gene_map)[1:5])
print("Sanity check on order (should be 0)")
sum(order(row.names(barcodes)) != order(row.names(gene_map)))

# add TCGA and write to file 
gene_map %>% as.data.frame() %>%
  mutate(TCGA_barcode = barcodes$submitter_id) %>%
    write.csv(output_fn,row.names=T, quote=F)