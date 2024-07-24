library(TCGAbiolinks)
library(optparse)
library(SummarizedExperiment)

option_list = list(
  make_option(c("-p", "--project_id"), type="character", default=NULL, 
              help="project_id", metavar="character"),
  make_option(c("-d", "--downloads_dir"), type="character", default=NULL, 
              help="downloads_dir", metavar="character"),
  make_option(c("--analysis_workflow_type"), type="character", default="ASCAT3", 
              help="", metavar="character"),  
  make_option(c("--output_table"), type="character", default=NULL, 
              help="", metavar="character"),
  make_option(c("--output_rds"), type="character", default=NULL, 
              help="", metavar="character"),
  make_option(c("--sample_list"), type="character", default=' ', 
              help="test file with sample names", metavar="character")    
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt)

project_id = opt$project_id # args[1]
analysis_workflow_type = opt$analysis_workflow_type # always methylation_beta_value at this point
print(paste('Downloading',project_id,',with workflow:',analysis_workflow_type))
output_table = opt$output_table
output_rds = opt$output_rds

#barcode = listSamples
if (nchar(opt$sample_list)>3){
    
  queried <- GDCquery(
    project = project_id, 
    data.category = "Copy Number Variation", 
    data.type = "Gene Level Copy Number",
    workflow.type = analysis_workflow_type)
  } else {
  queried <- GDCquery(
    project = project_id, 
    data.category = "Copy Number Variation", 
    data.type = "Gene Level Copy Number",
    workflow.type = analysis_workflow_type
  ) }


GDCdownload(queried,files.per.chunk = 50)
cn_rds <- GDCprepare(queried)

if (nchar(opt$sample_list)>3){
    barcodes = read.table(opt$sample_list, header = FALSE, sep = "", dec = ".")
    nnn = nchar(barcodes$V1[1])
    # Check that this is correct
    if (nnn>12){nnn=12}
    # select based on nnn
    barcodes = substr(barcodes$V1,1,nnn)
    cn_rds = cn_rds[,substr(cn_rds$sample,1,nnn) %in% barcodes]
}


saveRDS(cn_rds, output_rds)
write.csv( assays(cn_rds)$copy_number, file = output_table)