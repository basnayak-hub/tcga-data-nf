# get expression data from s3 and clean
library(tidyverse)
library(data.table)
library(SummarizedExperiment)
args = commandArgs(trailingOnly=TRUE)

expr_rds = args[1]
dragon_tf = args[2]
output_fn = args[3]

expr_rds = readRDS(expr_rds)

gene_name_map = data.frame("ENSG"=row.names(rowData(expr_rds)),
                           "gene_name"=rowData(expr_rds)$gene_name)

if (length(assays(expr_rds))==2){
  expr = assays(expr_rds)[[2]]
} else {
  expr = assays(expr_rds)[[1]]
}    


# format: genes in rows, samples in columns 

# for dragon they should be centered scaled 
summary(apply(expr[1:10,],1,mean))
summary(apply(expr[1:10,],1,sd))

standardize = function(x){return((x-mean(x))/sd(x))}

expr_std = apply(expr,1,standardize)

# this transposes the matrix
dim(expr_std)

# sanity check
summary(apply(expr_std,2,mean))
summary(apply(expr_std,2,sd))

expr_std_labeled = expr_std %>% as.data.frame()



# Read header of cnv/methylation data to match the names
df <- read.csv(file=dragon_tf,nrows=2)
tfs = colnames(df)[2:length(colnames(df))-1]
# Check if tfs are in symbol or ENSG format
if (tfs[1] %>% str_detect("^ENSG")) {
  # If it does, we need to remove the "ENSG" part
  print("Genes are in ENSG format, no need to change them now")
} else {
  # Select the genes that are in the tf_list
  names(expr_std_labeled) = gene_name_map$gene_name
  print("Genes are in symbol format, translating them")
}

expr_std_labeled <- expr_std_labeled %>% dplyr::select(any_of(tfs))

expr_std_labeled$TCGAbarcode = substr(expr_rds@colData$tcga.tcga_barcode, 1, 16)

write.csv(expr_std_labeled,file=output_fn,row.names=T,quote=F)