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

names(expr_std_labeled) = gene_name_map$gene_name

df <- read.csv(file=dragon_tf,nrows=2)  
tfs = colnames(df)[2:length(colnames(df))-1]
expr_std_labeled <- expr_std_labeled %>% dplyr::select(any_of(tfs))

expr_std_labeled$TCGAbarcode = substr(expr_rds@colData$tcga.tcga_barcode, 1, 16)

write.csv(expr_std_labeled,file=output_fn,row.names=T,quote=F)