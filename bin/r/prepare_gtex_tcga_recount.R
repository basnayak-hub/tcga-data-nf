library("recount3")
library("recount")
library("data.table")
library(readr)

args<-commandArgs(TRUE)

# args[1] = input for thga
# args[2] = output for tcha
# args[3] = input for gtex
# args[4] = output for gtex
# args[5] = normalization (TPM)
# args[6] = threshold for removal, defaults to 0.2 (20%)

in_tcga_fn = args[1]
out_tcga_fn = args[2]
in_gtex_fn = args[3]
out_gtex_fn = args[4]

if( is.null(args[5]) ){
  norm<-'TPM'} else{
   norm = args[5]
}

if(is.null(args[6])){threshold_remove<-0.2} else {
   threshold_remove = as.numeric(args[6])
}

# Create temp file location
#tmp <- file.path('.', log_fn)
# Open log
#lf <- log_open(tmp)

print('----------------------------------------------')
print('# Inputs:')
print(paste0("# GTEX-TCGA data"))
print('## Input arguments:')
print(paste0("TCGA input file: ",in_tcga_fn))
print(paste0('GTEx input file: ',in_gtex_fn))
print(paste0("TCGA output file: ",out_tcga_fn))
print(paste0('GTEx output file: ',out_gtex_fn))
print(paste0("##Parameters:"))
print(paste0("Norm: ",norm))
print(paste0("threshold remove: ",threshold_remove))
print('----------------------------------------------')

## Load the data that have been download already
print(paste0("Reading rds data..."))
GTEx_lung = readRDS(in_gtex_fn)
TCGA_lung = readRDS(in_tcga_fn)

### ## Once you have your RSE object, you can transform the raw coverage
### ## base-pair coverage counts using transform_counts().
### ## For RPKM, TPM or read outputs, check the details in transform_counts().

### ## Compute TPMs
## #TODO: here we use TPM as the final output, substitute with a more general name
if (norm=='TPM') {
  print(paste0("Computing TPM for gtex..."))
  assays(GTEx_lung)$counts <- recount3::transform_counts(GTEx_lung)
  assays(GTEx_lung)$TPM <- recount::getTPM(GTEx_lung)
  #colSums(assay(GTEx_lung, "TPM")) / 1e6 ## Should all be equal to 1
  
  print(paste0("Computing TPM for tcga..."))
  assays(TCGA_lung)$counts <- recount3::transform_counts(TCGA_lung)
  assays(TCGA_lung)$TPM <- recount::getTPM(TCGA_lung)
  #colSums(assay(TCGA_lung, "TPM")) / 1e6 ## Should all be equal to 1
  } else {
   print(paste0("Computing raw counts for gtex..."))
   assays(GTEx_lung)$TPM <- recount3::transform_counts(GTEx_lung)

  print(paste0("Computing raw counts for tcga..."))
  assays(TCGA_lung)$TPM <- recount3::transform_counts(TCGA_lung)
  }

# Check if GTEx & TCGA have the same set of genes
# proportion of rows where gene ids match in GTEx & TCGA
print(paste0("Checking integrity..."))
sum(GTEx_lung@rowRanges$gene_id == TCGA_lung@rowRanges$gene_id)/length(GTEx_lung@rowRanges$gene_id == TCGA_lung@rowRanges$gene_id)

print(paste0("Removing below threshold ..."))
# Remove genes with counts < 1TPM in at least 20% samples (ans. 25170 genes)
expression_all = cbind(assays(GTEx_lung)$TPM,assays(TCGA_lung)$TPM)
expression_cutoff = 1 # 1 TPM
minSamples = threshold_remove*ncol(expression_all) # at least 20% of samples
keep = rowSums(expression_all > expression_cutoff) >= minSamples
table(keep)
expression_filtered = expression_all[keep,]

print(paste0(nrow(expression_all), ":Total number of genes overlapped between TCGA and GTEx"))
print(paste0(nrow(expression_filtered), ":Number of genes after filtering ", expression_cutoff," cpm in ", minSamples, " samples"))

print("Filtering by TPM...")
# Filtered TPM matrices
GTEx_TPM_filtered = expression_filtered[,1:ncol(assays(GTEx_lung)$TPM)]
rownames(GTEx_TPM_filtered) = rownames(expression_filtered)
colnames(GTEx_TPM_filtered) = colnames(assays(GTEx_lung)$TPM)
TCGA_TPM_filtered = expression_filtered[,-(1:ncol(assays(GTEx_lung)$TPM))]
rownames(TCGA_TPM_filtered) = rownames(expression_filtered)
colnames(TCGA_TPM_filtered) = colnames(assays(TCGA_lung)$TPM)

print("Writing to output...")

write.csv(GTEx_TPM_filtered[1:5,1:5],file=out_gtex_fn)
write.csv(TCGA_TPM_filtered[1:5,1:5],file=out_tcga_fn)
