options(warn = -1)
library("optparse")
library(dplyr)
library(stringr)
library(huge)

#args<-commandArgs(TRUE)

# args[1] = project, tcga_luad
# args[2] = clinical data, patient.csv
# args[3] = expression data, expression.rds
# args[4] = output rds file
# args[5] = output txt file
# args[6] = purity threshold
# args[7] = minimum tpm
# args[8] = fraction of samples for filtering
# args[9] = tissue to be retrieved, 'Primary Tumor'/ 'Solid Tissue Normal'
# args[10] = normalization strategy, 'logTPM', 'TPM', 'counts', 
# args[11] = batch correction variable, 'logTPM', 'TPM', 'counts', 
# args[10] = adjustment variable, 'logTPM', 'TPM', 'counts', 

#/Users/violafanfani/Documents/uni-harvard/projects/tcga-data-supplement/data/processed/batch-coad-subtype-20240510/tcga_coad_cms1/data_download/cnv/tcga_coad_cms1.csv

option_list = list(
  make_option(c("-p", "--project"), type="character", default=NULL, 
              help="project name (e.g. tcga_luad)", metavar="project_name"),
  make_option(c("-c", "--cnv"), type="character", default=NULL, 
              help="cnv data csv", metavar="cnv_cvs"),
  make_option(c("-r", "--cnv_rds"), type="character", default="", 
              help="cnv data rds", metavar="cnv_rds"), 
  make_option(c("-o", "--output"), type="character", default='output_cnv.csv', 
              help="output filename for csv", metavar="output_csv"),
  make_option(c("--th_std"), type="double", default=0.0,
              help="std threshold, removes everything with std lower than threshold. By default removes everything that has zero std. Default:0.0", metavar="number")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

project = opt$project #args[1]
cnv_fn = opt$cnv#args[2]
output_fn = opt$output#args[3]
th_std = as.numeric(opt$th_std)

project_name <- toupper(substring(project, 6))
#patient_data <- "clinical_patient_luad.csv"
#exp_data <- "tcga_luad.rds"
#mut_data <- "tcga_luad_mutations.txt"
#mut_pivot_data <- "tcga_luad_mutations_pivot.csv"
#meth_data <- "tcga_luad_methylations.txt"

print('We start pre-processing the CNV')
print(paste0('Project:', project_name))

# read cnv data
cnv = read.csv(cnv_fn,check.names = F, row.names=1)

cnv_t = cnv %>%
  t() %>%
  data.frame()

# We first need to check the barcodes of the samples and make sure 
# all samples are matched to the same patient
cnv_t_labeled = cnv_t%>%
  mutate(tcga_patient1_barcode = substr(str_split(row.names(.),";",simplify=T)[,1],start=1,stop=12) ) %>%
  relocate(tcga_patient1_barcode) 

cnv_t_labeled = cnv_t_labeled%>%
  mutate(tcga_patient2_barcode = substr(str_split(row.names(.),";",simplify=T)[,2],start=1,stop=12) ) %>%
  relocate(tcga_patient2_barcode)
  
# Filter only samples with matched patient
cnv_t_labeled = cnv_t_labeled %>%
  dplyr::filter(tcga_patient1_barcode == tcga_patient2_barcode) %>%
  dplyr::select(-tcga_patient1_barcode, -tcga_patient2_barcode)%>%
  dplyr::mutate(tcga_sample_barcode = rownames(.)) %>%
  relocate(tcga_sample_barcode)


# For each patient we select the cancer sample

cnv_t_labeled = cnv_t_labeled %>%
  mutate(
    p1 = substr(str_split(row.names(.),";",simplify=T)[,1],start=1,stop=16),  # Split by semicolon
    p2 = substr(str_split(row.names(.),";",simplify=T)[,2],start=1,stop=16),  # Split by semicolon
    code1 = substr(str_split(row.names(.),";",simplify=T)[,1],start=14,stop=16),  # Split by semicolon
    code2 = substr(str_split(row.names(.),";",simplify=T)[,2],start=14,stop=16),  # Split by semicolon
    code1_num = as.numeric(substr(code1, 1, 2)),           # Extract the numeric part of the first code
    code2_num = as.numeric(substr(code2, 1, 2))            # Extract the numeric part of the second code
  ) %>%
  mutate(
    tcga_sample_barcode = if_else(code1_num < code2_num, p1, p2)) %>% # Select the code with the smaller numeric part) 
  relocate(tcga_sample_barcode)  # Drop intermediate columns

# And then we check if there are any duplicated samples
# If there are, we take the mean of the values for each gene across the samples
# This is a sanity check, but it should not happen if the data is clean

if (anyDuplicated(cnv_t_labeled$tcga_sample_barcode)) {
  print('WARNING: duplicated samples')
  cnv_t_labeled <- cnv_t_labeled %>%
    group_by(tcga_sample_barcode) %>%
    summarize(across(everything(), ~ mean(.x, na.rm = TRUE), .groups = "drop"))
}


# Remove any gene with missingness (NO imputation done) and follow this by removing any gene with zero standard deviation (same copy number across all samples), as that breaks DRAGON. 

cnv_t_no_miss = cnv_t_labeled %>% 
  dplyr::select(-which(apply(.,2,function(x){sum(is.na(x))>0}))) %>%
  dplyr::select(-which(apply(.,2,function(x){sd(x)<=th_std})))


# Apply nonparanormal transformation
cnv_t_npn = huge.npn(cnv_t_no_miss)

# sanity check that npn is operating on marginal level
# it requires more than one variable, not sure why - looks like
# just something with the vectorization of apply
a = max(abs(cnv_t_npn[,1:2] -huge.npn(cnv_t_no_miss[,1:2])))
stopifnot(a == 0)

# Convert rownames to TCGA barcode,by splitting row names and shortening to base TCGA barcode (removing replicate and sample type info).

# get final datasets
print('Saving table...')
print(paste('LOG:',"Table dimensions:",dim(cnv_t_npn)[1],"x",dim(cnv_t_npn)[2], sep = " ", ""))

write.table(cnv_t_npn, file = output_fn, sep = ',', quote = F, col.names=NA)
print('DONE!')

