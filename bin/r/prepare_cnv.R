options(warn = -1)
library("optparse")
library(dplyr)
library(stringr)
library(huge)
library(SummarizedExperiment)

#/Users/violafanfani/Documents/uni-harvard/projects/tcga-data-supplement/data/processed/batch-coad-subtype-20240510/tcga_coad_cms1/data_download/cnv/tcga_coad_cms1.csv

option_list = list(
  make_option(c("-p", "--project"), type="character", default=NULL, 
              help="project name (e.g. tcga_luad)", metavar="project_name"),
  make_option(c("-c", "--cnv"), type="character", default=NULL, 
              help="cnv data rds", metavar="cnv_rds"),
  #make_option(c("-r", "--cnv_rds"), type="character", default="", 
  #            help="cnv data rds", metavar="cnv_rds"), 
  make_option(c("-o", "--output"), type="character", default='output_cnv.csv', 
              help="output filename for csv", metavar="output_csv"),
  make_option(c("--th_std"), type="double", default=0.0,
              help="std threshold, removes everything with std lower than threshold. By default removes everything that has zero std. Default:0.0", metavar="number"),
  make_option(c("--tf_list"), type="character", default=" ", 
              help="TF list filename. Pass a text file with TF names to filter the output", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

project = opt$project #args[1]
cnv_fn = opt$cnv#args[2]
output_fn = opt$output#args[3]
th_std = as.numeric(opt$th_std)
tf_list_fn =  opt$tf_list


project_name <- toupper(substring(project, 6))
#patient_data <- "clinical_patient_luad.csv"
#exp_data <- "tcga_luad.rds"
#mut_data <- "tcga_luad_mutations.txt"
#mut_pivot_data <- "tcga_luad_mutations_pivot.csv"
#meth_data <- "tcga_luad_methylations.txt"

print('We start pre-processing the CNV')
print(paste0('Project:', project_name))



# read cnv data
#cnv = read.csv(cnv_fn,check.names = F, row.names=1)
cnv_rds = readRDS(cnv_fn)

if (tf_list_fn!=' '){
  # Read TF list
  print(paste(c('Reading TF list', tf_list_fn)), collapse = " ")
  tf_list = read.table(tf_list_fn)[,1] #read.table(paste(baseDir,"ext/TF_names_v_1.01.txt",sep="/"))[,1]
  print(head(tf_list))

  # Check if entry of tf_list starts with "ENSG"
  if (tf_list[1] %>% str_detect("^ENSG")) {
    # If it does, we need to remove the "ENSG" part
    print("Genes are in ENSG format, no need to change them now")
    # Select the genes that are in the tf_list
    cnv_rds = cnv_rds[cnv_rds@rowRanges$gene_id %in% tf_list]
    # Print the new length
    print(paste(c('There are', length(rownames(cnv_rds)), 'genes that are going to be used'), collapse = ' '))   
  } else {
    # Select the genes that are in the tf_list
    cnv_rds = cnv_rds[cnv_rds@rowRanges$gene_name %in% tf_list]
    rownames(cnv_rds) = cnv_rds@rowRanges$gene_name
    # Print the new length
    print(paste(c('There are', length(rownames(cnv_rds)), 'genes that are going to be used'), collapse = ' '))
  }
}
cnv = assays(cnv_rds)$copy_number

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
  relocate(tcga_sample_barcode)  %>% 
  dplyr::select(-p1,-p2,-code1,-code2,-code1_num, -code2_num)# Drop intermediate columns

# Let's actually check that the code is tumor
# We can do this by checking the last two characters of the code
# If they are less than 10 we can assume it's a tumor sample

cnv_t_labeled = cnv_t_labeled %>%
  mutate(
    code1 = substr(tcga_sample_barcode,start=14,stop=16),  # Split by semicolon
    code1_num = as.numeric(substr(code1, 1, 2)),           # Extract the numeric part of the first code
  ) 

q = cnv_t_labeled %>% dplyr::filter(code1_num>9)

if (nrow(q)>0){
  print("WARNING: one or more samples discarded because non-tumor")
  print( q$tcga_sample_barcode )
  
  cnv_t_labeled = cnv_t_labeled %>%
  dplyr::filter(code1_num<10)
}

cnv_t_labeled = cnv_t_labeled%>% 
  dplyr::select(-code1, -code1_num)# Drop intermediate columns


# And then we check if there are any duplicated samples
# If there are, we take the mean of the values for each gene across the samples
# This is a sanity check, but it should not happen if the data is clean
if (anyDuplicated(cnv_t_labeled$tcga_sample_barcode)) {
  print('WARNING: duplicated samples')
  cnv_t_labeled <- cnv_t_labeled %>%
    group_by(tcga_sample_barcode) %>%
    summarize(across(everything(), ~ mean(.x, na.rm = TRUE), .groups = "drop"))
}

# replace the rownames with the tcga_sample_barcode
rownames(cnv_t_labeled) = cnv_t_labeled$tcga_sample_barcode
# remove the TCGAbarcode column
cnv_t_labeled = cnv_t_labeled %>%
  dplyr::select(-tcga_sample_barcode)



# Remove any gene with missingness (NO imputation done) and follow this by removing any gene with zero standard deviation (same copy number across all samples), as that breaks DRAGON. 

cnv_t_no_miss = cnv_t_labeled %>% 
  dplyr::select(-which(apply(.,2,function(x){sum(is.na(x))>0}))) %>%
  dplyr::select(-which(apply(.,2,function(x){sd(x)<=th_std})))
  
# Check that there are no NAs in the data (although we should remove them before)
stopifnot(sum(is.na(cnv_t_no_miss))==0)

# Apply nonparanormal transformation
cnv_t_npn = huge.npn(cnv_t_no_miss)

# Create a data frame with the transformed data
cnv_final = cnv_t_npn %>%
  data.frame()
cnv_final$TCGAbarcode = rownames(cnv_final)

# get final datasets
print('Saving table...')
print(paste('LOG:',"Table dimensions:",dim(cnv_final)[1],"x",dim(cnv_final)[2], sep = " ", ""))

write.table(cnv_final, file = output_fn, sep = ',', quote = F, col.names=T, row.names = F)
print('DONE!')

