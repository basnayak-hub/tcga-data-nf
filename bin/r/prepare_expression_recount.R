library("recount3")
library('NetSciDataCompanion')
library("optparse")

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



option_list = list(
  make_option(c("-p", "--project"), type="character", default=NULL, 
              help="project name (e.g. luad)", metavar="project_name"),
  make_option(c("-c", "--clinical"), type="character", default=NULL, 
              help="patient data filename", metavar="patient_data"),
  make_option(c("-e", "--expression"), type="character", default=NULL, 
              help="expression data rds", metavar="expression_rds"),
  make_option(c("-r", "--output_rds"), type="character", default='filtered_expression.rds', 
              help="output filename for rds", metavar="output_rds"),
  make_option(c("-t", "--output_table"), type="character", default='filtered_expression.txt', 
              help="output filename for table", metavar="output_txt"),
  make_option(c("--th_purity"), type="double", default=0.7,
              help="purity threshold. Default:0.7", metavar="number"),
  make_option(c("--min_tpm"), type="double", default=1.0, 
              help="min tpm for filtering. Default:1", metavar="number"),
  make_option(c("--frac_samples"), type="double", default=0.1, 
              help="fraction of samples that need to have a min tpm. Default:0.1", metavar="number"),
make_option(c("--tissue_type"), type="character", default='all', 
              help="tissue type to be extracted. Default:all", metavar="character"),
make_option(c("--normalization"), type="character", default='logtpm', 
              help="normalization strategy, 'logtpm', 'tpm', 'counts'. Default:logtpm", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

project = opt$project #args[1]
patient_data = opt$clinical#args[2]
exp_data = opt$expression#args[3]
output_rds = opt$output_rds#args[4]
output_table = opt$output_table#args[5]
th_purity = as.numeric(opt$th_purity)
min_tpm = as.numeric(opt$min_tpm)
frac_samples = as.numeric(opt$frac_samples)
tissue_type = opt$tissue_type 
normalization = opt$normalization


project_name <- toupper(substring(project, 6))
print(paste0('Preparing Project:', project_name))
#patient_data <- "clinical_patient_luad.csv"
#exp_data <- "tcga_luad.rds"
#mut_data <- "tcga_luad_mutations.txt"
#mut_pivot_data <- "tcga_luad_mutations_pivot.csv"
#meth_data <- "tcga_luad_methylations.txt"


# Definitions
tumor_tissues=c("Primary Solid Tumor",
                  "Recurrent Solid Tumor",
                  "Primary Blood Derived Cancer - Peripheral Blood",
                  "Recurrent Blood Derived Cancer - Bone Marrow",
                  "Additional - New Primary",
                  "Metastatic","Additional Metastatic",
                  "Human Tumor Original Cells",
                  "Primary Blood Derived Cancer - Bone Marrow")

normal_tissues=c("Solid Tissue Normal",
                  "Buccal Cell Normal",
                  "EBV Immortalized Normal",
                  "Bone Marrow Normal",
                  "sample type 15",
                  "sample type 16")

# Create Companion Object
print('Creating Object from...')
print(paste0('Patient data:',patient_data))
obj <- CreateNetSciDataCompanionObject(clinical_patient_file = patient_data,
                                       project_name = project_name)

print(paste0('Reading expression:',exp_data))
# Read RDS expression data
test_exp_rds <- readRDS(exp_data)

# Extract included meta info from RDS expression object
rds_info <- obj$extractSampleAndGeneInfo(test_exp_rds)
rds_sample_info <- rds_info$rds_sample_info
rds_gene_info <- rds_info$rds_gene_info

# Normalize data
print('Normalize...')
test_exp_all <- obj$logTPMNormalization(test_exp_rds)
test_exp_count <- test_exp_all$counts
test_exp_tpm <- test_exp_all$TPM
test_exp_logtpm <- test_exp_all$logTPM


# Map column names to TCGA barcodes
print('Map Barcodes...')
newcolnames <- obj$mapUUIDtoTCGA(colnames(test_exp_logtpm))
colnames(test_exp_count) <- newcolnames[,2]
colnames(test_exp_tpm) <- newcolnames[,2]
colnames(test_exp_logtpm) <- newcolnames[,2]

# Get indices of nonduplicates
idcs_nonduplicate <- obj$filterDuplicatesSeqDepth(expression_count_matrix = test_exp_count)

# Get indices of genes that have a minimum TPM in the data.
# Here we are filtering by gene
test_exp_tpm_df <- data.frame(test_exp_tpm)
idcs_genes_mintpm <- obj$filterGenesByTPM(test_exp_tpm_df, min_tpm, frac_samples)

# Select the normalization strategy to be saved
if (normalization=='logtpm'){
    test_exp_final = test_exp_logtpm
    } else if (normalization=='tpm') {
   test_exp_final = test_exp_tpm
} else if (normalization=='count') {
   test_exp_final = test_exp_count
} else {
  stop("ERROR: normalization unknown")
} 

# Get indices of samples passing purity filtering
idcs_purity <- obj$filterPurity(colnames(test_exp_final), threshold = th_purity)

print("Tissue type: ")
print(tissue_type)

print('Get final data...')
# Get indices of samples, we have to distinguish between tumor, normal, all samples
if (tissue_type=='all'){
  # filter samples in specific tissue
  idcs_tumor <- obj$filterTumorSamples(colnames(test_exp_final))
  idcs_normal <- obj$filterNormalSamples(colnames(test_exp_final))
  # select the tumor samples that pass purity thresholds
  idcs_tumor_pure <- intersect(idcs_tumor, idcs_purity)
  # Select samples that are: nonduplicate, normal and tumor filtered by purity
  idcs_final <- intersect( idcs_nonduplicate, union(idcs_tumor_pure, idcs_normal))
} else if (tissue_type %in% tumor_tissues){
  # filter samples in specific tissue
  idcs_sample <- obj$filterTumorType(colnames(test_exp_final), tissue_type, rds_sample_info)
  # select the tumor samples that pass purity thresholds
  idcs_final <- intersect( idcs_nonduplicate, intersect(idcs_sample, idcs_purity))
  } else if (tissue_type %in% normal_tissues) {
  # filter samples in specific normal tissue
  idcs_sample <- obj$filterTumorType(colnames(test_exp_final), tissue_type, rds_sample_info)
  idcs_final <- intersect( idcs_nonduplicate, idcs_sample)
} else if  (tissue_type == 'tumor') {
  # select all tumor samples
  idcs_sample <- obj$filterTumorSamples(colnames(test_exp_final))
  idcs_final <- intersect( idcs_nonduplicate, intersect(idcs_sample, idcs_purity))
} else if (tissue_type == 'normal') {
  # select all normal samples
  idcs_sample <- obj$filterNormalSamples(colnames(test_exp_final))
  idcs_final <- intersect( idcs_nonduplicate, idcs_sample)
} else {
   stop('tissue type not recognised')
}

  # get final datasets
  final_table = test_exp_final[idcs_genes_mintpm, idcs_final]
  final_rds = test_exp_rds[idcs_genes_mintpm, idcs_final]

print('Saving RDS...')
saveRDS(final_rds, file = output_rds)

print('Saving table...')
write.table(final_table, file = output_table, sep = '\t', quote = F, col.names=NA)
print('DONE!')