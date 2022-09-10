library("recount3")
library('NetSciDataCompanion')

args<-commandArgs(TRUE)

# args[1] = project for gtex
# args[2] = project for tcga
# args[3] = organism
# args[4] = gencode version
# args[5] = type

project = args[1]
patient_data = args[2]
exp_data = args[3]
output_rds = args[4]
output_table = args[5]
th_purity = as.numeric(args[6])
min_tpm = as.numeric(args[7])
frac_samples = as.numeric(args[8])
tissue_type = args[9] 


print(args)
print('Downloading...')

project_name <- toupper(substring(project, 6))
print(project_name)
#patient_data <- "clinical_patient_luad.csv"
#exp_data <- "tcga_luad.rds"
#mut_data <- "tcga_luad_mutations.txt"
#mut_pivot_data <- "tcga_luad_mutations_pivot.csv"
#meth_data <- "tcga_luad_methylations.txt"


# Create Companion Object
print('Creating Object...')
print(patient_data)
print(project_name)
obj <- CreateNetSciDataCompanionObject(clinical_patient_file = patient_data,
                                       project_name = project_name)


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

# Get indices of samples that are primary tumor
idcs_sample <- obj$filterTumorType(colnames(test_exp_count), tissue_type, rds_sample_info)

# Get indices of nonduplicates
idcs_nonduplicate <- obj$filterDuplicatesSeqDepth(expression_count_matrix = test_exp_count[,idcs_sample])

# Get indices of genes that have a minimum TPM in the data
test_exp_tpm_df <- data.frame(test_exp_tpm)
idcs_genes_mintpm <- obj$filterGenesByTPM(test_exp_tpm_df, min_tpm, frac_samples)

print('Get final data...')
# Get indices of samples passing purity filtering
if (tissue_type=='Primary Tumor'){
  idcs_purity <- obj$filterPurity(colnames(test_exp_count), threshold = th_purity)
  final_table = test_exp_logtpm[,idcs_sample][idcs_genes_mintpm, intersect(idcs_nonduplicate, idcs_purity[idcs_sample])]
  final_rds = test_exp_rds[,idcs_sample][idcs_genes_mintpm, intersect(idcs_nonduplicate, idcs_purity[idcs_sample])]
} else if (tissue_type=='Solid Tissue Normal') {
  final_table = test_exp_logtpm[,idcs_sample][idcs_genes_mintpm, idcs_nonduplicate]
  final_rds = test_exp_rds[,idcs_sample][idcs_genes_mintpm, idcs_nonduplicate]
} else {
   print('tissue type not recognised')
}


print('Saving RDS...')
saveRDS(final_rds, file = output_rds)

print('Saving table...')
write.table(final_table, file = output_table, sep = '\t', quote = F, col.names=NA)
print('DONE!')