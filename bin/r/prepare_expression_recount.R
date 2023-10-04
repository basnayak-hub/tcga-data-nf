library("recount3")
library('NetSciDataCompanion')
library("optparse")
library("limma")
library("ggplot2") 
library("cowplot")
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



option_list = list(
  make_option(c("-p", "--project"), type="character", default=NULL, 
              help="project name (e.g. tcga_luad)", metavar="project_name"),
  make_option(c("-c", "--clinical"), type="character", default=NULL, 
              help="patient data filename", metavar="patient_data"),
  make_option(c("-e", "--expression"), type="character", default=NULL, 
              help="expression data rds", metavar="expression_rds"),
  make_option(c("-r", "--output_rds"), type="character", default='filtered_expression.rds', 
              help="output filename for rds", metavar="output_rds"),
  make_option(c("-t", "--output_table"), type="character", default='filtered_expression.txt', 
              help="output filename for table", metavar="output_txt"),
  make_option(c("-f", "--output_pca"), type="character", default='output pca ', 
              help="output filename for the pca figure", metavar="output_txt"),
  make_option(c("--th_purity"), type="double", default=0.7,
              help="purity threshold. Default:0.7", metavar="number"),
  make_option(c("--min_tpm"), type="double", default=1.0, 
              help="min tpm for filtering. Default:1", metavar="number"),
  make_option(c("--frac_samples"), type="double", default=0.1, 
              help="fraction of samples that need to have a min tpm. Default:0.1", metavar="number"),
make_option(c("--tissue_type"), type="character", default='all', 
              help="tissue type to be extracted. Default:all", metavar="character"),
make_option(c("--normalization"), type="character", default='logtpm', 
              help="normalization strategy, 'logtpm', 'tpm', 'counts', 'logCPM'. Default:logtpm", metavar="character"),
make_option(c("--batch_correction"), type="character", default='', 
              help="batch correction variable, By default is none. Default: '' ", metavar="character"),
make_option(c("--adjustment_variable"), type="character", default='', 
              help="adjustment variable for batch correction. By default is empty. Default: '' ", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

project = opt$project #args[1]
patient_data = opt$clinical#args[2]
exp_data = opt$expression#args[3]
output_rds = opt$output_rds#args[4]
output_table = opt$output_table#args[5]
output_pca = opt$output_pca#args[5]
th_purity = as.numeric(opt$th_purity)
min_tpm = as.numeric(opt$min_tpm)
frac_samples = as.numeric(opt$frac_samples)
tissue_type = opt$tissue_type 
normalization = opt$normalization
to_batch_correct_nominal_name = opt$batch_correction
adjustment_variable = opt$adjustment_variable

project_name <- toupper(substring(project, 6))
print(paste0('Preparing Project:', project_name))
#patient_data <- "clinical_patient_luad.csv"
#exp_data <- "tcga_luad.rds"
#mut_data <- "tcga_luad_mutations.txt"
#mut_pivot_data <- "tcga_luad_mutations_pivot.csv"
#meth_data <- "tcga_luad_methylations.txt"

print('Hereeeeee')

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

with_purity = c("ACC","BLCA","BRCA","CESC","COAD","GBM",
              "HNSC","KIRC","KIRP","KICH","LGG","LIHC",
              "LUAD","LUSC","OV","PRAD","READ","SKCM",
              "THCA","UCEC","UCS")



gtex_projects =c('LUNG', 'BRAIN', 'SKIN', 'ESOPHAGUS', 'BLOOD_VESSEL', 'ADIPOSE_TISSUE',
 'BLOOD', 'HEART', 'MUSCLE', 'COLON', 'THYROID', 'NERVE', 'LUNG', 'BREAST', 'TESTIS', 
 'STOMACH', 'PANCREAS', 'PITUITARY', 'ADRENAL_GLAND', 'PROSTATE', 'SPLEEN', 'LIVER', 
 'BONE_MARROW', 'OVARY', 'SMALL_INTESTINE', 'SALIVARY_GLANDS', 'VAGINA', 'UTERUS', 
 'KIDNEY', 'BLADDER', 'CERVIX_UTERI', 'FALLOPIAN_TUBE')

# Create Companion Object
print('Creating Object from...')
print(paste0('Patient data:',patient_data))
print(project_name %in% gtex_projects)

if (project_name %in% gtex_projects){
  obj <- CreateNetSciDataCompanionObject()
} else {
obj <- CreateNetSciDataCompanionObject( project_name = project_name)
}


print(paste0('Reading expression:',exp_data))
# Read RDS expression data
test_exp_rds <- readRDS(exp_data)

# Extract included meta info from RDS expression object
rds_info <- obj$extractSampleAndGeneInfo(test_exp_rds)
rds_sample_info <- rds_info$rds_sample_info
rds_gene_info <- rds_info$rds_gene_info

# Map column names to TCGA barcodes
print('Map Barcodes...')
if (normalization %in% c('count','tpm','logtpm') ){
  # Normalize data
  print('Normalize with TPM...')
  # Here we call it xpm, so it's the same for tpm and cpm
  test_exp_all <- obj$logTPMNormalization(test_exp_rds)
  test_exp_count <- test_exp_all$counts
  test_exp_xpm <- test_exp_all$TPM
  test_exp_logxpm <- test_exp_all$logTPM
} else if (normalization %in% c('cpm','logcpm')) {
  # Normalize data
  print('Normalize with CPM...')
  # Here we call it xpm, so it's the same for tpm and cpm
  test_exp_all <- obj$logCPMNormalization(test_exp_rds)
  test_exp_count <- test_exp_all$counts
  test_exp_xpm <- test_exp_all$CPM
  test_exp_logxpm <- test_exp_all$logCPM
} else {
    stop("ERROR: normalization unknown")
  } 
  
  # assign correct names
  print(paste('LOG:',"There are",length(colnames(test_exp_logxpm)),"intial columns", sep = " ", ""))

  if (project_name %in% gtex_projects){
    print('Processing GTEX, no need to remap sample names')
  } else {
    print(colnames(test_exp_logxpm)[1:5])
    newcolnames <- test_exp_rds@colData$tcga.tcga_barcode
    #newcolnames <- obj$mapUUIDtoTCGA(colnames(test_exp_logxpm))
    print(newcolnames)

    print(paste('LOG:',"There are",dim(newcolnames),"mapped columns", sep = " ", " "))
    colnames(test_exp_count) <- newcolnames#[,2]
    colnames(test_exp_xpm) <- newcolnames#[,2]
    colnames(test_exp_logxpm) <- newcolnames#[,2]
  }

  #### Filters
  # Get indices of nonduplicates
  # This works for gtex too, because the first 15 characters are the sample name
  idcs_nonduplicate <- obj$filterDuplicatesSeqDepth(expression_count_matrix = test_exp_count)
  print(paste('LOG:',"There are",length(idcs_nonduplicate),"non duplicate samples", sep = " ", ""))

  ## We need to remove only duplicates before batch correction

  test_exp_rds = test_exp_rds[, idcs_nonduplicate]
  #test_exp_all = test_exp_all[, idcs_nonduplicate]
  test_exp_count = test_exp_count[, idcs_nonduplicate]
  test_exp_xpm = test_exp_xpm[, idcs_nonduplicate]
  test_exp_logxpm = test_exp_logxpm[, idcs_nonduplicate]

### Batch correction through LIMMA
#### most probablyfor these steps we need to have the indices of samples that survived the purity and mintpt filtering
### also I pretend coad_exp contains the normalized expression
### Here we check that the variable is passed and is in the colData
### If to_batch_correct_nominal_variable is empty, this part of code is skipped
if ((nchar(to_batch_correct_nominal_name)>1)&(to_batch_correct_nominal_name %in% names(colData(test_exp_rds)))){
    ## get the variable to batch correct
    to_batch_correct_nominal_variable <- colData(test_exp_rds)[[to_batch_correct_nominal_name]]

    ## count how many values (if 1, there is nothing to correct)
    num_values <- length(unique(to_batch_correct_nominal_variable))
    if(num_values > 1)
      {
      # First we plot the pca colored by correction variable
        ## Now we plot the first and second component of the PCA
        first_pca_res <- prcomp(t(test_exp_logxpm))
        dtp <- data.frame( first_pca_res$x[,1:2]) # the first two componets are selected (NB: you can also select 3 for 3D plottings or 3+)
        p1 <- ggplot(data = dtp) + 
            geom_point(aes(x = PC1, y = PC2, col = to_batch_correct_nominal_variable))+
            labs(color=to_batch_correct_nominal_name) 

      print(paste('Batch correction for ', to_batch_correct_nominal_name))
      test_exp_logxpm <- removeBatchEffect(test_exp_logxpm, batch = to_batch_correct_nominal_variable) 
      
      # We plot the pca colored by correction variable after the batch correction
        first_pca_res <- prcomp(t(test_exp_logxpm))
        dtp <- data.frame( first_pca_res$x[,1:2]) # the first two componets are selected (NB: you can also select 3 for 3D plottings or 3+)
        p2 <- ggplot(data = dtp) + 
            geom_point(aes(x = PC1, y = PC2, color = to_batch_correct_nominal_variable)) 
        pgrid = plot_grid(p1, p2, labels = c('original', 'corrected'))+
            labs(color=to_batch_correct_nominal_name) 
    }else{
      print('Variable to batch correct for has only one value, so proceeding without batch correcting')
      
      first_pca_res <- prcomp(t(test_exp_logxpm))
      dtp <- data.frame( first_pca_res$x[,1:2]) # the first two componets are selected (NB: you can also select 3 for 3D plottings or 3+)
      p1 <- ggplot(data = dtp) + 
            geom_point(aes(x = PC1, y = PC2, col = to_batch_correct_nominal_variable)) +
            labs(color=to_batch_correct_nominal_name) 
      pgrid = plot_grid(p1, p1, labels = c('original', 'uncorrected'))
    }
}else{ if (to_batch_correct_nominal_name==''){
      print("WARNING: no batch correction variable is passed")
  }else{
    sprintf("WARNING:Provided variable name: '%s' to batch correct for is not part of the rds colData", to_batch_correct_nominal_name)
  }
  first_pca_res <- prcomp(t(test_exp_logxpm))
    dtp <- data.frame( first_pca_res$x[,1:2]) # the first two componets are selected (NB: you can also select 3 for 3D plottings or 3+)
    p1 <- ggplot(data = dtp) + 
            geom_point(aes(x = PC1, y = PC2)) + 
            theme_minimal() 
      pgrid = plot_grid(p1, p1, labels = c('original', 'original'))#), label_size = 12)
  
}
  
  # We save the PCA plots
  ggsave(file=output_pca, plot=pgrid, width = 12,height = 3)

  # Get indices of genes that have a minimum TPM in the data.
  # Here we are filtering by gene
  test_exp_xpm_df <- data.frame(test_exp_xpm)
  idcs_genes_minxpm <- obj$filterGenesByNormExpression(test_exp_xpm_df, min_tpm, frac_samples)
  print(paste('LOG:',"There are",length(idcs_genes_minxpm),"genes that pass mintpm filtering", sep = " ", ""))


  # Select the normalization strategy to be saved
  if (normalization=='logtpm'){
      print('Using logtpm...')
      test_exp_final = test_exp_logxpm
      } else if (normalization=='tpm') {
        print('Using tpm...')
    test_exp_final = test_exp_xpm
  } else if (normalization=='count') {
    print('Using counts...')
    test_exp_final = test_exp_count
  } else if (normalization=='cpm') {
    print('Using cpm...')
    test_exp_final = test_exp_xpm
  } else if (normalization=='logcpm') {
    print('Using logcpm...')
    test_exp_final = test_exp_logxpm
  } else {
    stop("ERROR: normalization unknown")
  } 

if (project_name %in% gtex_projects){

  # get final datasets
  final_table = test_exp_final[idcs_genes_minxpm, ]
  print('Saving table...')
  print(paste('LOG:',"Table dimensions:",dim(final_table)[1],"x",dim(final_table)[2], sep = " ", ""))
  write.table(final_table, file = output_table, sep = '\t', quote = F, col.names=NA)
  print('DONE!')

  final_rds = test_exp_rds[idcs_genes_minxpm,]

  # we add the normalized data as a new assay named like the normalization strategy
  # we need to rename the columns with the sample id as they are in the rds
  renamed_table = final_table
  colnames(renamed_table) = rownames(colData(final_rds))
  # we can now add the normalized data
  assays(final_rds)[[normalization]] = renamed_table

  print('Saving RDS...')
  saveRDS(final_rds, file = output_rds)

} else {
## TCGA samples
idcs_nonduplicate = seq(1:ncol(test_exp_final))
# Filter by purity
# we need a case where we get all tumors, regardless of purity (we include the NAs)
if (th_purity==0){
  print('Not filtering by purity. Pass a threshold>0')
  idcs_purity = idcs_nonduplicate
} else {
   # Get indices of samples passing purity filtering
    if (project_name %in% with_purity){
      print(paste("Purity threshold: ",th_purity, sep = " ", " "))
    idcs_purity <- obj$filterPurity(colnames(test_exp_final), threshold = th_purity, method="CPE")
    } else {
      print(paste("WARNING: no purity for", project_name))
      idcs_purity = idcs_nonduplicate
}}
if (length(idcs_purity)==0){
print(paste('WARNING:',"There are no samples that pass purity filtering (no pure samples:",length(idcs_purity),")", sep = " ", ""))
} else {
   print(paste('LOG:',"There are",length(idcs_purity),"samples that pass purity filtering", sep = " ", ""))
}

print(paste("Tissue type: ",tissue_type, sep = " ", " "))

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
  final_table = test_exp_final[idcs_genes_minxpm, idcs_final]
  print('Saving table...')
  print(paste('LOG:',"Table dimensions:",dim(final_table)[1],"x",dim(final_table)[2], sep = " ", ""))
  write.table(final_table, file = output_table, sep = '\t', quote = F, col.names=NA)
  print('DONE!')

  final_rds = test_exp_rds[idcs_genes_minxpm, idcs_final]

  # we add the normalized data as a new assay named like the normalization strategy
  # we need to rename the columns with the sample id as they are in the rds
  renamed_table = final_table
  colnames(renamed_table) = rownames(colData(final_rds))
  # we can now add the normalized data
  assays(final_rds)[[normalization]] = renamed_table

  print('Saving RDS...')
  saveRDS(final_rds, file = output_rds)

}