library(data.table)
library(stringr)
library("recount3")
library(NetworkDataCompanion)
library("optparse")
library(dplyr)
library(ggplot2)

plotPCA = function(df, output_fn)
{
  print("Starting PCA diagnostic plot")

  # Ensure all columns are numeric
  df <- df[, sapply(df, is.numeric)]
  
  print("Preview of numeric data:")
  print(df[1:5, 1:min(5, ncol(df))])

  # Remove constant (zero variance) columns
  df <- df[, apply(df, 2, var, na.rm = TRUE) != 0]

  # Remove rows with any NAs
  df_clean <- na.omit(df)

  # Check if any NAs remain
  print(paste("Remaining NAs:", sum(is.na(df_clean))))
  
  # Check if enough data remains
  if (nrow(df_clean) < 2 || ncol(df_clean) < 2) {
    stop("Not enough data to perform PCA.")
  }

  # Perform PCA (scale the data)
  pca <- prcomp(df_clean, scale. = TRUE)

  # Get PCA results
    pca_df <- as.data.frame(pca$x)
    pca_df$Sample <- rownames(pca_df)

    # Plot PCA
    p <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample)) +
      geom_point(color = "steelblue", size = 3) +
      #geom_text(vjust = 1.5, hjust = 1.2, size = 3) +
      xlab(paste0("PC1 (", round(summary(pca)$importance[2,1] * 100, 1), "%)")) +
      ylab(paste0("PC2 (", round(summary(pca)$importance[2,2] * 100, 1), "%)")) +
      theme_minimal()

    # Save the plot
    ggsave(output_fn, plot = p, width = 8, height = 6)
  }





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
              help = "Filter duplicate or missing values [default: %default]"),
make_option(c("--diagnostic_pca"), type = "character", default = NULL,
              help = "Optional output filename to save the figure (e.g., 'plot.png'). If not provided, figure is not saved.")
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
diagnostic_pca = opt$diagnostic_pca #args[8]

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

# Choose name of diagnostic PCA output
if (!is.null(diagnostic_pca)) {
  print("PCA diagnostic plot")
  print(class(meth_df))
  
  result <- tryCatch({
    plotPCA(meth_df, diagnostic_pca)
  }, error = function(e) {
    message("plotPCA failed: ", e$message)
    return(NULL)
  })
  
  if (is.null(result)) {
    message("No PCA plot generated.")
  }
  
} else {
  message("No PCA diagnostic.")
}

# filter methylation beta to only probes that are in the tf list
# Manage list of TF
if (tf_list_fn!=' '){

  # Read TF list
  print(paste(c('Reading TF list', tf_list_fn)), collapse = " ")
  tf_list = read.table(tf_list_fn)[,1] 
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
# Ensure gene_map is a data frame
gene_map <- as.data.frame(gene_map)

# Add UUID as a column (row names of gene_map)
gene_map$UUID <- row.names(gene_map)
# Add UUID as rownames to barcodes, matching gene_map rownames
barcodes_df <- barcodes
barcodes_df$UUID <- row.names(gene_map)  # Only valid if row order matches

# Now do the join
final_df <- gene_map %>%
  mutate(UUID = row.names(gene_map)) %>%
  inner_join(barcodes_df[, c("UUID", "submitter_id")], by = "UUID") %>%
  rename(TCGA_barcode = submitter_id) %>%
  dplyr::select(-UUID)

write.csv(final_df, output_fn, row.names = TRUE, quote = FALSE)
print("Done with get gene level methylation data")