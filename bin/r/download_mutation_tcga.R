library(TCGAbiolinks)
library(tidyverse)
args<-commandArgs(TRUE)

# args[1] = project for tcga
# args[2] = data category, should be "Simple Nucleotide Variation"
# args[3] = data type, should be  "Masked Somatic Mutation"
# args[4] = directory, where to download the data
# args[5] = out_maf, save the table with complete mutation data
# args[6] = out_pivot, save the pivot table with gene name, sample name, and number of mutations

project = args[1]
data_category = args[2]
data_type = args[3]
directory = args[4]
out_table = args[5]
out_pivot = args[6]

print('Arguments:')
print(args)
print('Downloading...')
## Download data and create the RSE object
query1 <- GDCquery( project = project, data.category = data_category, data.type = data_type, legacy=F)
GDCdownload(query1, directory = directory)
print('Preparing the mutation data table...')
muts <- GDCprepare(query1, directory = directory)
print(paste0('Saving the mutation data table to: ',out_table))
write.table(muts, out_table)

# muts is a table where each row is a mutation entry. For each mutation the gene it hits
# the location on the genome are reported. Also both the tumor and matched normal sample
# are reported.
muts$mutations=1
print('Preparing the pivot table...')
all_muts= muts %>% 
            group_by(Entrez_Gene_Id,Hugo_Symbol, Tumor_Sample_Barcode) %>%
            summarise(abundance = sum(mutations)) %>%
            #group_by(Entrez_Gene_Id) %>% 
            select(Entrez_Gene_Id,Hugo_Symbol, Tumor_Sample_Barcode,abundance) %>%
            #mutate(row = row_number()) %>%
            pivot_wider(names_from = Tumor_Sample_Barcode, values_from = abundance, values_fill = 0)

print(paste0('Saving the mutation pivot table to: ',out_pivot))
write_csv(all_muts, out_pivot)