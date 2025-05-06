library(netZooR)
library("optparse")
library(tidyverse)


################################################
############ HELPER FUNCTIONS ################
################################################

create_alpaca_table = function(alpaca_obj)
{
# Membership table
membership <-  as.vector(alpaca_obj[[1]])
names(membership) <- names(alpaca_obj[[1]])
membership_tibble <- tibble(node = names(membership), module = membership)

#Modularity scores table
scores <-  as.vector(alpaca_obj[[2]])
names(scores) <- names(alpaca_obj[[2]])

scores_tibble <- tibble(node = names(scores), modularity = scores)

# Join the tables and create a unique table
alpaca_tibble = full_join(membership_tibble, scores_tibble, by = "node")

  return(alpaca_tibble)
}


create_top_genes_table = function(top_genes)
{
# create a table with cluster name and top genes
  top_genes_tibble <- tibble(cluster = top_genes[[2]], genes = top_genes[[1]])

 top_genes_tibble <-  top_genes_tibble %>%
  mutate(genes = map_chr(genes, ~ paste(.x, collapse = ", ")))
  return(top_genes_tibble)
}

################################################
######### OPTPARSE #############################
################################################

option_list = list(
  make_option(c("-b", "--baseline"), type="character", default=NULL, 
              help="first network to be compared", metavar="baseline"),
  make_option(c("-p", "--perturbed"), type="character", default=NULL, 
              help="second network to be compared", metavar="perturbed"),
 make_option(c("-o", "--output_table"), type="character", default='alpaca.csv', 
              help="output filename for table", metavar="output_table"),
make_option(c("-a", "--output_alpaca"), type="character", default='alpaca.Rdata', 
              help="output filename for alpaca object", metavar="output_alpaca"),
 make_option(c("-t", "--top_genes"), type="character", default='alpaca_top_genes.csv', 
              help="output filename for alpaca top genes", metavar="top_genes"),
 make_option(c("--size"), type="numeric", default=10, 
              help="size of output genes", metavar="size")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# network 1
baseline_fn = opt$baseline #args[1]
# network 2
perturbed_fn = opt$perturbed #args[2]
# outputs
output_table_fn = opt$output_table #args[3]
output_alpaca_fn = opt$output_alpaca #args[3]
top_genes_fn = opt$top_genes #args[3]
size = opt$size #args[4]

# Read first network
print(paste('LOG: Reading baseline network:',baseline_fn))
baseline = read_table(baseline_fn)
baseline <- baseline    %>%
    pivot_longer(cols = -tf, names_to = "Gene", values_to = "baseline")

# Read second network
print(paste('LOG: Reading perturbed network:',perturbed_fn))
perturbed = read_table(perturbed_fn)
perturbed <- perturbed    %>%
    pivot_longer(cols = -tf, names_to = "Gene", values_to = "perturbed")

# now we merge the two dataframes based on the tf and Gene columns
print('LOG: Merging networks')
merged_networks <- merge(baseline, perturbed, by = c("tf", "Gene"))

# Run alpaca
print('LOG: Running alpaca')
alpaca_obj <- alpaca(merged_networks,NULL,verbose=F)
# Save the object
save(alpaca_obj, file = output_alpaca_fn)
# Create table and save it
print('LOG: Creating and saving alpaca table')
alpaca_tibble = create_alpaca_table(alpaca_obj) 
alpaca_tibble %>% write_csv(output_table_fn)


# Count nodes ending in _B
count_B <- alpaca_tibble %>%
  filter(grepl("_B$", node)) %>%
  count(module, name = "count_gene")

# extract top genes for each cluster
print('LOG: Extracting top genes')
top_genes <- alpacaExtractTopGenes(alpaca_obj,c(size))
create_top_genes_table(top_genes) %>% write_csv(top_genes_fn)

