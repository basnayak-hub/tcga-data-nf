library("optparse")
library(GENIE3)

check_k <- function(k) {
  if (k == "sqrt" || k == "all") {
    return(k)  # Valid if it's "sqrt" or "all"
  } else if (is.numeric(as.numeric(k)) && k == as.integer(k) && k >= 1) {
    return(as.numeric(k))  # Valid if it's an integer between 1 and p
  } else {
    stop('k needs to be either sqrt, all, or an integer') # Invalid otherwise
  }
}

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Gene expression table", metavar="input"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output network from WGCNA", metavar="output"),
 make_option(c("--tf_list"), type="character", default=" ", 
              help="TF list filename. Pass a text file with TF names to filter the output", metavar="tf_list"),
 make_option(c("--k"), type="character", default="sqrt", 
              help="number of candidate regulators", metavar="k"),
 make_option(c("--tree_method"), type="character", default="RF", help="tree method", metavar="tree_method"),
 make_option(c("--n_trees"), type="numeric", default=1000, 
              help="number of trees", metavar="n_trees"),
 make_option(c("--n_cores"), type="numeric", default=1, 
              help="number of cores", metavar="n_cores")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# gene expression
input_fn = opt$input #args[1]
# output network
output_fn = opt$output #args[2]
k = opt$k #args[3]
tree_method = opt$tree_method #args[3]
n_trees = opt$n_trees #args[3]
n_cores = opt$n_cores #args[3]
# These: if none, keep everything
tf_list_fn =  opt$tf_list #args[5]


print('Reading input file')
df <-data.frame(read.csv(input_fn,header = TRUE,
                                row.names = 1,sep = "\t"))

print(class(df))
print('Convert to matrix')
# Step 2: Convert the data frame to a matrix
expr_mat <- as.matrix(df)

# Step 3: Assign row and column names if they aren't preserved during conversion
# (Usually, row names and column names are preserved, but it's good to double-check)
rownames(expr_mat) <- rownames(df)
colnames(expr_mat) <- colnames(df)

print("Check parameter")
k = check_k(k)
print(paste('K:',as.character(k)))

set.seed(123) # For reproducibility of results
# Manage list of TF, if list is passed the TF will be given as regulators
if (tf_list_fn!=' '){
  # Read TF list
  print(paste(c('Reading TF list', tf_list_fn)), collapse = " ")
  tf_list = read.table(tf_list_fn)[,1] #read.table(paste(baseDir,"ext/TF_names_v_1.01.txt",sep="/"))[,1]
  print(head(tf_list))

  tf_list_filtered = intersect(tf_list,rownames(expr_mat))

  # check that the TFs are in the matrix
    if (length(tf_list_filtered)==0){
        stop('No TFs found in the matrix')
    } else {
       print(paste(c('TFs found in the matrix:', length(tf_list_filtered)), collapse = " "))
    }


} else {
   tf_list_filtered = NULL
}

weightMat <- GENIE3(expr_mat, 
            regulators=tf_list_filtered, 
            treeMethod=tree_method, 
            K=k, 
            nTrees=n_trees,
            nCores=n_cores, 
            verbose = TRUE)

# Save the weight matrix
write.table(weightMat, file=output_fn, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)





