library("optparse")
library(WGCNA)
library(doParallel)
registerDoParallel(cores=7)
library(ggplot2)


option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Gene expression table", metavar="input"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Output network from WGCNA", metavar="output"),
 make_option(c("--diagnostic_plot"), type="character", default='wgcna_diagnostic_power_plot.png', 
              help="output filename for the figure of the power", metavar="diagnostic_plot"),
 make_option(c("--power"), type="numeric", default=10, 
              help="WGCNA power", metavar="power")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# gene expression
input_fn = opt$input #args[1]
# output network
output_fn = opt$output #args[2]
# outputs
diagnostic_plot_fn = opt$diagnostic_plot #args[3]
power = opt$power #args[3]

print('Reading input file')
df_xprs <-data.frame(read.csv(input_fn,header = TRUE,
                                row.names = 1,sep = "\t"))
print('Pick soft threshold')
soft_th <- pickSoftThreshold(df_xprs) 
print('Draw diagnostic plot')
# We draw the diagnostic plot
jpeg(file=diagnostic_plot_fn)
par(mar=c(1,1,1,1))
plot(soft_th$fitIndices[,1], soft_th$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="l",
     main = paste("Mean connectivity")) 
text(soft_th$fitIndices[,1], soft_th$fitIndices[,5], labels= soft_th$fitIndices[,1],col="red")
dev.off()
# We get the adjacency matrix
print('Get the adjacency matrix')
adjacency_xprs <- round(adjacency(t(df_xprs),power = power),3)
print(adjacency_xprs[1:5,1:5])

print('Saving adjacency matrix')
write.csv(adjacency_xprs,output_fn, row.names = TRUE)