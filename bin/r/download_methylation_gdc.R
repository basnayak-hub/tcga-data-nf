library(GenomicDataCommons)
library(optparse)
library(NetworkDataCompanion)

# args<-commandArgs(TRUE)

# gdc_cases.project.project_id = args[1]
# gdc_type = args[2]
# gdc_platform = args[3]
# downloads_dir = args[4]
# output_dir = args[5]

option_list = list(
  make_option(c("-p", "--project_id"), type="character", default=NULL, 
              help="project_id", metavar="character"),
  make_option(c("-t", "--type"), type="character", default=NULL, 
              help="gdc_type", metavar="character"),
  make_option(c("--platform"), type="character", default=NULL, 
              help="gdc_platform", metavar="character"),
  make_option(c("-d", "--downloads_dir"), type="character", default=NULL, 
              help="downloads_dir", metavar="character"),
  make_option(c("--manifest_outpath"), type="character", default=NULL, 
              help="", metavar="character"),  
  make_option(c("--pathlist_outpath"), type="character", default=NULL, 
              help="", metavar="character"), 
  make_option(c("--header_outpath"), type="character", default=NULL, 
              help="", metavar="character"),
  make_option(c("--sample_list"), type="character", default=' ', 
              help="test file with sample names", metavar="character")    
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt)

gdc_cases.project.project_id = opt$project_id # args[1]
gdc_type = opt$type # always methylation_beta_value at this point
gdc_platform = opt$platform # which methylation platform
downloads_dir = opt$downloads_dir # where to temporarily cache the GDC data
manifest_outpath = opt$manifest_outpath # where to save the manifest 
pathlist_outpath = opt$pathlist_outpath # where to save the list of sample files
header_outpath = opt$header_outpath # where to save the header for the final methylation data


gdc_set_cache(directory = paste(c(downloads_dir,"tcga_methylation",gdc_cases.project.project_id),collapse="/"))

ge_manifest = files() %>%
  GenomicDataCommons::filter( cases.project.project_id == gdc_cases.project.project_id) %>%
  GenomicDataCommons::filter( type == gdc_type ) %>%
  GenomicDataCommons::filter( platform == gdc_platform ) %>% 
  manifest()
  
dim(ge_manifest)
head(ge_manifest)

# Adding support to filter the sample list
if (nchar(opt$sample_list)>3){
  ong = NetworkDataCompanion::CreateNetworkDataCompanionObject()
  submitters = read.table(opt$sample_list, header = FALSE, sep = "", dec = ".")
  print(submitters$V1)
  # Selecting the submitter IDs
  nnn = nchar(submitters$V1[1])
  ge_manifest = ge_manifest[substr(ong$mapUUIDtoTCGA(ge_manifest$id)$submitter,1,nnn) %in% submitters$V1,]
  print(ge_manifest)
}

# save manifest # TCGA_methylation_manifest.txt"
write.table(ge_manifest,file = manifest_outpath, sep = "\t", row.names = FALSE, quote = FALSE)

# write header column for methylation data
outstring = paste(c("probeID",ge_manifest$id),collapse=" ") # delimiter is " " to match the bash join results
# uncomment the line below for development
# outstring = paste(c("probeID",ge_manifest$id[1:5]),collapse=" ") # delimiter is " " to match the bash join results

write.table(outstring, file = header_outpath,row.names=FALSE,quote=FALSE,col.names=FALSE)

fullpath_list = list()
for(i in 1:nrow(ge_manifest)){
# uncomment the line below for development
#for(i in 1:3){
      options(warn=2)
      print(paste("Processing file:",i))
        
      download_finished <- FALSE

      while(!download_finished) {

          tmp <- tryCatch({

                  fullpath = gdcdata(ge_manifest$id[[i]])
                  print(ge_manifest$id[[i]])
                  print(paste("fullpath:",fullpath))
                  download_finished <- TRUE

                },
                error = function(e) {Sys.sleep(800)})
        }
      
      # This should be the line that actually downloads the file 

      fullpath_list[[i]] = fullpath
      dirname = names(fullpath) 
}

# KS: pathdf is written to a convenience file for the bash join
pathdf = data.frame("file"=unlist(fullpath_list)) 
write.table(pathdf,file = pathlist_outpath, sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
print("DONE: all files downloaded")