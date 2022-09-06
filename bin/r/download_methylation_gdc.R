# 20220620
# Pulls methylation data from GDC and writes to S3
# Format of the GDC data is one file per participant, two columsn: probe\tbeta
# To get the participant ID, you have to look at the file name which has a UUID
# You can then map the UUID back to the TCGA id

# 20220811 Placeholder file for new repo.

# With VF: you will find suggestions from viola. 
# Follow as you see fit
# KS: notes re VF

# VF: you won't need to explicitly save things on s3. Remove 'put'
# KS: restructuring: this iterates through all of the files and saves them one at a time
# on s3 rather than outputting a single file. I think we need to save them
# locally within the nextflow filesystem somewhere, then use bash for
# a join, then output the joined file to s3. 

library(GenomicDataCommons)

args<-commandArgs(TRUE)

gdc_cases.project.project_id = args[1]
gdc_type = args[2]
gdc_platform = args[3]
downloads_dir = args[4]
output_dir = args[5]

# Make output dir if it doesn't exist
# KS @VF: how do we avoid this?

#if(!file.exists(output_dir))
#   dir.create(output_dir)

gdc_set_cache(directory = paste(c(downloads_dir,"tcga_methylation",gdc_cases.project.project_id),collapse="/"))

ge_manifest = files() %>%
  GenomicDataCommons::filter( cases.project.project_id == gdc_cases.project.project_id) %>%
  GenomicDataCommons::filter( type == gdc_type ) %>%
  GenomicDataCommons::filter( platform == gdc_platform ) %>% 
  manifest()
  
dim(ge_manifest)
head(ge_manifest)

# write.table(apply(ge_manifest,1,paste,sep="/"),file = paste(output_dir,"TCGA_methylation_manifest.txt",sep="/"), sep = "\t", row.names = FALSE, quote = FALSE)

fullpath_list = list()
for(i in 1:nrow(ge_manifest)){
      options(warn=2)
      print(paste("Processing file:",i))
      
      # This should be the line that actually downloads the file 
      fullpath = gdcdata(ge_manifest$id[[i]])
      print(paste("fullpath:",fullpath))
      fullpath_list[[i]] = fullpath
      dirname = names(fullpath) 
}

# KS: pathdf is written to a convenience file for the bash join
# We could probably do it directly from the manifest but
# this seemed more straightforward

pathdf = data.frame("file"=unlist(fullpath_list))
write.table(pathdf,file = "TCGA_methylation_paths.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
