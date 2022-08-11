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

resultsDir = args[1]
gdc_cases.project.project_id = args[2]
gdc_type = args[3]
gdc_platform = args[4]

gdc_set_cache(directory = paste(c(resultsDir,"tcga_methylation",gdc_cases.project.project_id),collapse="/"))

ge_manifest = files() %>%
  GenomicDataCommons::filter( cases.project.project_id == gdc_cases.project.project_id) %>%
  GenomicDataCommons::filter( type == gdc_type ) %>%
  GenomicDataCommons::filter( platform == gdc_platform ) %>% 
  manifest()
  
dim(ge_manifest)
head(ge_manifest)

write.table(ge_manifest,file = paste(resultsDir,"TCGA_methylation_manifest.txt",sep="/"), sep = "\t", row.names = FALSE, quote = FALSE)

for(i in 1:nrow(ge_manifest)){

      options(warn=2)
      print(paste("Processing file:",i))
      
      fullpath = gdcdata(ge_manifest$id[[i]])
      print(paste("fullpath:",fullpath))
      dirname = names(fullpath)
      filename = tail(strsplit(fullpath,"/")[[1]],n=1)
      print(paste("Filename:",filename))
      
      put_object(file=fullpath,
	object = filename,
     	bucket = paste(c(s3_bucket,gdc_cases.project.project_id,"/",dirname),collapse=""), 
       	region="us-east-2",
       	multipart=F) # multipart=F because there is a warning message with multipart=T
      	
      # remove file from system

      splitpath = strsplit(fullpath,"/")[[1]]
      folderpath = paste(splitpath[-length(splitpath)],collapse="/")
      system2(command="rm",args=c("-r",folderpath))

}
