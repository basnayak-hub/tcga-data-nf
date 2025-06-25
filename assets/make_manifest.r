# Minimum working example of manifest creation
# This sources the Illumina 450k manifest, annotated to the
# HG38 genome build with gencode v36.
# This is from the public manifest resources here: https://zwdzwd.github.io/InfiniumAnnotation

library(NetworkDataCompanion)
my_friend = NetworkDataCompanion::CreateNetworkDataCompanionObject()

my_manifest = my_friend$mapProbesToGenes(probelist = NULL, # this is the default, and maps every probe in the manifest
                                         rangeUp = 200,
                                         rangeDown = 0,
                                         localManifestPath="/home/ubuntu/tcga-ov-subtype/data/config/HM27.hg38.manifest.tsv")

write.csv(my_manifest,file="27k_manifest_TSS200_TSS0_one_probe_to_many_genes.csv")
