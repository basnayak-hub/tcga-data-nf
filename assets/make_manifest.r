# Minimum working example of manifest creation
# This sources the Illumina 450k manifest, annotated to the
# HG38 genome build with gencode v36.
# This is from the public manifest resources here: https://zwdzwd.github.io/InfiniumAnnotation

library(NetSciDataCompanion)
my_friend = NetSciDataCompanion::CreateNetSciDataCompanionObject()

my_manifest = my_friend$mapProbesToGenes(probelist = NULL, # this is the default, and maps every probe in the manifest
                                         rangeUp = 200,
                                         rangeDown = 0)

write.csv(my_manifest,file="450k_manifest_TSS200_TSS0_one_probe_to_many_genes.csv")
