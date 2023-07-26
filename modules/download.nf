process downloadRecount{

    publishDir "${params.resultsDir}/${uuid}/recount3/", pattern: "${uuid}.rds", mode: 'copy', overwrite: true
    
    input:
        tuple val(uuid),val(project),val(project_home),val(organism),val(annotation),val(type)
    output:
        tuple val(uuid),val(project),val(project_home),val(organism),val(annotation),val(type),file("${uuid}.rds"),file("${uuid}_recount3_metatada.csv")
    script:
        """
            Rscript '${baseDir}/bin/r/download_expression_recount.R' ${project} ${project_home} ${organism} ${annotation} ${type} "${uuid}.rds";
            echo "uuid,project,project_home,organism,annotation,type,output_rds" > "${uuid}_recount3_metatada.csv";
            echo "${uuid},${project},${project_home},${organism},${annotation},${type},${params.resultsDir}/${uuid}/recount3/${uuid}.rds" >> "${uuid}_recount3_metatada.csv"
        """
    stub:
        """
        touch "${uuid}.rds"
        touch "${uuid}_recount3_metatada.csv"
        """
}


process mergeRecountMetadata{

    publishDir "${params.resultsDir}/", pattern: "downloaded_recount_metadata.csv", mode: 'copy', overwrite: true
    
    input:
        val(results)
    output:
        tuple val(results), path("downloaded_recount_metadata.csv")

    script:
	values = results.flatten().collect{"$it"}.join(',')
    """
    fit.py merge-tables -o downloaded_recount_metadata.csv -t ${values}
    """

}


process downloadMutations{

    publishDir "${params.resultsDir}/${uuid}/mutations/", pattern: "${uuid}_mutations*", mode: 'copy', overwrite: true
    
    input:
        tuple val(uuid),val(project),val(data_category),val(data_type),val(download_dir)
    output:
        tuple val(uuid),val(project),val(data_category),val(data_type),val(download_dir),file("${uuid}_mutations.txt"),file("${uuid}_mutations_pivot.csv"),file("${uuid}_mutations_metadata.csv")
    script:
        """
            Rscript '${baseDir}/bin/r/download_mutation_tcga.R' ${project}  "${data_category}" "${data_type}" "${download_dir}" "${uuid}_mutations.txt" "${uuid}_mutations_pivot.csv";
            echo "uuid,project,data_category,data_type,download_dir,mutation_table,pivot_table" > "${uuid}_mutations_metadata.csv";
            echo "${uuid},${project},${data_category},${data_type},${download_dir},${params.resultsDir}/${uuid}/mutations/${uuid}_mutations.txt,${params.resultsDir}/${uuid}/mutations/${uuid}_mutations_pivot.csv" >> "${uuid}_mutations_metadata.csv"
        """
    stub:
        """
        touch "${uuid}_mutations.txt"
        touch "${uuid}_mutations_pivot.csv"
        touch "${uuid}_mutations_metadata.csv"
        """

} 

process mergeMutationsMetadata{

    publishDir "${params.resultsDir}/", pattern: "downloaded_mutation_metadata.csv", mode: 'copy', overwrite: true
    
    input:
        val(results)
    output:
        tuple val(results), path("downloaded_mutation_metadata.csv")

    script:
	values = results.flatten().collect{"$it"}.join(',')
    """
    fit.py merge-tables -o downloaded_mutation_metadata.csv -t ${values}
    """

}


process downloadMethylation{
    publishDir "${params.resultsDir}/${uuid}/methylation/", pattern: "${uuid}_methylation*", mode: 'copy', overwrite: true
    
    input:
        tuple val(uuid),val(project),val(gdc_type),val(gdc_platform),val(download_dir)
    output:
        tuple val(uuid),val(project),val(gdc_type),val(gdc_platform),val(download_dir),file("${uuid}_methylation_manifest.txt"), file("${uuid}_methylations.txt"), file("${uuid}_methylation_metadata.csv")
    script:
        """
            Rscript '${baseDir}/bin/r/download_methylation_gdc.R' -p '${project}'  -t '${gdc_type}' --platform '${gdc_platform}' -d '${download_dir}' --manifest_outpath '${uuid}_methylation_manifest.txt' --pathlist_outpath '${uuid}_methylation_paths.txt' --header_outpath '${uuid}_methylation_header.txt'
            bash '${baseDir}/bin/bash/join_methylation_gdc.sh'  "${uuid}_methylations.txt" "${uuid}_methylation_paths.txt"
	        cat  '${uuid}_methylation_header.txt' "${uuid}_methylations.txt" > "${uuid}_methylations_labeled.txt"
            mv "${uuid}_methylations_labeled.txt" "${uuid}_methylations.txt";
            echo "uuid,project,gdc_type,gdc_platform,download_dir,methylation_manifest,methylation_table" > "${uuid}_methylation_metadata.csv";
            echo "${uuid},${project},${gdc_type},${gdc_platform},${download_dir},${params.resultsDir}/${uuid}/methylation/${uuid}_methylation_manifest.txt,${params.resultsDir}/${uuid}/methylation/${uuid}_methylations.txt" >> "${uuid}_methylation_metadata.csv"
        
        """

    stub:
        """
        touch "${uuid}_methylation_manifest.txt"
        touch "${uuid}_methylations.txt"
        touch "${uuid}_methylation_metadata"
        """


}

process mergeMethylationMetadata{

    publishDir "${params.resultsDir}/", pattern: "downloaded_methylation_metadata.csv", mode: 'copy', overwrite: true
    
    input:
        val(results)
    output:
        tuple val(results), path("downloaded_methylation_metadata.csv")

    script:
	values = results.flatten().collect{"$it"}.join(',')
    """
    fit.py merge-tables -o downloaded_methylation_metadata.csv -t ${values}
    """

}


process downloadClinical{
    publishDir "${params.resultsDir}/${uuid}/clinical", pattern: "*.csv", mode: 'copy', overwrite: true
    
    input:
        tuple val(uuid),val(project),val(data_category),val(data_type),val(data_format)
    output:
        tuple val(uuid),val(project),val(data_category),val(data_type),val(data_format),path("*")
    script:
        """
            Rscript '${baseDir}/bin/r/download_clinical_tcga.R' ${project}  "${data_category}" "${data_type}" "${data_format}" "."
        """

    stub:
        """
        touch clinical.csv
        """
}


workflow downloadWf{
    main:
        modalities = params.download_metadata.keySet()

         if (modalities.contains('expression_recount3')){
            // Channel for recount3 data download
            channelRecount = Channel.from(params.download_metadata.expression_recount3.entrySet())
                                        .map{
                                            item -> tuple(
                                                item.getKey(),
                                                item.getValue().project,
                                                item.getValue().project_home,
                                                item.getValue().organism,
                                                item.getValue().annotation,
                                                item.getValue().type
                                            )
                                        }.view()
            
            dr = downloadRecount(channelRecount)
            mergeRecountMetadata(dr.map{it -> it[-1]}.collect())
            
        }

        if (modalities.contains('mutation_tcgabiolinks')){
        // Add here channels for mutations and methylation
        channelMutation = Channel.from(params.download_metadata.mutation_tcgabiolinks.entrySet())
                                    .map{
                                        item -> tuple(
                                            item.getKey(),
                                            item.getValue().project,
                                            item.getValue().data_category,
                                            item.getValue().data_type,
                                            item.getValue().download_dir,
                                        )
                                    }.view()

            dmu = downloadMutations(channelMutation)
            mergeMutationsMetadata(dmu.map{it -> it[-1]}.collect())
        } 
        if (modalities.contains('methylation_gdc')){
        channelMethylation = Channel.from(params.download_metadata.methylation_gdc.entrySet())
                                    .map{
                                        item -> tuple(
                                            item.getKey(),
                                            item.getValue().project,
                                            item.getValue().gdc_type,
                                            item.getValue().gdc_platform,
                                            item.getValue().download_dir,
                                        )
                                    }.view()

            dme = downloadMethylation(channelMethylation)
            mergeMethylationMetadata(dme.map{it -> it[-1]}.collect())
        }

        if (modalities.contains('clinical_tcgabiolinks')){
        channelClinical = Channel.from(params.download_metadata.clinical_tcgabiolinks.entrySet())
                                    .map{
                                        item -> tuple(
                                            item.getKey(),
                                            item.getValue().project,
                                            item.getValue().data_category,
                                            item.getValue().data_type,
                                            item.getValue().data_format,
                                        )
                                    }.view()

            dcli = downloadClinical(channelClinical)

    }
}
