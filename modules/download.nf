process downloadRecount{

    publishDir "${params.resultsDir}/${uuid}/recount3/", pattern: "${uuid}.rds", mode: 'copy', overwrite: true
    
    input:
        tuple val(uuid),val(project),val(project_home),val(organism),val(annotation),val(type), val(samples)
    output:
        tuple val(uuid),val(project),val(project_home),val(organism),val(annotation),val(type), val(samples) ,file("${uuid}.rds"),file("${uuid}_recount3_metatada.csv")
    script:
        """
            Rscript '${baseDir}/bin/r/download_expression_recount.R' ${project} ${project_home} ${organism} ${annotation} ${type} ${samples} "${uuid}.rds" > "${uuid}_recount3_downloads.log";
            echo "uuid,project,project_home,organism,annotation,type,samples,output_rds" > "${uuid}_recount3_metatada.csv";
            echo "${uuid},${project},${project_home},${organism},${annotation},${type},${samples},${params.resultsDir}/${uuid}/recount3/${uuid}.rds" >> "${uuid}_recount3_metatada.csv"
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
        tuple val(uuid),val(project),val(data_category),val(data_type),val(download_dir), val(samples)
    output:
        tuple val(uuid),val(project),val(data_category),val(data_type),val(download_dir), val(samples),file("${uuid}_mutations.txt"),file("${uuid}_mutations_pivot.csv"),file("${uuid}_mutations_metadata.csv")
    script:
        """
            Rscript '${baseDir}/bin/r/download_mutation_tcga.R' ${project}  "${data_category}" "${data_type}" "${download_dir}" "${samples}" "${uuid}_mutations.txt" "${uuid}_mutations_pivot.csv";
            echo "uuid,project,data_category,data_type,download_dir,samples,mutation_table,pivot_table" > "${uuid}_mutations_metadata.csv";
            echo "${uuid},${project},${data_category},${data_type},${download_dir},${samples},${params.resultsDir}/${uuid}/mutations/${uuid}_mutations.txt,${params.resultsDir}/${uuid}/mutations/${uuid}_mutations_pivot.csv" >> "${uuid}_mutations_metadata.csv"
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
        tuple val(uuid),val(project),val(gdc_type),val(gdc_platform),val(download_dir),val(samples)
    output:
        tuple val(uuid),val(project),val(gdc_type),val(gdc_platform),val(download_dir),val(samples),file("${uuid}_methylation_manifest.txt"), file("${uuid}_methylations.txt"), file("${uuid}_methylation_metadata.csv")
    script:
        """
            Rscript '${baseDir}/bin/r/download_methylation_gdc.R' -p '${project}'  -t '${gdc_type}' --platform '${gdc_platform}' -d '${download_dir}' --manifest_outpath '${uuid}_methylation_manifest.txt' --pathlist_outpath '${uuid}_methylation_paths.txt' --header_outpath '${uuid}_methylation_header.txt' --sample_list ${samples}
            bash '${baseDir}/bin/bash/join_methylation_gdc.sh'  "${uuid}_methylations.txt" "${uuid}_methylation_paths.txt"
	        cat  '${uuid}_methylation_header.txt' "${uuid}_methylations.txt" > "${uuid}_methylations_labeled.txt"
            mv "${uuid}_methylations_labeled.txt" "${uuid}_methylations.txt";
            echo "uuid,project,gdc_type,gdc_platform,download_dir,samples,methylation_manifest,methylation_table" > "${uuid}_methylation_metadata.csv";
            echo "${uuid},${project},${gdc_type},${gdc_platform},${download_dir},${samples},${params.resultsDir}/${uuid}/methylation/${uuid}_methylation_manifest.txt,${params.resultsDir}/${uuid}/methylation/${uuid}_methylations.txt" >> "${uuid}_methylation_metadata.csv"
        
        """

    stub:
        """
        touch "${uuid}_methylation_manifest.txt"
        touch "${uuid}_methylations.txt"
        touch "${uuid}_methylation_metadata.csv"
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

workflow downloadRecount3Wf{
    take: channelRecount//parData
    main:
        // channelRecount = Channel.from(parData.entrySet())
        //                                 .map{
        //                                     item -> tuple(
        //                                         item.getKey(),
        //                                         item.value.project,
        //                                         item.value.project_home,
        //                                         item.value.organism,
        //                                         item.value.annotation,
        //                                         item.value.type,
        //                                         item.value.samples,
        //                                     )
        //                                 }.view()
            dr = downloadRecount(channelRecount)
            mergeRecountMetadata(dr.map{it -> it[-1]}.collect())
    emit: dr
}

workflow downloadMethylationWf{
    take: channelMethylation
    main:
        // channelMethylation = Channel.from(parData.entrySet())
        //                             .map{
        //                                 item -> tuple(
        //                                     item.getKey(),
        //                                     item.value.project,
        //                                     item.value.gdc_type,
        //                                     item.value.gdc_platform,
        //                                     item.value.download_dir,
        //                                     item.value.samples,
        //                                 )
        //                             }.view()

        dme = downloadMethylation(channelMethylation)
        mergeMethylationMetadata(dme.map{it -> it[-1]}.collect())
    emit: dme
}

workflow downloadMutationsWf{
    take: channelMutation
    main:
        // channelMutation = Channel.from(parData.entrySet())
        //                             .map{
        //                                 item -> tuple(
        //                                     item.getKey(),
        //                                     item.value.project,
        //                                     item.value.data_category,
        //                                     item.value.data_type,
        //                                     item.value.download_dir,
        //                                     item.value.samples,
        //                                 )
        //                             }.view()

            dmu = downloadMutations(channelMutation)
            mergeMutationsMetadata(dmu.map{it -> it[-1]}.collect())
    emit: dmu
}

workflow downloadClinicalWf{
    take: channelClinical
    main:
        // channelClinical = Channel.from(parData.entrySet())
        //                             .map{
        //                                 item -> tuple(
        //                                     item.getKey(),
        //                                     item.value.project,
        //                                     item.value.data_category,
        //                                     item.value.data_type,
        //                                     item.value.data_format,
        //                                 )
        //                             }.view()
        dcli = downloadClinical(channelClinical)
    emit: dcli
}



workflow downloadWf{
    main:        
    
        // Read the metadata file
        expression_recount3 = Channel.fromPath(params.download_metadata).splitJson(path: "expression_recount3")
        mutation_tcgabiolinks = Channel.fromPath(params.download_metadata).splitJson(path: "mutation_tcgabiolinks")
        clinical_tcgabiolinks = Channel.fromPath(params.download_metadata).splitJson(path: "clinical_tcgabiolinks")
        methylation_gdc = Channel.fromPath(params.download_metadata).splitJson(path: "methylation_gdc")

        // dataCh.branch {
        // expression_recount3: (it.key == 'expression_recount3')
        //     return it.value.entrySet()
        // mutation_tcgabiolinks: (it.key=='mutation_tcgabiolinks')
        //     return it.value.entrySet()
        // clinical_tcgabiolinks: (it.key=='clinical_tcgabiolinks')
        //     return it.value.entrySet()
        // methylation_gdc: (it.key=='methylation_gdc')
        //     return it.value.entrySet()
        // other: true
        //     return it.value.entrySet()
        // }.set{branchCh}




        // Process recount3 data
        if (expression_recount3!=[null]){
            channelRecount = expression_recount3.map{
                                            item -> tuple(
                                                item.key,
                                                item.value.project,
                                                item.value.project_home,
                                                item.value.organism,
                                                item.value.annotation,
                                                item.value.type,
                                                item.value.samples,
                                            )
                                        }.view()
        
            downloadRecount3Wf(channelRecount)
        }

        //process mutations
        if (mutation_tcgabiolinks!=[null]){
            channelMutation = mutation_tcgabiolinks.map{
                                            item -> tuple(
                                                item.key,
                                                item.value.project,
                                                item.value.data_category,
                                                item.value.data_type,
                                                item.value.download_dir,
                                                item.value.samples
                                            )
                                        }.view()
        
            downloadMutationsWf(channelMutation)
        } 
        // Process methylation
        if (methylation_gdc!=[null]){
            channelMethylation = methylation_gdc.map{
                                            item -> tuple(
                                                item.key,
                                                item.value.project,
                                                item.value.gdc_type,
                                                item.value.gdc_platform,
                                                item.value.download_dir,
                                                item.value.samples,
                                            )
                                        }.view()
            downloadMethylationWf(channelMethylation)
        }
        // Process clinical
        if (clinical_tcgabiolinks!=[null]){
            channelClinical = clinical_tcgabiolinks.map{
                                            item -> tuple(
                                                item.key,
                                                item.value.project,
                                                item.value.data_category,
                                                item.value.data_type,
                                                item.value.data_format,
                                            )
                                        }.view()
            downloadClinicalWf(channelClinical)
    }
}


workflow fullDownloadWf{
    take:
        dataCh
    main:
        //println(dataCh.value.getClass())
        // Empty channels
        dr = Channel.empty()
        dmu = Channel.empty()
        dme = Channel.empty()
        dc = Channel.empty()
        // DOWNLOAD RECOUNT3

        modalities = dataCh.map{it -> it.value.keySet()}.collect().view()
        //if (dataCh.containsKey('expression_recount3')){
            dChRe = dataCh.map{it -> tuple(
                                                it.key,
                                                it.value.get('expression_recount3').project,
                                                it.value.get('expression_recount3').project_home,
                                                it.value.get('expression_recount3').organism,
                                                it.value.get('expression_recount3').annotation,
                                                it.value.get('expression_recount3').type,
                                                it.value.get('expression_recount3').samples)
                                                }.view()
            //dChRe = Channel.from(dataCh.key).combine(Channel.from(dataCh.value.expression_recount3.project))
            //                        .combine(Channel.from(dataCh.value.expression_recount3.project_home))
            //                        .combine(Channel.from(dataCh.value.expression_recount3.organism))
            //                        .combine(Channel.from(dataCh.value.expression_recount3.annotation))
            //                        .combine(Channel.from(dataCh.value.expression_recount3.type))
            //                        .combine(Channel.from(dataCh.value.expression_recount3.samples)).view()
            // dr is the download recount channel
            // output: tuple val(uuid),val(project),val(project_home),val(organism),val(annotation),val(type),val(samples),file("${uuid}.rds"),file("${uuid}_recount3_metatada.csv")
            dr = downloadRecount3Wf(dChRe)
        //}
        // DOWNLOAD MUTATIONS
        //if (dataCh.map{it -> it.value.keySet()}.contains('mutation_tcgabiolinks')){
            dChMu = dataCh.map{it -> tuple(
                                                it.key,
                                                it.value.get('mutation_tcgabiolinks').project,
                                                it.value.get('mutation_tcgabiolinks').data_category,
                                                it.value.get('mutation_tcgabiolinks').data_type,
                                                it.value.get('mutation_tcgabiolinks').download_dir,
                                                it.value.get('mutation_tcgabiolinks').samples)
                                                }.view()
            //dChMu = Channel.from(dataCh.key).combine(Channel.from(dataCh.value.mutation_tcgabiolinks.project))
            //                        .combine(Channel.from(dataCh.value.mutation_tcgabiolinks.data_category))
            //                        .combine(Channel.from(dataCh.value.mutation_tcgabiolinks.data_type))
            //                        .combine(Channel.from(dataCh.value.mutation_tcgabiolinks.download_dir))
            //                        .combine(Channel.from(dataCh.value.mutation_tcgabiolinks.samples)).view()
            // dr is the download recount channel
            // output: tuple val(uuid),val(project),val(data_category),val(data_type),val(download_dir),file("${uuid}_mutations.txt"),file("${uuid}_mutations_pivot.csv"),file("${uuid}_mutations_metadata.csv")
            dmu = downloadMutationsWf(dChMu)
        //}

        // DOWNLOAD METHYLATION
        //if (dataCh.map{it -> it.value.keySet()}.contains('methylation_gdc')){
            dChMe = dataCh.map{it -> tuple(
                                                it.key,
                                                it.value.get('methylation_gdc').project,
                                                it.value.get('methylation_gdc').gdc_type,
                                                it.value.get('methylation_gdc').gdc_platform,
                                                it.value.get('methylation_gdc').download_dir,
                                                it.value.get('methylation_gdc').samples)
                                                }.view()
        //     //dChMe = Channel.from(dataCh.key).combine(Channel.from(dataCh.value.methylation_gdc.project))
        //     //                        .combine(Channel.from(dataCh.value.methylation_gdc.gdc_type))
        //     //                        .combine(Channel.from(dataCh.value.methylation_gdc.gdc_platform))
        //     //                        .combine(Channel.from(dataCh.value.methylation_gdc.download_dir))
        //     //                        .combine(Channel.from(dataCh.value.methylation_gdc.samples)).view()
        //     // dr is the download recount channel
        //     // output: tuple val(uuid),val(project),val(data_category),val(data_type),val(download_dir),file("${uuid}_mutations.txt"),file("${uuid}_mutations_pivot.csv"),file("${uuid}_mutations_metadata.csv")
        dme = downloadMethylationWf(dChMe)
        //}

        // // DOWNLOAD CLINICAL
        // if (dataCh.map{it -> it.value.keySet()}.contains('clinical_tcgabiolinks')){
        //     dChCl = dataCh.map{it -> tuple(
        //                                         it.key,
        //                                         it.value.get('clinical_tcgabiolinks').project,
        //                                         it.value.get('clinical_tcgabiolinks').data_category,
        //                                         it.value.get('clinical_tcgabiolinks').data_type,
        //                                         it.value.get('clinical_tcgabiolinks').data_format)
        //                                         }.view()
        //     //dChCl = Channel.from(dataCh.key).combine(Channel.from(dataCh.value.clinical_tcgabiolinks.project))
        //     //                        .combine(Channel.from(dataCh.value.clinical_tcgabiolinks.data_category))
        //     //                        .combine(Channel.from(dataCh.value.clinical_tcgabiolinks.data_type))
        //     //                        .combine(Channel.from(dataCh.value.clinical_tcgabiolinks.data_format)).view()
        //     // dr is the download recount channel
        //     // output: tuple val(uuid),val(project),val(data_category),val(data_type),val(download_dir),val(samples),file("${uuid}_mutations.txt"),file("${uuid}_mutations_pivot.csv"),file("${uuid}_mutations_metadata.csv")
        //     dc = downloadClinicalWf(dChCl)
        //}
    emit:
        dr
        dmu
        dme
        dc
    }


