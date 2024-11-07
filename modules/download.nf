process downloadRecount{

    label "r_download"
    publishDir "${params.resultsDir}/${params.batchName}/${uuid}/data_download/recount3/", pattern: "${uuid}.rds", mode: 'copy', overwrite: true
    
    input:
        tuple val(uuid),val(project),val(project_home),val(organism),val(annotation),val(type), path(samples)
    output:
        tuple val(uuid),val(project),val(project_home),val(organism),val(annotation),val(type), path(samples) ,file("${uuid}.rds"),file("${uuid}_recount3_metatada.csv")
    script:
        log.info "Downloading recount for $uuid, $project"
        """
            if ! test -f ${samples}; then
            touch ${samples}
            fi
            Rscript '${baseDir}/bin/r/download_expression_recount.R' ${project} ${project_home} ${organism} ${annotation} ${type} ${samples} "${uuid}.rds" > "${uuid}_recount3_downloads.log";
            echo "uuid,project,project_home,organism,annotation,type,samples,output_rds" > "${uuid}_recount3_metatada.csv";
            echo "${uuid},${project},${project_home},${organism},${annotation},${type},${samples},${params.resultsDir}/${params.batchName}/${uuid}/data_download/recount3/${uuid}.rds" >> "${uuid}_recount3_metatada.csv"


        """
    stub:
        """
        touch "${uuid}.rds"
        touch "${uuid}_recount3_metatada.csv"
        """
}

// Concatenate all the metadata recount table in one file
process mergeRecountMetadata{

    //conda 'containers/conda_envs/merge_tables.yml'
    label 'merge_tables'

    publishDir "${params.resultsDir}/${params.batchName}/", pattern: "downloaded_recount_metadata.csv", mode: 'copy', overwrite: true
    
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
    // Conda environment is managed by label
    label "r_download"

    publishDir "${params.resultsDir}/${params.batchName}/${uuid}/data_download/mutations/", pattern: "${uuid}_mutations*", mode: 'copy', overwrite: true
    
    input:
        tuple val(uuid),val(project),val(data_category),val(data_type),val(download_dir), path(samples)
    output:
        tuple val(uuid),val(project),val(data_category),val(data_type),val(download_dir), path(samples),file("${uuid}_mutations.txt"),file("${uuid}_mutations_pivot.csv"),file("${uuid}_mutations_metadata.csv")
    script:
    log.info "Downloading mutations for $uuid, $project"
        """
            if ! test -f ${samples}; then
            touch ${samples}
            fi
            Rscript '${baseDir}/bin/r/download_mutation_tcga.R' ${project}  "${data_category}" "${data_type}" "${download_dir}" "${samples}" "${uuid}_mutations.txt" "${uuid}_mutations_pivot.csv" > "${uuid}_download_mutations.log";
            echo "uuid,project,data_category,data_type,download_dir,samples,mutation_table,pivot_table" > "${uuid}_mutations_metadata.csv";
            echo "${uuid},${project},${data_category},${data_type},${download_dir},${samples},${params.resultsDir}/${params.batchName}/${uuid}/data_download/mutations/${uuid}_mutations.txt,${params.resultsDir}/${params.batchName}/${uuid}/data_download/mutations/${uuid}_mutations_pivot.csv" >> "${uuid}_mutations_metadata.csv"
        """
    stub:
        """
        touch "${uuid}_mutations.txt"
        touch "${uuid}_mutations_pivot.csv"
        touch "${uuid}_mutations_metadata.csv"
        """

} 

// Concatenate all the metadata mutation table in one file
process mergeMutationsMetadata{

    //conda 'containers/conda_envs/merge_tables.yml' is managed by label
    label 'merge_tables'

    publishDir "${params.resultsDir}/${params.batchName}/", pattern: "downloaded_mutation_metadata.csv", mode: 'copy', overwrite: true
    
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

    label "r_download"

    publishDir "${params.resultsDir}/${params.batchName}/${uuid}/data_download/methylation/", pattern: "${uuid}_methylation*", mode: 'copy', overwrite: true
    
    input:
        tuple val(uuid),val(project),val(gdc_type),val(gdc_platform),val(download_dir),path(samples)
    output:
        tuple val(uuid),val(project),val(gdc_type),val(gdc_platform),val(download_dir),path(samples),file("${uuid}_methylation_manifest.txt"), file("${uuid}_methylations.txt"), file("${uuid}_methylation_metadata.csv")
    script:
    log.info "Downloading methulation for $uuid, $project"
        """
            if ! test -f ${samples}; then
            touch ${samples}
            fi
            Rscript '${baseDir}/bin/r/download_methylation_gdc.R' -p '${project}'  -t '${gdc_type}' --platform '${gdc_platform}' -d '${download_dir}' --manifest_outpath '${uuid}_methylation_manifest.txt' --pathlist_outpath '${uuid}_methylation_paths.txt' --header_outpath '${uuid}_methylation_header.txt' --sample_list ${samples}
            bash '${baseDir}/bin/bash/join_methylation_gdc.sh'  "${uuid}_methylations.txt" "${uuid}_methylation_paths.txt"
	        cat  '${uuid}_methylation_header.txt' "${uuid}_methylations.txt" > "${uuid}_methylations_labeled.txt"
            mv "${uuid}_methylations_labeled.txt" "${uuid}_methylations.txt";
            echo "uuid,project,gdc_type,gdc_platform,download_dir,samples,methylation_manifest,methylation_table" > "${uuid}_methylation_metadata.csv";
            echo "${uuid},${project},${gdc_type},${gdc_platform},${download_dir},${samples},${params.resultsDir}/${params.batchName}/${uuid}/data_download/methylation/${uuid}_methylation_manifest.txt,${params.resultsDir}/${params.batchName}/${uuid}/data_download/methylation/${uuid}_methylations.txt" >> "${uuid}_methylation_metadata.csv"
        
"""

    stub:
        """
        touch "${uuid}_methylation_manifest.txt"
        touch "${uuid}_methylations.txt"
        touch "${uuid}_methylation_metadata.csv"
        """


}

process mergeMethylationMetadata{
    // label uses conda environment
    label 'merge_tables'

    publishDir "${params.resultsDir}/${params.batchName}/", pattern: "downloaded_methylation_metadata.csv", mode: 'copy', overwrite: true
    
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

    label "r_download"

    publishDir "${params.resultsDir}/${params.batchName}/${uuid}/data_download/clinical", pattern: "*.csv", mode: 'copy', overwrite: true
    
    input:
        tuple val(uuid),val(project),val(data_category),val(data_type),val(data_format)
    output:
        tuple val(uuid),val(project),val(data_category),val(data_type),val(data_format),path("*")
    script:
    log.info "Downloading clinical data for $uuid, $project"
        """
            Rscript '${baseDir}/bin/r/download_clinical_tcga.R' ${project}  "${data_category}" "${data_type}" "${data_format}" "." 
        """

    stub:
        """
        touch clinical.csv
        """
}



process downloadCNV{

    label "r_download"
    publishDir "${params.resultsDir}/${params.batchName}/${uuid}/data_download/cnv/", pattern: "${uuid}.*", mode: 'copy', overwrite: true
    
    input:
        tuple val(uuid),val(project),val(workflow_type),path(samples)
    output:
        tuple val(uuid),val(project),val(workflow_type), path(samples) ,file("${uuid}.rds"),file("${uuid}.csv"),file("${uuid}_cnv_metatada.csv")
    script:
        """
        log.info "Downloading CNV for $uuid, $project"
            if ! test -f ${samples}; then
            touch ${samples}
            fi
            Rscript '${baseDir}/bin/r/download_cnv_tcga.R' -p ${project} --analysis_workflow_type ${workflow_type} --sample_list ${samples} --output_rds "${uuid}.rds" --output_table ${uuid}.csv > "${uuid}_cnv_downloads.log";
            echo "uuid,project,workflow_type,samples,output_rds,output_table" > "${uuid}_cnv_metatada.csv";
            echo "${uuid},${project},${workflow_type},${samples},${params.resultsDir}/${params.batchName}/${uuid}/data_download/cnv/${uuid}.rds,${params.resultsDir}/${params.batchName}/${uuid}/data_download/cnv/${uuid}.csv"  >> "${uuid}_cnv_metatada.csv"
        """
    stub:
        """
        touch "${uuid}.rds"
        touch "${uuid}.csv"
        touch "${uuid}_cnv_metatada.csv"
        """
}

process mergeCNVMetadata{
    // label uses conda environment
    label 'merge_tables'

    publishDir "${params.resultsDir}/${params.batchName}/", pattern: "downloaded_cnv_metadata.csv", mode: 'copy', overwrite: true
    
    input:
        val(results)
    output:
        tuple val(results), path("downloaded_cnv_metadata.csv")

    script:
	values = results.flatten().collect{"$it"}.join(',')
    """
    fit.py merge-tables -o downloaded_cnv_metadata.csv -t ${values}
    """

}



workflow downloadCNVWf{
    take: channelCNV
    main:
            dcnv = downloadCNV(channelCNV)
            mergeCNVMetadata(dcnv.map{it -> it[-1]}.collect())
    emit: dcnv
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
        //                                 }
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
        //                             }

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
        //        cond                             item.value.samples,
        //                                 )
        //                             }

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
        //                             }
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
        cnv_tcgabiolinks = Channel.fromPath(params.download_metadata).splitJson(path: "cnv_tcgabiolinks")


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
                                                file(item.value.samples),
                                            )
                                        }
        
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
                                                file(item.value.samples)
                                            )
                                        }
        
            downloadMutationsWf(channelMutation)
        }
        
        if (cnv_tcgabiolinks!=[null]){
            channelCNV = cnv_tcgabiolinks.map{
                                            item -> tuple(
                                                item.key,
                                                item.value.project,
                                                item.value.workflow_type,
                                                file(item.value.samples)
                                            )
                                        }
        
            downloadCNVWf(channelCNV)
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
                                                file(item.value.samples),
                                            )
                                        }
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
                                        }
            downloadClinicalWf(channelClinical)
    }
}

def get_expression_recount3 ( x ) {
    if (x.value.get('expression_recount3')){println('expression_recount3 exists')

    def filePath = params.inputFile ?: null  // Assuming you have a parameter named inputFile
    //println(x.value.get('expression_recount3').samples)
    //println(x.value.get('expression_recount3').samples==null)
    if (x.value.get('expression_recount3').samples != null) {
        fileObj = file(filePath)
    } else {
        // Handle null case appropriately
        println "File path is null."
        fileObj = 'NA'
    }


    a = tuple(x.key,
                                                        x.value.get('expression_recount3').project,
                                                        x.value.get('expression_recount3').project_home,
                                                        x.value.get('expression_recount3').organism,
                                                        x.value.get('expression_recount3').annotation,
                                                        x.value.get('expression_recount3').type,
                                                        fileObj)
    
    return a

    } else {return null}
}

def get_mutation_tcgabiolinks ( x ) {
    if (x.value.get('mutation_tcgabiolinks')){println('mutation_tcgabiolinks exists')
    a = tuple(x.key,
        x.value.get('mutation_tcgabiolinks').project,
        x.value.get('mutation_tcgabiolinks').data_category,
        x.value.get('mutation_tcgabiolinks').data_type,
        x.value.get('mutation_tcgabiolinks').download_dir,
        file(x.value.get('mutation_tcgabiolinks').samples))
    return a
    } else {return null}
}

def get_cnv_tcgabiolinks ( x ) {
    if (x.value.get('cnv_tcgabiolinks')){println('cnv_tcgabiolinks exists')
    a = tuple( x.key,
            x.value.get('cnv_tcgabiolinks').project,
            x.value.get('cnv_tcgabiolinks').workflow_type,
            file(x.value.get('cnv_tcgabiolinks').samples))
    return a
    } else {return null}
}

def get_methylation_gdc ( x ) {
    if (x.value.get('methylation_gdc')){println('methylation_gdc exists')

                a = tuple(  x.key,
                                x.value.get('methylation_gdc').project,
                                x.value.get('methylation_gdc').gdc_type,
                                x.value.get('methylation_gdc').gdc_platform,
                                x.value.get('methylation_gdc').download_dir,
                                file(x.value.get('methylation_gdc').samples))
                return a
                } else {return null}
            }

def get_clinical_tcgabiolinks ( x ) {
    if (x.value.get('clinical_tcgabiolinks')){println('clinical_tcgabiolinks exists')
        a = tuple( x.key,
                    x.value.get('clinical_tcgabiolinks').project,
                    x.value.get('clinical_tcgabiolinks').data_category,
                    x.value.get('clinical_tcgabiolinks').data_type,
                    x.value.get('clinical_tcgabiolinks').data_format )
        return a
    } else {return null}
    }

workflow fullDownloadWf{
    take:
        dataCh
    main:
        // Empty channels
        dr = Channel.empty()
        dmu = Channel.empty()
        dme = Channel.empty()
        dc = Channel.empty()

        // DOWNLOAD RECOUNT3
        dChRe = dataCh.map{
            it -> if(it.value.keySet().contains('expression_recount3')){
                tuple(
                    it.key,
                    it.value.get('expression_recount3').project,
                    it.value.get('expression_recount3').project_home,
                    it.value.get('expression_recount3').organism,
                    it.value.get('expression_recount3').annotation,
                    it.value.get('expression_recount3').type,
                    file(it.value.get('expression_recount3').samples)
                )
            }
        }
        dr = downloadRecount3Wf(dChRe)
        
        // DOWNLOAD MUTATIONS
        dChMu = dataCh.map{
            it -> if(it.value.keySet().contains('mutation_tcgabiolinks')){
                tuple(
                    it.key,
                    it.value.get('mutation_tcgabiolinks').project,
                    it.value.get('mutation_tcgabiolinks').data_category,
                    it.value.get('mutation_tcgabiolinks').data_type,
                    it.value.get('mutation_tcgabiolinks').download_dir,
                    file(it.value.get('mutation_tcgabiolinks').samples)
                )
            }
        }
        dmu = downloadMutationsWf(dChMu)
        
        //DOWNLOAD CNV   
        dChCNV = dataCh.map{
            it -> if(it.value.keySet().contains('cnv_tcgabiolinks')){
                tuple(
                    it.key,
                    it.value.get('cnv_tcgabiolinks').project,
                    it.value.get('cnv_tcgabiolinks').workflow_type,
                    file(it.value.get('cnv_tcgabiolinks').samples)
                )
            }
        }
        dcnv = downloadCNVWf(dChCNV)

        // DOWNLOAD METHYLATION
        dChMe = dataCh.map{
            it -> if(it.value.keySet().contains('methylation_gdc')){
                tuple(
                    it.key,
                    it.value.get('methylation_gdc').project,
                    it.value.get('methylation_gdc').gdc_type,
                    it.value.get('methylation_gdc').gdc_platform,
                    it.value.get('methylation_gdc').download_dir,
                    file(it.value.get('methylation_gdc').samples)
                )
            }
        }
        dme = downloadMethylationWf(dChMe)

        // DOWNLOAD CLINICAL
        dChCl = dataCh.map{
            it -> if(it.value.keySet().contains('clinical_tcgabiolinks')){
                tuple(
                    it.key,
                    it.value.get('clinical_tcgabiolinks').project,
                    it.value.get('clinical_tcgabiolinks').data_category,
                    it.value.get('clinical_tcgabiolinks').data_type,
                    it.value.get('clinical_tcgabiolinks').data_format
                )
            }
        }
        dc = downloadClinicalWf(dChCl)

    emit:
        dr
        dmu
        dme
        dc
        dcnv
    }


