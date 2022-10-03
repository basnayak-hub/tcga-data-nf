process downloadRecount{

    publishDir "${params.resultsDir}/${uuid}/recount3/", pattern: "${uuid}.rds", mode: 'copy', overwrite: true
    
    input:
        tuple val(uuid),val(project),val(project_home),val(organism),val(annotation),val(type)
    output:
        tuple val(uuid),val(project),val(project_home),val(organism),val(annotation),val(type),file("${uuid}.rds")
    script:
        """
            Rscript '${baseDir}/bin/r/download_expression_recount.R' ${project} ${project_home} ${organism} ${annotation} ${type} "${uuid}.rds"
        """

}

process downloadMutations{

    publishDir "${params.resultsDir}/${uuid}/mutations/", pattern: "${uuid}_mutations*", mode: 'copy', overwrite: true
    
    input:
        tuple val(uuid),val(project),val(data_category),val(data_type),val(download_dir)
    output:
        tuple val(uuid),val(project),val(data_category),val(data_type),val(download_dir),file("${uuid}_mutations.txt"),file("${uuid}_mutations_pivot.csv")
    script:
        """
            Rscript '${baseDir}/bin/r/download_mutation_tcga.R' ${project}  "${data_category}" "${data_type}" "${download_dir}" "${uuid}_mutations.txt" "${uuid}_mutations_pivot.csv"
        """

} 

process downloadMethylation{
    publishDir "${params.resultsDir}/${uuid}/methylation/", pattern: "${uuid}_methylation*", mode: 'copy', overwrite: true
    
    input:
        tuple val(uuid),val(project),val(gdc_type),val(gdc_platform),val(download_dir)
    output:
        tuple val(uuid),val(project),val(gdc_type),val(gdc_platform),val(download_dir),file("TCGA_methylation_paths.txt"),file("${uuid}_methylations.txt")
    script:
        """
            Rscript '${baseDir}/bin/r/download_methylation_gdc.R' ${project}  "${gdc_type}" "${gdc_platform}" "${download_dir}" "${params.resultsDir}methylation"
            bash '${baseDir}/bin/bash/join_methylation_gdc.sh'  "${uuid}_methylations.txt"
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
}
