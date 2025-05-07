

workflow prepareRecountWf {
    take:
    prepareRecountCh

    main:

    //prepareRecountCh.view{"hello: ${it}"}
    // Tissues channel
    channelTissues = Channel.from(params.tissues.entrySet())
        .map { item ->
            tuple(
                item.getKey(),
                item.getValue(),
            )
        }
        .transpose()
    // Batch correction channel

    channelBatchCorrection = Channel.from(params.batch_correction.entrySet())
        .map { item ->
            tuple(
                item.getKey(),
                item.getValue(),
            )
        }
        .transpose()

    // Combine all channels
    prepareCh = prepareRecountCh
        .combine(Channel.from(params.recount.norm))
        .combine(Channel.from(params.recount.min_tpm))
        .combine(Channel.from(params.recount.frac_samples))
        .combine(Channel.from(params.recount.th_purity))
        .combine(channelTissues, by: 0)
        .combine(channelBatchCorrection, by: 0)
        .combine(channelBatchCorrection, by: 0)

    // prepareTCGARecount
    readyRecountCh = prepareTCGARecount(prepareCh)

    emit:
    readyRecountCh
}


workflow prepareMethylationWf {
    take:
    prepareMethylationCh

    main:
    // Tissue channel
    channelTissues = Channel.from(params.tissues.entrySet())
        .map { item ->
            tuple(
                item.getKey(),
                item.getValue(),
            )
        }
        .transpose()


    promoterMethCh = GetGeneLevelPromoterMethylation(prepareMethylationCh, file(params.methylation.probe_map), file(params.methylation.tf_list))
    promoterMethCh
    readyMethCh = CleanMethylationData(promoterMethCh[0].combine(channelTissues, by: 0))

    emit:
    readyMethCh[0]
}

workflow prepareCNVWf {
    take:
    prepareCNVCh

    main:

    readyCNVCh = prepareCNV(prepareCNVCh, file(params.cnv.tf_list))

    emit:
    readyCNVCh
}

workflow prepareWf {
    println("Prepare Workflow, params.profile: ${params.profileName}")
    // Prepare Recount
    if (params.recount.metadata_prepare != '') {
        // Data channel

        if (params.profileName == 'testPrepare') {
            println('Test Prepare')
            prepareRecountCh = Channel.fromPath(params.recount.metadata_prepare)
                .splitCsv(header: true)
                .map { row ->
                    def expressionFilePath = "${workflow.projectDir}/${row.output_rds}"
                    tuple(row.uuid, row.project, file(expressionFilePath))
                }
                .view()
        }
        else {

            prepareRecountCh = Channel.fromPath(params.recount.metadata_prepare)
                .splitCsv(header: true)
                .map { row -> tuple(row.uuid, row.project, file(row.output_rds)) }
        }

        prepareRecountWf(prepareRecountCh)
    }
    // Prepare Methylation
    if (params.methylation.metadata_prepare != '') {
        // Data channel

        if (params.profileName == 'testPrepare') {
            prepareMethylationCh = Channel.fromPath(params.methylation.metadata_prepare)
                .splitCsv(header: true)
                .map { row ->
                    def methylationFilePath = "${workflow.projectDir}/${row.methylation_table}"
                    tuple(row.uuid, row.project, file(methylationFilePath))
                }
        }
        else {

            prepareMethylationCh = Channel.fromPath(params.methylation.metadata_prepare)
                .splitCsv(header: true)
                .map { row -> tuple(row.uuid, row.project, file(row.methylation_table)) }
        }

        prepareMethylationWf(prepareMethylationCh)
    }
}
process prepareTCGARecount {

    label 'prepare_expression', 'process_medium'

    publishDir "${params.resultsDir}/${params.batchName}/${uuid}/data_prepared/recount3/", pattern: "recount3_${uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.*", mode: 'copy', overwrite: true
    publishDir "${params.resultsDir}/${params.batchName}/${uuid}/data_prepared/recount3/", pattern: "recount3_${uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.log", mode: 'copy', overwrite: true
    publishDir "${params.resultsDir}/${params.batchName}/${uuid}/data_prepared/recount3/", pattern: "recount3_pca_${uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.png", mode: 'copy', overwrite: true

    input:
    tuple val(uuid), val(tcga_project), path(tcga_expression_fn), val(norm), val(min_tpm), val(frac_samples), val(th_purity), val(tissue_type), val(batch_correction), val(adjustment_variable)

    output:
    tuple val(uuid), val(tcga_project), path(tcga_expression_fn), val(norm), val(min_tpm), val(frac_samples), val(th_purity), val(tissue_type), val(batch_correction), val(adjustment_variable), path("recount3_${uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.rds"), path("recount3_${uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.txt"), path("recount3_pca_${uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.png"), path("recount3_${uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.log")

    script:
    log.info("... Preparing recount for ${uuid},${tcga_project},${tcga_expression_fn}")
    """
        Rscript '${baseDir}/bin/r/prepare_expression_recount.R' -p ${tcga_project}\
            -e ${tcga_expression_fn} \
                -r recount3_${uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.rds \
                -t recount3_${uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.txt \
                -f recount3_pca_${uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.png \
                --normalization ${norm} \
                --th_purity ${th_purity} \
                --min_tpm ${min_tpm} \
                --frac_samples ${frac_samples} \
                --batch_correction ${batch_correction} \
                --adjustment_variable ${adjustment_variable} \
                --tissue_type ${tissue_type} \
                ${params.recount.extra} >> "recount3_${uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.log"
"""
}


process GetGeneLevelPromoterMethylation {

    label 'prepare_methylation', 'process_medium'

    publishDir "${params.resultsDir}/${params.batchName}/${uuid}/data_prepared/methylation/", mode: 'copy', pattern: "${uuid}_tf_promoter_methylation_raw*", overwrite: true

    input:
    tuple val(uuid), val(project), path(methdata)
    path probe_map
    path tf_list

    output:
    tuple val(uuid), val(project), path(methdata), path("${uuid}_tf_promoter_methylation_raw.csv")
    path("${uuid}_tf_promoter_methylation_clean_raw_pca.png"), optional: true

    script:
    log.info("... Getting gene level promoter methylation ${uuid},${project}")
    def filter = tf_list.name != 'NO_FILE' ? "--tf_list ${tf_list}" : ""
    """
            if ${params.methylation.diagnostic_pca}; then
            save_opt="--diagnostic_pca ${uuid}_tf_promoter_methylation_raw_pca.png"
        else
            save_opt=""
        fi
        Rscript ${baseDir}/bin/r/get_gene_level_methylation.r -p ${project} -m ${methdata} -o "${uuid}_tf_promoter_methylation_raw.csv" --probemap ${probe_map} ${filter} \${save_opt} > "${uuid}_tf_promoter_methylation_raw.log" 
    """
}


process CleanMethylationData {

    label 'prepare_methylation', 'process_medium'

    publishDir "${params.resultsDir}/${params.batchName}/${uuid}/data_prepared/methylation/", mode: 'copy', pattern: "${uuid}_tf_promoter_methylation_clean_${tissueType}*", overwrite: true

    input:
    tuple val(uuid), val(project), path(methdata), path(rawMethylation), val(tissueType)

    output:
    tuple val(uuid), val(project), path(methdata), path(rawMethylation), val(tissueType), path("${uuid}_tf_promoter_methylation_clean_${tissueType}.csv"), path("${uuid}_tf_promoter_methylation_clean_${tissueType}.log")
    path("${uuid}_tf_promoter_methylation_clean_${tissueType}_pca.png"), optional: true
    

    script:
    log.info("... Cleaning methylation data ${uuid},${project}")
    """
        if ${params.methylation.diagnostic_pca}; then
            save_opt="--diagnostic_pca ${uuid}_tf_promoter_methylation_clean_${tissueType}_pca.png"
        else
            save_opt=""
        fi
        Rscript ${baseDir}/bin/r/clean_methylation_data.r -p ${project} -m  ${rawMethylation} --tissue_type "${tissueType}" -o "${uuid}_tf_promoter_methylation_clean_${tissueType}.csv" --th_missing ${params.methylation.th_missing} \${save_opt} > "${uuid}_tf_promoter_methylation_clean_${tissueType}.log"
    """
}

process prepareCNV {

    label 'prepare_cnv', 'process_medium'

    publishDir "${params.resultsDir}/${params.batchName}/${uuid}/data_prepared/cnv/", mode: 'copy', pattern: "${uuid}_cnv_clean.*", overwrite: true

    input:
    tuple val(uuid), val(project), path(cnvdata)
    path tf_list

    output:
    tuple val(uuid), val(project), path(cnvdata), path("${uuid}_cnv_clean.csv"), path("${uuid}_cnv_clean.log"), path("${uuid}_remove.txt")

    script:
    log.info("... Getting gene level promoter methylation ${uuid},${project}")
    def filter = tf_list.name != 'NO_FILE' ? "--tf_list ${tf_list}" : ""
    """
        Rscript ${baseDir}/bin/r/prepare_cnv.r -p ${project} -c ${cnvdata} -o "${uuid}_cnv_clean.csv"  ${filter} > "${uuid}_cnv_clean.log" 
        touch "${uuid}_remove.txt"
    """

    stub:
    """
        touch "${uuid}_cnv_clean.log" 
        """
}
