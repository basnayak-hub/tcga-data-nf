process prepareTCGARecount{

    publishDir "${params.resultsDir}/recount3/${tcga_uuid}/", pattern: "recount3_${tcga_uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.*", mode: 'copy', overwrite: true
    publishDir "${params.resultsDir}/recount3/${tcga_uuid}/", pattern: "recount3_${tcga_uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.log", mode: 'copy', overwrite: true
    publishDir "${params.resultsDir}/recount3/${tcga_uuid}/", pattern: "recount3_pca_${tcga_uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.png", mode: 'copy', overwrite: true
        
    input:
        tuple val(tcga_uuid),val(tcga_project),path(tcga_expression_fn),val(norm), val(min_tpm), val(frac_samples), val(th_purity), val(tissue_type), val(batch_correction), val(adjustment_variable)
    output:
        tuple val(tcga_uuid),val(tcga_project),path(tcga_expression_fn),val(norm), val(min_tpm), val(frac_samples), val(th_purity), val(tissue_type), val(batch_correction), val(adjustment_variable),\
                path("recount3_${tcga_uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.rds"),\
                path("recount3_${tcga_uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.txt"),\
                path("recount3_pca_${tcga_uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.png"),\
                path("recount3_${tcga_uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.log")

    script:
        """
        Rscript '${baseDir}/bin/r/prepare_expression_recount.R' -p ${tcga_project}\
            -e ${tcga_expression_fn} \
                -r recount3_${tcga_uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.rds \
                -t recount3_${tcga_uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.txt \
                -f recount3_pca_${tcga_uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.png \
                --normalization ${norm} \
                --th_purity ${th_purity} \
                --min_tpm ${min_tpm} \
                --frac_samples ${frac_samples} \
                --batch_correction ${batch_correction} \
                --adjustment_variable ${adjustment_variable} \
                --tissue_type ${tissue_type} >> "recount3_${tcga_uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.log"
"""

}


process prepareGTEXRecount{

    publishDir "${params.resultsDir}/recount3/${tcga_uuid}/", pattern: "recount3_${tcga_uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.*", mode: 'copy', overwrite: true
    publishDir "${params.resultsDir}/recount3/${tcga_uuid}/", pattern: "recount3_${tcga_uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.log", mode: 'copy', overwrite: true
    publishDir "${params.resultsDir}/recount3/${tcga_uuid}/", pattern: "recount3_pca_${tcga_uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.png", mode: 'copy', overwrite: true
        
    input:
        tuple val(gtex_uuid),val(gtex_project),path(tcga_expression_fn),path(tcga_patient_fn),val(norm), val(min_tpm), val(frac_samples), val(th_purity), val(tissue_type), val(batch_correction), val(adjustment_variable)
    output:
        tuple val(tcga_uuid),val(tcga_project),path(tcga_expression_fn),path(tcga_patient_fn),val(norm), val(min_tpm), val(frac_samples), val(th_purity), val(tissue_type), val(batch_correction), val(adjustment_variable),\
                path("recount3_${tcga_uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.rds"),\
                path("recount3_${tcga_uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.txt"),\
                path("recount3_pca_${tcga_uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.png"),\
                path("recount3_${tcga_uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.log")

    script:
        """
        Rscript '${baseDir}/bin/r/prepare_expression_recount.R' -p ${tcga_project} -c ${tcga_patient_fn}\
            -e ${tcga_expression_fn} \
                -r recount3_${tcga_uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.rds \
                -t recount3_${tcga_uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.txt \
                -f recount3_pca_${tcga_uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.png \
                --normalization ${norm} \
                --th_purity ${th_purity} \
                --min_tpm ${min_tpm} \
                --frac_samples ${frac_samples} \
                --batch_correction ${batch_correction} \
                --adjustment_variable ${adjustment_variable} \
                --tissue_type ${tissue_type} >> "recount3_${tcga_uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.log"
"""
}

process GetGeneLevelPromoterMethylation {

    publishDir "${params.resultsDir}/methylation/", mode: 'copy', pattern: "${uuid}_tf_promoter_methylation_raw.*",  overwrite: true

    // input: path to file of probes (rows) x samples (columns)
    // assume local for starters 
    input:
       tuple val(uuid), val(project), path(methdata)

    // output: file of samples (rows) x genes (columns)
    output:
        tuple val(uuid), val(project), path(methdata), path("${uuid}_tf_promoter_methylation_raw.csv")

    // return gene-level methylation measurements
    
    """
        Rscript ${baseDir}/bin/r/getGeneLevelMethylation.r ${project} ${methdata} "${uuid}_tf_promoter_methylation_raw.csv" ${params.methylation.probe_map} ${params.methylation.tf_list} > "${uuid}_tf_promoter_methylation_raw.log" 
    """

}


process CleanMethylationData {

    publishDir "${params.resultsDir}/methylation/", mode: 'copy', pattern: "${uuid}_tf_promoter_methylation_clean_${tissueType}.*",  overwrite: true

    input:
	tuple val(uuid), val(project), path(methdata), path(rawMethylation), val(tissueType)
 
    // output: file of samples (rows) x genes (columns)
    
    output:
        tuple val(uuid), val(project), path(methdata), path(rawMethylation), val(tissueType), path("${uuid}_tf_promoter_methylation_clean_${tissueType}.csv")
    
    """
         Rscript ${baseDir}/bin/r/cleanMethylationData.r ${baseDir} ${project} "${uuid}_tf_promoter_methylation_clean_${tissueType}.csv" ${rawMethylation} ${tissueType} ${params.methylation.to_npn} ${params.methylation.to_mval} > ${uuid}_tf_promoter_methylation_clean_${tissueType}.log
    """
}


workflow prepareRecountWf{
    take:
        prepareRecountCh
    main:

    // Tissues channel
    channelTissues = Channel.from(params.tissues.entrySet())
                                        .map{
                                            item -> tuple(
                                                item.getKey(),
                                                item.getValue()
                                            )
                                        }.transpose()
    // Batch correction channel
    channelBatchCorrection = Channel.from(params.batch_correction.entrySet())
                                        .map{
                                            item -> tuple(
                                                item.getKey(),
                                                item.getValue()
                                            )
                                        }.transpose().view()

    // Combine all channels
    prepareCh = (prepareRecountCh
                    .combine(Channel.from(params.recount.norm))
                    .combine(Channel.from(params.recount.min_tpm))
                    .combine(Channel.from(params.recount.frac_samples))
                   .combine(Channel.from(params.recount.th_purity)))
                   .combine(channelTissues, by: 0)
                   .combine(channelBatchCorrection, by: 0)
                   .combine(channelBatchCorrection, by: 0).view()

    // prepareTCGARecount
    readyRecountCh = prepareTCGARecount(prepareCh)

    emit:
        readyRecountCh
}



workflow prepareMethylationWf{
    take:
        prepareMethylationCh
    main:
    // Tissue channel
    channelTissues = Channel.from(params.tissues.entrySet())
                                        .map{
                                            item -> tuple(
                                                item.getKey(),
                                                item.getValue()
                                            )
                                        }.transpose().view()

    println('Empty')
    promoterMethCh =  GetGeneLevelPromoterMethylation(prepareMethylationCh)
    promoterMethCh.view()
    readyMethCh = CleanMethylationData(promoterMethCh.combine(channelTissues, by: 0))

    emit:
        readyMethCh

}

workflow prepareWf{
    main:

    if (params.recount.metadata_prepare!='') {

    // Data channel
    prepareRecountCh = Channel
                .fromPath( params.recount.metadata_prepare)
                .splitCsv( header: true)
                .map { row -> tuple( row.tcga_uuid,row.tcga_project, file(row.tcga_expression_file)) }.view()

    prepareRecountWf(prepareRecountCh)

    //prepareTCGARecount(prepareCh)
    }

    if (params.methylation.metadata_prepare!=''){    
    // Data channel
    println(params.methylation.metadata_prepare)
    
    println('Methylation Channel')
    prepareMethylationCh = Channel
                .fromPath( params.methylation.metadata_prepare)
                .splitCsv( header: true)
                .map { row -> tuple( row.tcga_uuid,row.tcga_project, file(row.tcga_methylation_file) ) }.view()


    // channelTissues = Channel.from(params.tissues.entrySet())
    //                                     .map{
    //                                         item -> tuple(
    //                                             item.getKey(),
    //                                             item.getValue()
    //                                         )
    //                                     }.transpose().view()
    // println('Empty')
    // promoterMethCh =  GetGeneLevelPromoterMethylation(prepareMethylationCh)
    // promoterMethCh.view()
    // readyMethCh = CleanMethylationData(promoterMethCh.combine(channelTissues, by: 0))
    prepareMethylationWf(prepareMethylationCh)


    }

}