process prepareRecount{

    publishDir "${params.resultsDir}/recount3/${uuid}/", pattern: "*.csv", mode: 'copy', overwrite: true
    publishDir "${params.resultsDir}/recount3/${uuid}/", pattern: "${uuid}_${uuid_tcga}_${uuid_gtex}_recount_prepare.log", mode: 'copy', overwrite: true
    
    input:
        tuple val(uuid),val(uuid_tcga),path(tcga_fn),val(uuid_gtex),path(gtex_fn),val(norm),val(th)
    output:
        tuple val(uuid),val(uuid_tcga),path(tcga_fn),val(uuid_gtex),path(gtex_fn),val(norm),val(th),path("${uuid_tcga}_norm${norm}_th${th}.csv"),path("${uuid_gtex}_norm${norm}_th${th}.csv"),path("${uuid}_${uuid_tcga}_${uuid_gtex}_recount_prepare.log")
    script:
        """
            Rscript '${baseDir}/bin/R/prepare_gtex_tcga_recount.R' ${tcga_fn} "${uuid_tcga}_norm${norm}_th${th}.csv" ${gtex_fn} "${uuid_gtex}_norm${norm}_th${th}.csv" ${norm} ${th} >> "${uuid}_${uuid_tcga}_${uuid_gtex}_recount_prepare.log"
        """

}

