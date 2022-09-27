process prepareTCGARecount{

    publishDir "${params.resultsDir}/recount3/${tcga_uuid}/", pattern: "recount3_${tcga_uuid}_purity${th_purity}_mintpm${min_tpm}_fracsamples${frac_samples}_tissue${tissue_type}.*", mode: 'copy', overwrite: true
    publishDir "${params.resultsDir}/recount3/${tcga_uuid}/", pattern: "recount3_${tcga_uuid}_purity${th_purity}_mintpm${min_tpm}_fracsamples${frac_samples}_tissue${tissue_type}.log", mode: 'copy', overwrite: true
    
    input:
        tuple val(tcga_uuid),val(tcga_project),path(tcga_expression_fn),path(tcga_patient_fn),val(norm), val(min_tpm), val(frac_samples), val(th_purity), val(tissue_type)
    output:
        tuple val(tcga_uuid),val(tcga_project),path(tcga_expression_fn),path(tcga_patient_fn),val(norm), val(min_tpm), val(frac_samples), val(th_purity), val(tissue_type),\
                path("recount3_${tcga_uuid}_purity${th_purity}_mintpm${min_tpm}_fracsamples${frac_samples}_tissue${tissue_type}.rds"),\
                path("recount3_${tcga_uuid}_purity${th_purity}_mintpm${min_tpm}_fracsamples${frac_samples}_tissue${tissue_type}.txt"),\
                path("recount3_${tcga_uuid}_purity${th_purity}_mintpm${min_tpm}_fracsamples${frac_samples}_tissue${tissue_type}.log")

    script:
        """

        Rscript '${baseDir}/bin/R/prepare_expression_recount.R' -p ${tcga_project} -c ${tcga_patient_fn}\
            -e ${tcga_expression_fn} \
                -r recount3_${tcga_uuid}_purity${th_purity}_mintpm${min_tpm}_fracsamples${frac_samples}_tissue${tissue_type}.rds \
                -t recount3_${tcga_uuid}_purity${th_purity}_mintpm${min_tpm}_fracsamples${frac_samples}_tissue${tissue_type}.txt \
                --th_purity ${th_purity} \
                --min_tpm ${min_tpm} \
                --frac_samples ${frac_samples} \
                --tissue_type ${tissue_type} >> "recount3_${tcga_uuid}_purity${th_purity}_mintpm${min_tpm}_fracsamples${frac_samples}_tissue${tissue_type}.log"
                    """

}  
