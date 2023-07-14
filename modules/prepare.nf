process prepareTCGARecount{

    publishDir "${params.resultsDir}/recount3/${tcga_uuid}/", pattern: "recount3_${tcga_uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.*", mode: 'copy', overwrite: true
    publishDir "${params.resultsDir}/recount3/${tcga_uuid}/", pattern: "recount3_${tcga_uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.log", mode: 'copy', overwrite: true
    publishDir "${params.resultsDir}/recount3/${tcga_uuid}/", pattern: "recount3_pca_${tcga_uuid}_purity0${th_purity.toString().substring(2)}_norm${norm}_mintpm${min_tpm}_fracsamples0${frac_samples.toString().substring(2)}_tissue${tissue_type}_batch${batch_correction.minus('.').minus(' ').minus('-').minus('_')}_adj${adjustment_variable.minus('.').minus(' ').minus('-').minus('_')}.png", mode: 'copy', overwrite: true
        
    input:
        tuple val(tcga_uuid),val(tcga_project),path(tcga_expression_fn),path(tcga_patient_fn),val(norm), val(min_tpm), val(frac_samples), val(th_purity), val(tissue_type), val(batch_correction), val(adjustment_variable)
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
