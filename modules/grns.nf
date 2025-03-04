process runGENIE3 {

    label "r_genie3"

    publishDir "${params.resultsDir}/${params.batchName}/${uuid}/analysis/genie3/",  pattern: 'genie3*', mode: 'copy'

        input:
            tuple val(uuid), path(expression)

        output:
            tuple val(uuid), path(expression), path("genie3_${uuid}.txt"), path("genie3_${uuid}.log")
        
        script:
        log.info "... Running GENIE3 $uuid,$expression"
            """
                ${baseDir}/bin/bash/remove_dot_ensembl.sh ${expression} ${uuid}.nodot.txt;
                Rscript ${baseDir}/bin/r/run_genie3.r -i ${uuid}.nodot.txt -o genie3_${uuid}.txt --tf_list ${params.genie3.tf_list} --n_cores ${params.genie3.n_cores} --tree_method ${params.genie3.tree_method} --n_trees ${params.genie3.n_trees} --k ${params.genie3.k} > genie3_${uuid}.log
            """

        stub:
            """
            touch "${uuid}.nodot.txt"
            touch genie3_${uuid}.txt
            """
}