
process runTCGALioness {

publishDir "${params.resultsDir}/${params.batchName}/${uuid}/",  pattern: '{panda,lioness}*', mode: 'copy'

    input:
        tuple val(uuid), path(expression)

    output:
        tuple val(uuid), path(expression), path("panda.txt"), path("lioness/", type:'dir')
    
    shell:
        '''
            cat !{expression} | awk 'BEGIN { OFS=FS="\\t" } {  sub(/\\..*$/, "", \$1); print  }' >> !{uuid}.nodot.txt;
            netzoopy lioness -e !{uuid}.nodot.txt -m !{params.zoo.motif} -p !{params.zoo.ppi} -op panda.txt -ol lioness/ !{params.zoo.panda_lioness} >> panda_lioness_!{uuid}.log
        '''
    stub:
        """
        touch "${uuid}.nodot.txt"
        touch "panda_lioness_${uuid}.log"
        mkdir lioness
        touch panda.txt
        """
}


process runTCGAOtterLioness {

publishDir "${params.resultsDir}/${params.batchName}/${uuid}/", mode: 'copy'

    input:
        tuple val(uuid), path(expression)

    output:
        tuple val(uuid), path(expression), path("./otter/otter.h5"), path("./otter/lioness_otter/", type:'dir')
    
    shell:
        '''
            cat !{expression} | awk 'BEGIN { OFS=FS="\\t" } {  sub(/\\..*$/, "", \$1); print  }' >> !{uuid}.nodot.txt;
            netzoopy otterlioness -e !{uuid}.nodot.txt -m !{params.zoo.motif} -p !{params.zoo.ppi} -of otter/ !{params.zoo.otter_lioness} >> otter_lioness_!{uuid}.log
        '''

    stub:
        """
        mkdir otter
        mkdir otter/lioness_otter
        touch "${uuid}.nodot.txt"
        touch otter_lioness_${uuid}.log
        """
}




// add this back once there is panda output again
//output:
//        tuple val(uuid), path(expression), path(motif), path(ppi), val(start), val(end),path("panda.txt"), path("lioness/", type:'dir')
    
// 
process runTCGAPandaExplore {

publishDir "${params.resultsDir}/${params.batchName}/${uuid}/figures/",  pattern: '*.png', mode: 'copy'

    input:
        tuple val(uuid), path(expression), path(panda)

    output:
        tuple val(uuid), path(expression), path(panda),path('panda_scores.png')
    """
    explore.py plot-panda-scores ${panda} panda_scores.png --is_adj
    """
}


// add this back
// input:
//     tuple val(uuid), path(expression), path(motif), path(ppi), val(start), val(end), path(lioness)


// Fix this after 0.9.3 netzoopy
process runTCGALionessExplore {

publishDir "${params.resultsDir}/${params.batchName}/${uuid}/figures/",  pattern: '*.png', mode: 'copy'

    input:
        tuple val(uuid), path(expression), path(panda), path(lioness)

    output:
        tuple val(uuid), path(expression), path(panda),path(lioness),path('lioness_scores.png'),path('correlation_panda_lioness.png'),path('single_correlation.png')
    """
    explore.py plot-lioness-scores ${lioness} lioness_scores.png --panda_filename ${panda} --correlation correlation_panda_lioness.png --singles single_correlation.png
    """
}

/// The workflows

// I am adding a network generation workflow specific for our TCGA data. 

workflow LionessPandaTCGAWf {

    take:data
    main:
    lio = runTCGALioness(data)
    lio.take(3).view()
    runTCGAPandaExplore(lio.take(3))
    runTCGALionessExplore(lio)

}



workflow LionessOtterTCGAWf {

    take:data
    main:
    lio = runTCGAOtterLioness(data)

}