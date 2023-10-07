
process runTCGAPanda {

//conda "/Users/violafanfani/Documents/uni-harvard/workflows/tcga-data-nf/containers/env.netzoopy.yml"

publishDir "${params.resultsDir}/${params.batchName}/${uuid}/",  pattern: 'panda*', mode: 'copy'

    input:
        tuple val(uuid), path(expression)

    output:
        tuple val(uuid), path(expression), path("panda_${uuid}.txt"), path("panda_${uuid}.log")
    
    shell:
        '''
            cat !{expression} | awk 'BEGIN { OFS=FS="\\t" } {  sub(/\\..*$/, "", \$1); print  }' > !{uuid}.nodot.txt;
            which netzoopy;
            netzoopy panda -e !{uuid}.nodot.txt -m !{params.zoo.motif} -p !{params.zoo.ppi} -o panda_!{uuid}.txt !{params.zoo.panda} > panda_!{uuid}.log
        '''

    stub:
        """
        touch "${uuid}.nodot.txt"
        touch panda_${uuid}.txt
        """
}

process runTCGALioness {

//conda "/Users/violafanfani/Documents/uni-harvard/workflows/tcga-data-nf/containers/env.netzoopy.yml"


publishDir "${params.resultsDir}/${params.batchName}/${uuid}/",  pattern: '{panda,lioness}*', mode: 'copy'

    input:
        tuple val(uuid), path(expression)

    output:
        tuple val(uuid), path(expression), path("panda.txt"), path("lioness/", type:'dir')
    
    shell:
        '''
            cat !{expression} | awk 'BEGIN { OFS=FS="\\t" } {  sub(/\\..*$/, "", \$1); print  }' > !{uuid}.nodot.txt;
            netzoopy lioness -e !{uuid}.nodot.txt -m !{params.zoo.motif} -p !{params.zoo.ppi} -op panda.txt -ol lioness/ !{params.zoo.panda_lioness} > panda_lioness_!{uuid}.log
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

//conda "/Users/violafanfani/Documents/uni-harvard/workflows/tcga-data-nf/containers/env.netzoopy.yml"

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


process runTCGADragon {

    publishDir "${params.resultsDir}/${params.batchName}/${uuid}/", mode: 'copy', pattern:"${uuid}_dragon*",  overwrite: true

    input:
        tuple val(uuid),path(methylationData),path(expressionData)
    	// output: file of samples (rows) x genes (columns)
    output:
        tuple val(uuid),path(methylationData),path(expressionData),path("${uuid}_dragon_filtered_expression.csv"), path("${uuid}_dragon_mat.tsv"), path("${uuid}_dragon_input.tsv")
    
    """
    Rscript ${baseDir}/bin/r/get_dragon_expression_data.r "${expressionData}" "${methylationData}" "${uuid}_dragon_filtered_expression.csv";
    python -u ${baseDir}/bin/run_dragon.py ${methylationData} "${uuid}_dragon_filtered_expression.csv" "." ${uuid} > dragon_log.txt
    """
}


/// The workflows

// I am adding a network generation workflow specific for our TCGA data. 

workflow PandaTCGAWf {

    take:data
    main:
    lio = runTCGAPanda(data)
    //lio.take(3).view()
    //runTCGAPandaExplore(lio.take(3))

}


workflow LionessPandaTCGAWf {

    take:data
    main:
    lio = runTCGALioness(data)
    //lio.take(3).view()
    //runTCGAPandaExplore(lio.take(3))
    //runTCGALionessExplore(lio)

}



workflow LionessOtterTCGAWf {

    take:data
    main:
    lio = runTCGAOtterLioness(data)

}

workflow DragonTCGAWf {

    take:data
    main:
    dragon = runTCGADragon(data)

}




process RunLionessDragon{

    publishDir "${params.resultsDir}", mode: 'move',  overwrite: true

    input:
        tuple val(uuid),path(methylationData),path(expressionData)

    output:
        tuple val(uuid),path(uuid), path("results_ld.txt")

    """
    python -u ${baseDir}/bin/py/run_lioness_dragon.py ${methylationData} ${expressionData} ${uuid}/lioness ${uuid} > results_ld.txt
    """
}


workflow analyzeExpressionWf{
    take:
        data

    main:

    // defaults results directory
    //batchName = params.batchName ? params.batchName : "batch-${params.workflow}-null"

    // 
    //if (!params.metadata) exit 1, "requires a CSV metadata file."
    // Data channel
    // format uuid, file(network)

    

    zooAnimals = Channel.from(params.zoo.animals)

    data.combine(zooAnimals).branch {
                    panda: it[-1] == 'panda'
                    pandalioness: it[-1] == 'panda_lioness'
                    //otter: it[-1] == 'otter'
                    otterlioness: it[-1] == 'otter_lioness'    
                }.set { zooAnalysisCh }

    PandaTCGAWf(zooAnalysisCh.panda)

    LionessPandaTCGAWf(zooAnalysisCh.pandalioness)

    LionessOtterTCGAWf(zooAnalysisCh.otterlioness) 

    //DragonTCGAWf()

}

workflow analyzeDragonWf{
    take:
        data

    main:
    zooAnimals = Channel.from(params.zoo.animals)

    data.combine(zooAnimals).branch {
                    dragon: it[-1] == 'dragon'
                    dragonlioness: it[-1] == 'dragon_lioness'
                }.set { zooAnalysisCh }

    DragonTCGAWf(zooAnalysisCh.dragon)
    //DragonLionessTCGAWf(zooAnalysisCh.dragonlioness)


}