include { runGENIE3; runWGCNA } from './grns.nf'

process runTCGAPanda {

    label "netzoopy_panda"
    //conda "/Users/violafanfani/Documents/uni-harvard/workflows/tcga-data-nf/containers/env.netzoopy.yml"

    publishDir "${params.resultsDir}/${params.batchName}/${uuid}/analysis/panda/",  pattern: 'panda*', mode: 'copy'

        input:
            tuple val(uuid), path(expression)

        output:
            tuple val(uuid), path(expression), path("panda_${uuid}.txt"), path("panda_${uuid}.log")
        
        script:
        log.info "... Running PANDA $uuid,$expression"
            """
                ${baseDir}/bin/bash/remove_dot_ensembl.sh ${expression} ${uuid}.nodot.txt;
                netzoopy panda -e ${uuid}.nodot.txt -m "${params.zoo.motif}" -p ${params.zoo.ppi} -o panda_${uuid}.txt ${params.zoo.panda} > panda_${uuid}.log
            """

        stub:
            """
            touch "${uuid}.nodot.txt"
            touch panda_${uuid}.txt
            """
}

process runTCGALioness {

    //conda "/Users/violafanfani/Documents/uni-harvard/workflows/tcga-data-nf/containers/env.netzoopy.yml"
    label "netzoopy_pandalioness"

    publishDir "${params.resultsDir}/${params.batchName}/${uuid}/analysis/panda_lioness/",  pattern: '{panda,lioness}*', mode: 'copy'

        input:
            tuple val(uuid), path(expression)

        output:
            tuple val(uuid), path(expression), path("panda.txt"), path("lioness/", type:'dir')
        
        script:
        log.info "... Running PANDA LIONESS $uuid,$expression"
            """
                ${baseDir}/bin/bash/remove_dot_ensembl.sh ${expression} ${uuid}.nodot.txt;
                netzoopy lioness -e ${uuid}.nodot.txt -m ${params.zoo.motif} -p ${params.zoo.ppi} -op panda.txt -ol lioness/ ${params.zoo.panda_lioness} > panda_lioness_${uuid}.log
            """
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

    publishDir "${params.resultsDir}/${params.batchName}/${uuid}/analysis/otter_lioness/", mode: 'copy'

        input:
            tuple val(uuid), path(expression)

        output:
            tuple val(uuid), path(expression), path("./otter/otter.h5"), path("./otter/lioness_otter/", type:'dir')
        
        shell:
            log.info "... Running OTTER LIONESS $uuid,$expression"
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

process alignMethylationExpression {

    label 'r_base'
    publishDir "${params.resultsDir}/${params.batchName}/${uuid}/analysis/dragon/", mode: 'copy', pattern:"${uuid}_dragon*",  overwrite: true

    input:
        tuple val(uuid),path(methylationData),path(expressionData)
    	// output: file of samples (rows) x genes (columns)
    output:
        tuple val(uuid),path(methylationData),path(expressionData),path("${uuid}_dragon_filtered_expression.csv")
    
    script:
        log.info "... Align methylation and expression $uuid,$expressionData, $methylationData"
    """
        Rscript ${baseDir}/bin/r/get_dragon_expression_data.r "${expressionData}" "${methylationData}" "${uuid}_dragon_filtered_expression.csv";
    """
    stub:
        """
        touch "${uuid}_dragon_filtered_expression.csv"
        """
}

process runTCGADragon {

    label "netzoopy_dragon"
    publishDir "${params.resultsDir}/${params.batchName}/${uuid}/analysis/dragon/", mode: 'copy', pattern:"${uuid}_dragon*",  overwrite: true

    input:
        tuple val(uuid),path(methylationData),path(expressionData), path(expressionAlignedData)
    	// output: file of samples (rows) x genes (columns)
    output:
        tuple val(uuid),path(methylationData),path(expressionData),path(expressionAlignedData), path("${uuid}_dragon_mat.tsv"), path("${uuid}_dragon_input.tsv"), path("${uuid}_dragon.log")
    
    script:
        log.info "... Running DRAGON $uuid,$expressionAlignedData"
        """
        run_dragon.py dragon -m ${methylationData} -e ${expressionAlignedData} -i "${uuid}_dragon_input.tsv" -o "${uuid}_dragon_mat.tsv" > "${uuid}_dragon.log"
        """
    stub:
        """
        touch "${uuid}_dragon_input.tsv"
        touch "${uuid}_dragon_mat.tsv" 
        touch "${uuid}_dragon.log"
        """
}



process runTCGALionessDragon{

    label "netzoopy_dragonlioness"
    publishDir "${params.resultsDir}/${params.batchName}/${uuid}/analysis/lioness_dragon/", mode: 'copy', pattern:"*",  overwrite: true

    input:
        tuple val(uuid),path(methylationData),path(expressionData), path(expressionAlignedData)

    output:
        tuple val(uuid),path(methylationData),path(expressionData), path(expressionAlignedData), path("lioness_dragon/", type:'dir'), path("${uuid}_lioness_dragon.log")

    script:
    log.info "... Running LIONESS DRAGON $uuid,$expressionAlignedData"
    """
    run_dragon.py lioness-dragon -m ${methylationData} -e ${expressionAlignedData} -o lioness_dragon > "${uuid}_lioness_dragon.log"
    """

    stub:
        """
        mkdir lioness_dragon
        touch "${uuid}_lioness_dragon.log"
        """
}

process runTCGAPandaExplore {

publishDir "${params.resultsDir}/${params.batchName}/${uuid}/analysis/figures/",  pattern: '*.png', mode: 'copy'

    input:
        tuple val(uuid), path(expression), path(panda)

    output:
        tuple val(uuid), path(expression), path(panda),path('panda_scores.png')
    """
    explore.py plot-panda-scores ${panda} panda_scores.png --is_adj
    """
}


// Fix this after 0.9.3 netzoopy
process runTCGALionessExplore {

publishDir "${params.resultsDir}/${params.batchName}/${uuid}/analysis/figures/",  pattern: '*.png', mode: 'copy'

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

workflow PandaTCGAWf {

    take:data
    main:
    data
    lio = runTCGAPanda(data)
    //lio.take(3).view()
    //runTCGAPandaExplore(lio.take(3))

}


workflow LionessPandaTCGAWf {

    take:data
    main:
    //motif = channel.fromPath("${params.zoo.motif}")
    //ppi = channel.fromPath("${params.zoo.ppi}")
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
    dragonCh = alignMethylationExpression(data)
    dragon = runTCGADragon(dragonCh)

}

workflow DragonLionessTCGAWf {

    take:data
    main:
    dragonCh = alignMethylationExpression(data)
    dragon = runTCGALionessDragon(dragonCh)

}

workflow analyzeExpressionWf{
    take:
        data

    main:

    zooAnimals = Channel.from(params.zoo.animals)

    data.combine(zooAnimals).branch { 
                    panda: it[-1] == 'panda'
                    pandalioness: it[-1] == 'panda_lioness'
                    //otter: it[-1] == 'otter'
                    //otterlioness: it[-1] == 'otter_lioness'    
                }.set { zooAnalysisCh }
                
    

    PandaTCGAWf(zooAnalysisCh.panda.map{it -> tuple(it[0], it[1])})

    LionessPandaTCGAWf(zooAnalysisCh.pandalioness.map{it -> tuple(it[0], it[1])})

    //LionessOtterTCGAWf(zooAnalysisCh.otterlioness) 
    // Print genie3 parameter
    
    if (params.genie3.run_genie3){
        runGENIE3(data.map{it -> tuple(it[0], it[1])})
    }

    if (params.wgcna.run_wgcna){
        runWGCNA(data.map{it -> tuple(it[0], it[1])})
    }

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

    DragonTCGAWf(zooAnalysisCh.dragon.map{it -> tuple(it[0], it[1], it[2])})
    DragonLionessTCGAWf(zooAnalysisCh.dragonlioness.map{it -> tuple(it[0], it[1], it[2])})

}