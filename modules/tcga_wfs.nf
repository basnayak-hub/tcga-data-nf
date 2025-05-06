include { runGENIE3 ; runWGCNA } from './grns.nf'
/// The workflows
// I am adding a network generation workflow specific for our TCGA data. 

workflow PandaTCGAWf {
    take:
    data

    main:
    data
    panda = runTCGAPanda(data)

    emit:
    panda
}


workflow LionessPandaTCGAWf {
    take:
    data

    main:
    lio = runTCGALioness(data)

    emit:
    lio
}

workflow LionessOtterTCGAWf {
    take:
    data

    main:
    lio = runTCGAOtterLioness(data)
}

workflow DragonTCGAWf {
    take:
    data

    main:
    // If the second input is expression, the workflow will align the data, otherwise it will just assume the two files are aligned nad will run dragon on them.
    data
        .branch {
            toalign: it[3] == 'expression'
            aligned: it[3] != 'expression'
        }
        .set { dataAlignCh }

    // Align type1 and expression data
    toalignDragonCh = alignMethylationExpression(dataAlignCh.toalign)
    // Remove the original "expression" column
    toalignDragonCh.map { it -> tuple(it[0], it[1], it[2], it[3], it[5]) }

    // Concatenate the aligned data with the already aligned data
    dragonCh = dataAlignCh.aligned.concat(
        toalignDragonCh.map { it -> tuple(it[0], it[1], it[2], it[3], it[5]) }
    )
    // Run dragon on the aligned data (all of it)
    dragon = runTCGADragon(dragonCh)

    emit:
    dragon
}

workflow DragonLionessTCGAWf {
    take:
    data

    main:

    data
        .branch {
            toalign: it[3] == 'expression'
            aligned: it[3] != 'expression'
        }
        .set { dataAlignCh }

    toalignDragonCh = alignMethylationExpression(dataAlignCh.toalign)
    toalignDragonCh.map { it -> tuple(it[0], it[1], it[2], it[3], it[5]) }.view { "toalign: ${it}" }

    dragonCh = dataAlignCh.aligned
        .concat(
            toalignDragonCh.map { it -> tuple(it[0], it[1], it[2], it[3], it[5]) }
        )
        .view { "dragonCh: ${it}" }

    dragonLioness = runTCGALionessDragon(dragonCh)

    emit:
    dragonLioness
}

workflow AlpacaPandaTCGAWf {
    take:
    data

    main:
    alpacaCh = runTCGAAlpaca(data)

    emit:
    alpacaCh
}

workflow analyzeExpressionWf {
    take:
    data

    main:

    zooAnimals = Channel.from(params.zoo.animals)

    data
        .combine(zooAnimals)
        .branch {
            panda: it[-1] == 'panda'
            pandalioness: it[-1] == 'panda_lioness'
        }
        .set { zooAnalysisCh }



    pandaOut = PandaTCGAWf(zooAnalysisCh.panda.map { it -> tuple(it[0], it[1]) })

    // Check if ALPACA is in animals
    if (params.zoo.animals.contains('alpaca')) {
        // Generate all combinations (pairs) of elements from pandaOut using two separate channels
        channel1 = pandaOut
        channel2 = pandaOut

        // Cartesian product (every combination of two elements)
        pandaOutPairs = channel1
            .combine(channel2)
            .filter { it[0].compareTo(it[4]) < 0 }
        // Optional: Avoid self-pairs if needed

        // Run ALPACAWf with each pair of pandaOut
        alpacaOut = AlpacaPandaTCGAWf(pandaOutPairs)
    }



    lionessOut = LionessPandaTCGAWf(zooAnalysisCh.pandalioness.map { it -> tuple(it[0], it[1]) })

    //LionessOtterTCGAWf(zooAnalysisCh.otterlioness) 

    // Print genie3 parameter
    if (params.genie3.run_genie3) {
        runGENIE3(data.map { it -> tuple(it[0], it[1]) })
    }

    if (params.wgcna.run_wgcna) {
        runWGCNA(data.map { it -> tuple(it[0], it[1]) })
    }
}

workflow analyzeDragonWf {
    take:
    data

    main:
    zooAnimals = Channel.from(params.zoo.animals)

    data
        .combine(zooAnimals)
        .branch {
            dragon: it[-1] == 'dragon'
            dragonlioness: it[-1] == 'dragon_lioness'
        }
        .set { zooAnalysisCh }

    DragonTCGAWf(zooAnalysisCh.dragon.map { it -> tuple(it[0], it[1], it[2], it[3], it[4]) })
    DragonLionessTCGAWf(zooAnalysisCh.dragonlioness.map { it -> tuple(it[0], it[1], it[2], it[3], it[4]) })
}

process runTCGAPanda {

    label "netzoopy_panda", "process_panda"

    publishDir "${params.resultsDir}/${params.batchName}/${uuid}/analysis/panda/", pattern: 'panda*', mode: 'copy'

    input:
    tuple val(uuid), path(expression)

    output:
    tuple val(uuid), path(expression), path("panda_${uuid}.txt"), path("panda_${uuid}.log")

    script:
    log.info("... Running PANDA ${uuid},${expression}")
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

    label "netzoopy_pandalioness", "process_panda"

    publishDir "${params.resultsDir}/${params.batchName}/${uuid}/analysis/panda_lioness/", pattern: '{panda,lioness}*', mode: 'copy'

    input:
    tuple val(uuid), path(expression)

    output:
    tuple val(uuid), path(expression), path("panda.txt"), path("lioness/", type: 'dir')

    script:
    log.info("... Running PANDA LIONESS ${uuid},${expression}")
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

    label "netzoopy_pandalioness", "process_panda"

    publishDir "${params.resultsDir}/${params.batchName}/${uuid}/analysis/otter_lioness/", mode: 'copy'

    input:
    tuple val(uuid), path(expression)

    output:
    tuple val(uuid), path(expression), path("./otter/otter.h5"), path("./otter/lioness_otter/", type: 'dir')

    shell:
    log.info("... Running OTTER LIONESS ${uuid},${expression}")
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

    label 'r_base', 'process_medium'
    publishDir "${params.resultsDir}/${params.batchName}/${uuid}/analysis/dragon/", mode: 'copy', pattern: "${uuid}_${type1}_${type2}_dragon*", overwrite: true

    input:
    tuple val(uuid), val(type1), path(methylationData), val(type2), path(expressionData)

    output:
    tuple val(uuid), val(type1), path(methylationData), val(type2), path(expressionData), path("${uuid}_${type1}_${type2}_dragon_filtered_expression.csv")

    script:
    log.info("... Align ${uuid}: ${type2}:${expressionData} to ${type1}:${methylationData}")
    """
        Rscript ${baseDir}/bin/r/get_dragon_expression_data.r "${expressionData}" "${methylationData}" "${uuid}_${type1}_${type2}_dragon_filtered_expression.csv";
    """

    stub:
    """
        touch "${uuid}_${type1}_${type2}_dragon_filtered_expression.csv"
        """
}

process runTCGADragon {

    label "netzoopy_dragon", 'process_high'
    publishDir "${params.resultsDir}/${params.batchName}/${uuid}/analysis/dragon/", mode: 'copy', pattern: "${uuid}_${type1}_${type2}_dragon*", overwrite: true

    input:
    tuple val(uuid), val(type1), path(methylationData), val(type2), path(expressionAlignedData)

    output:
    tuple val(uuid), val(type1), path(methylationData), val(type2), path(expressionAlignedData), path("${uuid}_${type1}_${type2}_dragon_mat.tsv"), path("${uuid}_${type1}_${type2}_dragon_input.tsv"), path("${uuid}_${type1}_${type2}_dragon.log")

    script:
    log.info("... Running DRAGON ${uuid},${methylationData},${expressionAlignedData},${type1},${type2}")
    """
        run_dragon.py dragon -m ${methylationData} -e ${expressionAlignedData} -i "${uuid}_${type1}_${type2}_dragon_input.tsv" -o "${uuid}_${type1}_${type2}_dragon_mat.tsv" --type1 ${type1} --type2 ${type2} > "${uuid}_${type1}_${type2}_dragon.log"
        """

    stub:
    """
        touch "${uuid}_${type1}_${type2}_dragon_input.tsv"
        touch "${uuid}_${type1}_${type2}_dragon_mat.tsv" 
        touch "${uuid}_${type1}_${type2}_dragon.log"
        """
}



process runTCGALionessDragon {

    label "netzoopy_dragonlioness", 'process_high'
    publishDir "${params.resultsDir}/${params.batchName}/${uuid}/analysis/lioness_dragon_${type1}_${type2}/", mode: 'copy', pattern: "*", overwrite: true

    input:
    tuple val(uuid), val(type1), path(methylationData), val(type2), path(expressionAlignedData)

    output:
    tuple val(uuid), val(type1), path(methylationData), val(type2), path(expressionAlignedData), path("lioness_dragon_${type1}_${type2}/", type: 'dir'), path("${uuid}_${type1}_${type2}_lioness_dragon.log")

    script:
    log.info("... Running LIONESS DRAGON ${uuid},${expressionAlignedData}")
    """
    run_dragon.py lioness-dragon -m ${methylationData} -e ${expressionAlignedData} -o lioness_dragon_${type1}_${type2} --type1 ${type1} --type2 ${type2} > "${uuid}_${type1}_${type2}_lioness_dragon.log"
    """

    stub:
    """
        mkdir lioness_dragon_${type1}_${type2}
        touch "${uuid}_${type1}_${type2}_lioness_dragon.log"
        """
}

// Fix this after 0.9.3 netzoopy
process runTCGAPandaExplore {

    publishDir "${params.resultsDir}/${params.batchName}/${uuid}/analysis/figures/", pattern: '*.png', mode: 'copy'

    input:
    tuple val(uuid), path(expression), path(panda)

    output:
    tuple val(uuid), path(expression), path(panda), path('panda_scores.png')

    script:
    """
        explore.py plot-panda-scores ${panda} panda_scores.png --is_adj
        """
}


// Fix this after 0.9.3 netzoopy
process runTCGALionessExplore {

    publishDir "${params.resultsDir}/${params.batchName}/${uuid}/analysis/figures/", pattern: '*.png', mode: 'copy'

    input:
    tuple val(uuid), path(expression), path(panda), path(lioness)

    output:
    tuple val(uuid), path(expression), path(panda), path(lioness), path('lioness_scores.png'), path('correlation_panda_lioness.png'), path('single_correlation.png')

    script:
    """
        explore.py plot-lioness-scores ${lioness} lioness_scores.png --panda_filename ${panda} --correlation correlation_panda_lioness.png --singles single_correlation.png
        """
}


process runTCGAAlpaca {

    label "netzoor", 'process_ultra'


    publishDir "${params.resultsDir}/${params.batchName}/${uuid1}_${uuid2}/analysis/alpaca/", pattern: "*${uuid1}_${uuid2}*", mode: 'copy'

    input:
    tuple val(uuid1), path(expression1), path(panda1), path(log1), val(uuid2), path(expression2), path(panda2), path(log2)

    output:
    tuple val(uuid1), path(expression1), path(panda1), path(log1), val(uuid2), path(expression2), path(panda2), path(log2), path("membership_${uuid1}_${uuid2}.csv"), path("alpaca_${uuid1}_${uuid2}.Rds"), path("top_genes_${uuid1}_${uuid2}.txt"), path("alpaca_${uuid1}_${uuid2}.log")

    script:
    log.info("... Running ALPACA ${uuid1}-vs-${uuid2}")
    """
                Rscript ${baseDir}/bin/r/run_alpaca.r -b ${panda1} -p ${panda2} -o membership_${uuid1}_${uuid2}.csv -a alpaca_${uuid1}_${uuid2}.Rds -t top_genes_${uuid1}_${uuid2}.txt > alpaca_${uuid1}_${uuid2}.log
            """

    stub:
    """
            touch "membership_${uuid1}_${uuid2}.txt"
            touch "alpaca_${uuid1}_${uuid2}.Rds"
            touch "top_genes_${uuid1}_${uuid2}.txt"
            touch "alpaca_${uuid1}_${uuid2}.log"
            """
}
