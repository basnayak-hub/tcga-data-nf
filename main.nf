// enabling nextflow DSL v2
nextflow.enable.dsl=2

// import from modules
include { plotGGPlot } from './modules/rsteps.nf'


process PrintHelloOnFile {

    publishDir "${params.resultsDir}", pattern: "results.txt", mode: 'copy'

    input:
        path data

    output:
        path 'results.txt'

    """
        fit.py ${data} > results.txt
    """

}



workflow {
    if ("${params.pipeline}"=='R'){
        inFile = channel.from("${params.rInputFile}")
        plotGGPlot(inFile)}
    else
        channel.fromPath("${params.inputFile}") | PrintHelloOnFile
    
}