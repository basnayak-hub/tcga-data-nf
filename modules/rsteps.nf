process plotGGPlot {
    publishDir "${params.resultsDir}/figures/", pattern: "scatter.pdf", mode: 'copy', overwrite: true

    input: 
        val inFile
    output:
        tuple val(inFile), path("scatter.pdf")
    """
        Rscript '${baseDir}/bin/plot.r' ${inFile} scatter.pdf
    """

}