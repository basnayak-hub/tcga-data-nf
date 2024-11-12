
process CollectTableStats {
    tag "Final stats "

    publishDir "${resultsDir}/data", pattern: "stats.csv", mode: 'copy'

    input:
        val(uuid)
        val(results)

    output:
        path "stats.csv"

    script:
    id = uuid.collect{"$it"}.join(',')
	values = results.flatten().collect{"$it"}.join(',')
    """
        # collect stats
        collect_stats.py collect-generic stats.csv --data ${values} --uuid ${id}
    """
}


process preprocessExpressionMetadata {

    input:
        path(templateFile)

    output:
        file("processed_expression_metadata.csv")

    script:
    """
    echo ${templateFile}
    bash preprocess_metadata_test.sh ${templateFile} "processed_expression_metadata.csv" ${params.testDataFolder}
    """
}

process preprocessDragonMetadata {

    input:
        path(templateFile)

    output:
        file("processed_dragon_metadata.csv")

    script:
    """
    echo ${templateFile}
    bash preprocess_metadata_test.sh ${templateFile} "processed_dragon_metadata.csv" ${params.testDataFolder}
    """
}

process preprocessFullMetadata {
    input:
    path(templateFile)

    output:
    path "processed_full.config"

    script:
    """
    bash ${baseDir}/bin/preprocess_metadata_test.sh ${templateFile} processed_full.config ${params.testDataFolder}
    """
}