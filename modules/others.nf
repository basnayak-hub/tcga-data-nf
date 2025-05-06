// This is a Nextflow script that processes metadata files for different types of data.

process preprocessExpressionMetadata {

    label "process_low"

    input:
    path templateFile

    output:
    file "processed_expression_metadata.csv"

    script:
    """
    echo ${templateFile}
    bash preprocess_metadata_test.sh ${templateFile} "processed_expression_metadata.csv" ${params.testDataFolder}
    """
}

process preprocessDragonMetadata {

    label "process_low"

    input:
    path templateFile

    output:
    file "processed_dragon_metadata.csv"

    script:
    """
    echo ${templateFile}
    bash preprocess_metadata_test.sh ${templateFile} "processed_dragon_metadata.csv" ${params.testDataFolder}
    """
}

process preprocessFullMetadata {

    label "process_low"

    input:
    path templateFile

    output:
    path "processed_full.config"

    script:
    """
    bash ${baseDir}/bin/preprocess_metadata_test.sh ${templateFile} processed_full.config ${params.testDataFolder}
    """
}
