
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
