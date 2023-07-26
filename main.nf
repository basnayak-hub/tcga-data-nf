// enabling nextflow DSL v2
nextflow.enable.dsl=2


// import from modules
include { downloadWf; downloadRecount; mergeRecountMetadata; downloadMutations; mergeMutationsMetadata; downloadMethylation; mergeMethylationMetadata; downloadClinical} from './modules/download.nf'
include { prepareWf; prepareTCGARecount} from './modules/prepare.nf'
include { LionessPandaTCGAWf; LionessOtterTCGAWf; PandaTCGAWf  } from './modules/tcga_wfs.nf'


// printing message of the day
motd = """
--------------------------------------------------------------------------
tcga-data-nf ($workflow.manifest.version)
--------------------------------------------------------------------------
Session ID   : $workflow.sessionId
Pipeline: $params.pipeline

--------------------------------------------------------------------------
Environment information
--------------------------------------------------------------------------
Container    : $workflow.container
Config files : $workflow.configFiles
Project dir  : $workflow.projectDir
Work dir     : $workflow.workDir
Launch dir   : $workflow.launchDir
Command line : $workflow.commandLine
Repository   : $workflow.repository
CommitID     : $workflow.commitId
Revision     : $workflow.revision
--------------------------------------------------------------------------
"""

log.info motd



workflow devWf{
    main:
        modalities = params.download_metadata.keySet()

        if (modalities.contains('methylation_gdc')){
        channelMethylation = Channel.from(params.download_metadata.methylation_gdc.entrySet())
                                    .map{
                                        item -> tuple(
                                            item.getKey(),
                                            item.getValue().project,
                                            item.getValue().gdc_type,
                                            item.getValue().gdc_platform,
                                            item.getValue().download_dir,
                                        )
                                    }.view()

            dme = downloadMethylation(channelMethylation)
            mergeMethylationMetadata(dme.map{it -> it[-1]}.collect())
        
}}

// worflow analyzeDragonWf {
//     RunLionessDragon(methylation, expression)
// }

// workflow methylationWf{

//     // Channel for cancer names
//     channelPancancer = Channel.from(params.pancancer_ids.entrySet())
//         .map{
//              item -> tuple(
//                 item.getKey(),
//              item.getValue().project,
// 	     item.getValue().tissueType)
//              }.view()
	     

//     x = DownloadMethylation(channelPancancer) | FilterTissueType | GetGeneLevelPromoterMethylation | CleanMethylationData
//     y = DownloadExpression(channelPancancer) | GetGeneExpression
//     z = x.join(y)
//     RunLionessDragon(z) | UploadLionessDragons | CleanupLionessDragons
// }

workflow analyzeWf{
    main:

    // defaults results directory
    //batchName = params.batchName ? params.batchName : "batch-${params.workflow}-null"

    // 
    //if (!params.metadata) exit 1, "requires a CSV metadata file."
    // Data channel
    // format uuid, file(network)
    data = Channel
                .fromPath(params.metadata, checkIfExists: true)
                .splitCsv(header:true)
                .map { row -> tuple(row.uuid, file("${row.expression}"))}

    dataMethylation = Channel
                .fromPath(params.metadata, checkIfExists: true)
                .splitCsv(header:true)
                .map { row -> tuple(row.uuid, file("${row.methylation}"))}

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


process copyConfigFiles{
    publishDir "${params.resultsDir}",pattern:"${params.logInfoFile}", mode: 'copy', overwrite: true
    input:
        file(config_file)
        file(out_conf)
    output:
        path(out_conf)
    script:
        """
        cat $config_file >> $out_conf
        """
}

process copyMotd{
    output:
        path(params.logInfoFile)
    shell:
        """
        echo '${motd}' > ${params.logInfoFile}
        """
}

workflow saveConfig {
    mo = copyMotd()
    cf = Channel.from(workflow.configFiles)
    copyConfigFiles(cf,mo)
}

workflow {

    // First we copy the configuration files and motd into the resultsDir
    //copyConfigFiles(cf.map($it -> tuple(file($it), $it.getName())))
    //copyConfigFiles(Channel.of(file(workflow.configFiles))).view()
    saveConfig()

    // We separate pipelines for downloading data and preparing it.
    // This allows for separate management of raw data and intermediate clean data
    if (params.pipeline == 'download')
        downloadWf(

        )
    else if (params.pipeline == 'analyze')
        analyzeWf()
    else if (params.pipeline == 'prepare')
        prepareWf()
    else if (params.pipeline == 'dev')
        devWf()
    else
        downloadWf()
        
}