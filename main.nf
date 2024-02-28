// enabling nextflow DSL v2
nextflow.enable.dsl=2


// import from modules
//include { downloadWf; downloadRecount; downloadRecount3Wf; mergeRecountMetadata; downloadMutations; mergeMutationsMetadata; downloadMethylation; mergeMethylationMetadata; downloadClinical} from './modules/download.nf'
include {fullDownloadWf; downloadRecount;downloadRecount3Wf; downloadMutationsWf; downloadMethylationWf; downloadClinicalWf;  downloadWf} from './modules/download.nf'
include { prepareRecountWf; prepareWf; prepareTCGARecount; prepareMethylationWf} from './modules/prepare.nf'
include { LionessPandaTCGAWf; LionessOtterTCGAWf; PandaTCGAWf; analyzeExpressionWf; analyzeDragonWf} from './modules/tcga_wfs.nf'


// printing message of the day
motd = """
--------------------------------------------------------------------------
tcga-data-nf ($workflow.manifest.version)
--------------------------------------------------------------------------
Session ID   : $workflow.sessionId
Pipeline: $params.pipeline
Results dir  : $params.resultsDir
Batch name   : $params.batchName

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

            dme = downloadMethylationWf(channelMethylation)
            mergeMethylationMetadata(dme.map{it -> it[-1]}.collect())
        }   

         if (modalities.contains('expression_recount3')){
            // Channel for recount3 data download
            downloadRecount3Wf(params.download_metadata.expression_recount3)
    }
}

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



workflow fullWf{
    // DOWNLOAD
    // Data channel, metadata for Download
    dCh = Channel.fromPath(params.full_metadata).splitJson().view()

    //dCh.map{it -> it.value.keySet()}.toString()
    fullDownloadWf(dCh)
    // // Full download
    //fullDownCh = fullDownloadWf(dCh)

    // // PREPARE  (prepare expression and methylation)
    // // Recount prepare (from fullDownCh.dr)
    // fullDownCh.dr.map{it -> tuple(it[0], it[1], it[7], it[7])}
    // readyRecount = prepareRecountWf(fullDownCh.dr.map{it -> tuple(it[0], it[1], it[7])})
    // // Methylation prepare (from fullDownCh.dme)
    // fullDownCh.dme.map{it -> tuple(it[0], it[1], it[6], it[7])}
    // readyMethylation = prepareMethylationWf(fullDownCh.dme.map{it -> tuple(it[0], it[1], it[7])})

    // //fullDownloadWf(params.full_metadata.entrySet().each{it -> it.getKey()},dCh.map{it -> it.getValue()})
    // //Channel.from(params.full_metadata.keySet()).view() uuids = params.full_metadata.keySet()

    // // ANALYZE
    // readyRecount.map{it -> tuple(it[0],it[11])}
    // readyMethylation.map{it -> tuple(it[0],it[5])}

    // //analyzeWf(readyRecount.map{it -> tuple(it[0],it[11])}, readyMethylation.map{it -> tuple(it[0],it[5])})

}




workflow {

    // First we copy the configuration files and motd into the resultsDir
    //copyConfigFiles(cf.map($it -> tuple(file($it), $it.getName())))
    //copyConfigFiles(Channel.of(file(workflow.configFiles))).view()
    //remove comment!!
    //saveConfig()

    // We separate pipelines for downloading data and preparing it.
    // This allows for separate management of raw data and intermediate clean data
    if (params.pipeline == 'download'){

    downloadWf()
    } else if (params.pipeline == 'analyze'){
        data = Channel
                    .fromPath(params.metadata, checkIfExists: true)
                    .splitCsv(header:true)
                    .map { row -> tuple(row.uuid, file("${row.expression}"))}

        analyzeExpressionWf(data)
        dataDragon = Channel
                    .fromPath(params.metadata_dragon, checkIfExists: true)
                    .splitCsv(header:true)
                    .map { row -> tuple(row.uuid, file("${row.methylation}"), file("${row.expression}"))}

        analyzeDragonWf(dataDragon)
    } else if (params.pipeline == 'prepare')
        prepareWf()
    else if (params.pipeline == 'dev')
        downloadWf()
    else if (params.pipeline == 'full')
        fullWf()
    else
        downloadWf()
        
}