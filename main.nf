// enabling nextflow DSL v2
nextflow.enable.dsl=2


// import from modules
//include { downloadWf; downloadRecount; downloadRecount3Wf; mergeRecountMetadata; downloadMutations; mergeMutationsMetadata; downloadMethylation; mergeMethylationMetadata; downloadClinical} from './modules/download.nf'
include {fullDownloadWf; downloadRecount;downloadRecount3Wf; downloadMutationsWf; downloadMethylationWf; downloadClinicalWf;  downloadWf} from './modules/download.nf'
include { prepareRecountWf; prepareWf; prepareTCGARecount; prepareMethylationWf; prepareCNVWf} from './modules/prepare.nf'
include { LionessPandaTCGAWf; LionessOtterTCGAWf; PandaTCGAWf; analyzeExpressionWf; analyzeDragonWf} from './modules/tcga_wfs.nf'
include { preprocessExpressionMetadata; preprocessDragonMetadata; preprocessFullMetadata} from './modules/others.nf'

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

// COPY CONFIG FILES
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




// SAVE CONFIG, copy the motd and the config files into the results directory
workflow saveConfig {
    mo = copyMotd()
    cf = Channel.from(workflow.configFiles)
    println "Config files: ${workflow.configFiles}"
    copyConfigFiles(cf,mo)
}


// FULL WORKFLOW
workflow fullWf{
    // DOWNLOAD
    // Data channel, metadata for Download

    if (params.profileName=='test'){
        println "Testing"
            metadata_full = preprocessFullMetadata("${params.full_metadata}")
        
    } else {
        metadata_full = Channel.fromPath(params.full_metadata)
    }
    dCh = metadata_full.splitJson()

    //dCh.map{it -> it.value.keySet()}.toString()
    //fullDownloadWf(dCh)
    
    // Full download
    fullDownCh = fullDownloadWf(dCh)

    fullDownCh.dr.map{it -> tuple(it[0], it[1], it[7], it[7])}.view{"full: ${it}"}
    fullDownCh.dmu.map{it -> tuple(it[0], it[1], it[7], it[7])}.view{"full: ${it}"}

    // view fullDownCh
    // PREPARE  (prepare expression and methylation)
    // Recount prepare (from fullDownCh.dr)
    fullDownCh.dr.map{it -> tuple(it[0], it[1], it[7], it[7])}.view()
    readyRecount = prepareRecountWf(fullDownCh.dr.map{it -> tuple(it[0], it[1], it[7])})
    
    // Methylation prepare (from fullDownCh.dme)
    fullDownCh.dme.map{it -> tuple(it[0], it[1], it[6], it[7])}
    readyMethylation = prepareMethylationWf(fullDownCh.dme.map{it -> tuple(it[0], it[1], it[7])})

    // Prepare CNV
    fullDownCh.dcnv.map{it -> tuple(it[0], it[1], it[7], it[4])}.view()
    readyCNV = prepareCNVWf(fullDownCh.dcnv.map{it -> tuple(it[0], it[1], it[4])})


    // ANALYZE
    
    readyRecount.map{it -> tuple(it[0],it[11])}
    readyMethylation.map{it -> tuple(it[0],it[5])}

    analyzeExpressionWf(readyRecount.map{it -> tuple(it[0],it[11])})

    methCh = readyMethylation.map{it -> tuple(it[0],'methylation',it[5])}.combine(readyRecount.map{it -> tuple(it[0],"expression",it[10])}, by:0)
    

    cnvCh = readyCNV.map{it -> tuple(it[0],"cnv",it[3])}.combine(readyRecount.map{it -> tuple(it[0],"expression", it[10])}, by:0)

    dragonCh = methCh.concat(cnvCh).view{"dragon from full: ${it}"}
    analyzeDragonWf(dragonCh) 
}


workflow {

    // First we copy the configuration files and motd into the resultsDir
    saveConfig()

    println "Pipeline: ${params.pipeline}"

    // We separate pipelines for downloading, preparing, analyzing  the data
    // This allows for separate management of raw data and intermediate clean data
    // Alongside, we have a full pipeline that does all the steps in one go
    if (params.pipeline == 'download'){
            // DOWNLOAD
            downloadWf()
        } else if (params.pipeline == 'prepare'){
            // PREPARE
            prepareWf() 
        } else if (params.pipeline == 'analyze'){
            // ANALYZE

            if (params.profileName=='testAnalyze'){
                // For the test we preprocess the metadata
                println "Test analyze"
                metadata_expression = preprocessExpressionMetadata("${params.metadata_expression}")
                metadata_dragon = preprocessDragonMetadata("${params.metadata_dragon}")
                
                } else {
                // For the full pipeline we read the metadata from the file
                metadata_expression = Channel.fromPath(params.metadata_expression)
                metadata_dragon = Channel.fromPath(params.metadata_dragon)
                }

            // PIPELINE FOR EXPRESSION ONLY (GRNs)
            data = metadata_expression
                        .splitCsv(header:true)
                        .map { row -> tuple(row.uuid, file("${row.expression}"))}

            analyzeExpressionWf(data)
            // PIPELINE FOR DRAGON: EXPRESSION and METHYLATION
            
            //dataDragon = metadata_dragon
            //            .splitCsv(header:true)
            //            .map { row -> tuple(row.uuid, file("${row.methylation}"), file//("${row.expression}"))}

            dataDragon = metadata_dragon
                        .splitCsv(header:true)
                        .map { row -> tuple(row.uuid, row.type1, file("${row.methylation}"),row.type2, file("${row.expression}"))}

            analyzeDragonWf(dataDragon)
        } else if (params.pipeline == 'full')
            // FULL PIPELINE
            fullWf()
    else
        // 
        error "Error: this pipeline name doesn't exist. \nChoose one between download/prepare/analyze/full as pipeline parameter"
        
}