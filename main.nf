// enabling nextflow DSL v2
nextflow.enable.dsl=2


// import from modules
include { downloadRecount; downloadMutations; downloadMethylation; downloadClinical} from './modules/download.nf'
include { prepareTCGARecount} from './modules/prepare.nf'



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


workflow downloadWf{
    main:
        modalities = params.download_metadata.keySet()

         if (modalities.contains('expression_recount3')){
            // Channel for recount3 data download
            channelRecount = Channel.from(params.download_metadata.expression_recount3.entrySet())
                                        .map{
                                            item -> tuple(
                                                item.getKey(),
                                                item.getValue().project,
                                                item.getValue().project_home,
                                                item.getValue().organism,
                                                item.getValue().annotation,
                                                item.getValue().type
                                            )
                                        }.view()
            
            downloadRecount(channelRecount)
        }

        if (modalities.contains('mutation_tcgabiolinks')){
        // Add here channels for mutations and methylation
        channelMutation = Channel.from(params.download_metadata.mutation_tcgabiolinks.entrySet())
                                    .map{
                                        item -> tuple(
                                            item.getKey(),
                                            item.getValue().project,
                                            item.getValue().data_category,
                                            item.getValue().data_type,
                                            item.getValue().download_dir,
                                        )
                                    }.view()

            downloadMutations(channelMutation)
        } 
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

            downloadMethylation(channelMethylation)
        }

        if (modalities.contains('clinical_tcgabiolinks')){
        channelClinical = Channel.from(params.download_metadata.clinical_tcgabiolinks.entrySet())
                                    .map{
                                        item -> tuple(
                                            item.getKey(),
                                            item.getValue().project,
                                            item.getValue().data_category,
                                            item.getValue().data_type,
                                            item.getValue().data_format,
                                        )
                                    }.view()

            downloadClinical(channelClinical)
        }


}


workflow prepareWf{
    main:
    // Tissue channels
    channelTissues = Channel.from(params.tissues.entrySet())
                                        .map{
                                            item -> tuple(
                                                item.getKey(),
                                                item.getValue()
                                            )
                                        }.transpose()
    // Batch correction channel
    channelBatchCorrection = Channel.from(params.batch_correction.entrySet())
                                        .map{
                                            item -> tuple(
                                                item.getKey(),
                                                item.getValue()
                                            )
                                        }.transpose()
    
    // Data channel
    prepareRecountCh = Channel
                .fromPath( params.recount.metadata_prepare)
                .splitCsv( header: true)
                .map { row -> tuple( row.tcga_uuid,row.tcga_project, file(row.tcga_expression_file),file(row.tcga_patient_file) ) }

    // Parameter channel
    prepareCh = (prepareRecountCh
                    .combine(Channel.from(params.recount.norm))
                    .combine(Channel.from(params.recount.min_tpm))
                    .combine(Channel.from(params.recount.frac_samples))
                   .combine(Channel.from(params.recount.th_purity)))
                   .combine(channelTissues, by: 0)
                   .combine(channelBatchCorrection, by: 0)
                   .combine(channelBatchCorrection, by: 0).view()

    prepareTCGARecount(prepareCh)
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
        downloadWf()
    else if (params.pipeline == 'prepare')
        prepareWf()
    else
        downloadWf()
        
}