// enabling nextflow DSL v2
nextflow.enable.dsl=2


// import from modules
include { downloadRecount; downloadMutations; downloadMethylation; downloadClinical} from './modules/download.nf'
include { prepareRecount} from './modules/prepare.nf'



// printing message of the day
motd = """
--------------------------------------------------------------------------
recount3-access-nf ($workflow.manifest.version)
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

        prepareRecountCh = Channel
                .fromPath( params.recount.metadata_prepare)
                .splitCsv( header: true)
                .map { row -> tuple( row.uuid, row.tcga_uuid, file(row.tcga_file),row.gtex_uuid, file(row.gtex_file) ) }.view()
        


        prepareRecount(prepareRecountCh.combine(params.recount.norm).combine(params.recount.th))
}

workflow {
    // We separate pipelines for downloading data and preparing it.
    // This allows for separate management of raw data and intermediate clean data
    if (params.pipeline == 'download')
        downloadWf()
    else if (params.pipeline == 'prepare')
        prepareWf()
    else
        downloadWf()
}