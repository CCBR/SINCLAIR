#!/usr/bin/env nextflow
/*
===================================================================
    CCBR/SINCLAIR
===================================================================
    GitHub: https://github.com/CCBR/SINCLAIR
-------------------------------------------------------------------

*/

// Using DSL-2
nextflow.enable.dsl=2

log.info """\
        S I N C L A I R   P I P E L I N E
        ===================================
        NF version   : $nextflow.version
        runName      : $workflow.runName
        username     : $workflow.userName
        configs      : $workflow.configFiles
        profile      : $workflow.profile
        cmd line     : $workflow.commandLine
        start time   : $workflow.start
        projectDir   : $workflow.projectDir
        launchDir    : $workflow.launchDir
        workDir      : $workflow.workDir
        homeDir      : $workflow.homeDir
        outDir       : $params.outdir
        """
        .stripIndent()

/*
===================================================================
    Import workflows
===================================================================
*/
include { PREPROCESS_EXQC                           } from './workflows/pre_process'
include { GEX_EXQC                                  } from './workflows/gex'
include { ATAC_EXQC                                 } from './workflows/atac'
// include { VDJ_EXQC                            } from '.workflows/sCRNA_vdj'
// include { CITE_EXQC                           } from '.workflows/sCRNA_cite'

/*
===================================================================
    Set handlers
===================================================================
*/
workflow.onComplete {
    if (!workflow.stubRun && !workflow.commandLine.contains('-preview')) {
        println "Running spooker"
        def message = Utils.spooker(workflow)
        if (message) {
            println message
        }
    }
}
workflow.onError {
    if (!workflow.stubRun && !workflow.commandLine.contains('-preview')) {
        println "Running spooker (failed)"
        def message = Utils.spooker(workflow)
        if (message) {
            println message
        }
    }
}

/*
===================================================================
    Define workflows
===================================================================
*/
//
// WORKFLOW: Run initialization on input samples, check manifests files
//
workflow {
    main:
        PREPROCESS_EXQC ()
        GEX_EXQC (
            PREPROCESS_EXQC.out.ch_fqdir_h5,
            PREPROCESS_EXQC.out.group_samplesheet,
        )

}
