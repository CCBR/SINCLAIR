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
//

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
workflow GEX {
    main:
        PREPROCESS_EXQC ()
        GEX_EXQC (
            PREPROCESS_EXQC.out.ch_meta,
            PREPROCESS_EXQC.out.group_samplesheet,
            PREPROCESS_EXQC.out.h5
        )

}

workflow ATAC {
    main:
        PREPROCESS_EXQC ()
        ATAC_EXQC ()
    emit:
        samplesheet         = PREPROCESS_EXQC.out.samplesheet
}

// //
// // WORKFLOW: Run file sample type 1
// //
// workflow GEX {
//     main:
//         GEX_EXQC()
//     emit:
//         file1                   = GEX_EXQC.out.file1
//         file2                   = GEX_EXQC.out.file2
// }


// //
// // WORKFLOW: Run file sample type 2
// //
// workflow ATAC {
//     main:
//         ATAC_EXQC()
//     emit:
//         file1                   = SUBWK_2_EXQC.out.file1
//         file2                   = SUBWK_2_EXQC.out.file2
// }
