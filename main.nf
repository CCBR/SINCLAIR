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

// Plugins
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

workflow.onComplete {
    if (!workflow.stubRun && !workflow.commandLine.contains('-preview')) {
        def message = Utils.spooker(workflow)
        if (message) {
            println message
        }
    }
}


workflow {
    LOG()
    validateParameters()
    PREPROCESS_EXQC ()
    GEX_EXQC (
        PREPROCESS_EXQC.out.ch_fqdir_h5,
        PREPROCESS_EXQC.out.group_samplesheet,
    )

}


workflow LOG {
    log.info """\
            SINCLAIR $workflow.manifest.version
            ===================================
            cmd line     : $workflow.commandLine
            start time   : $workflow.start
            NF outdir    : $params.outdir
            """
            .stripIndent()

    log.info paramsSummaryLog(workflow)
}
