/*
=======================================================================================================
Validate inputs
=======================================================================================================
*/

/*
=======================================================================================================
Assign local subworkflows
=======================================================================================================
*/
include { INPUT_CHECK_ATAC                               } from '../subworkflows/input_check'

/*
=======================================================================================================
Assign Local Modules
=======================================================================================================
*/
include { SAMPLESHEET_CHECK                              } from '../modules/local/samplesheet_check.nf'

/*
=======================================================================================================
RUN MAIN WORKFLOW
=======================================================================================================
*/
workflow ATAC_EXQC {
    main:
        ch_input = file(params.input, checkIfExists: true)

        // Set output path to relative
        outdir_path = Channel.fromPath(params.outdir,relative:true)

        // Read in samplesheet
        INPUT_CHECK_ATAC (
            ch_input
        )

        // Run cellranger count

    emit:
        samplesheet        = INPUT_CHECK_ATAC.out.atac_samplesheet
}
