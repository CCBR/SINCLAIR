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
include { INPUT_CHECK_GEX                               } from '../subworkflows/input_check'

/*
=======================================================================================================
Assign Local Modules
=======================================================================================================
*/
include { SAMPLESHEET_CHECK                             } from '../modules/local/samplesheet_check.nf'
include { CELLRANGER_COUNT                              } from '../modules/local/cellranger_count_gex.nf'
/*
=======================================================================================================
RUN MAIN WORKFLOW
=======================================================================================================
*/
workflow PREPROCESS_EXQC {
    main:
        ch_input = file(params.input, checkIfExists: true)
        ch_contrast = file(params.contrast, checkIfExists: true)

        // Read in samplesheet, contrast manifest
        INPUT_CHECK_GEX (
            ch_input,
            ch_contrast,
            params.run_cellranger
        )

        // create metadata
        // creates GEX mapped input of sample:gex_input_dir
        ch_meta = INPUT_CHECK_GEX.out.gex_samplesheet
            .splitCsv( header:true, sep:',', strip:true )
            .map { row ->
                    def id = row["sample"]
                    def inDir = row["gex_input_dir"]
                    return [ id, inDir ]
                }

        // Run cellranger count
        if (params.run_cellranger) {

            CELLRANGER_COUNT (
                ch_meta,
                params.genome_dir,
                params.save_cellranger_extra_files
            )
            // [id, fastq input dir, h5 file]
            ch_fqdir_h5 = ch_meta.join(CELLRANGER_COUNT.out.h5)
        } else {
            ch_fqdir_h5 = ch_meta.map { id, inDir ->
                [id, '', file("${inDir}/*.h5", checkIfExists: true)]
            }
        }
    emit:
        group_samplesheet   = INPUT_CHECK_GEX.out.group_samplesheet
        ch_fqdir_h5 = ch_fqdir_h5
}
