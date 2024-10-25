/*
=======================================================================================================
Validate inputs
=======================================================================================================
*/

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
//input on command line
if (params.input)    { ch_input    = file(params.input)    } else { exit 1, 'Input samplesheet not specified!' }
if (params.contrast) { ch_contrast = file(params.contrast) } else { exit 1, 'Contrast samplesheet not specified!' }

/*
=======================================================================================================
Assign local subworkflows
=======================================================================================================
*/

/*
=======================================================================================================
Assign Local Modules
=======================================================================================================
*/
include { SEURAT_PREPROCESS                             } from '../modules/local/seurat_preprocess.nf'
include { SEURAT_MERGE                                  } from '../modules/local/seurat_merge.nf'
include { BATCH_CORRECT_HARMONY                         } from '../modules/local/batch_correction_harmony.nf'
include { BATCH_CORRECT_RPCA                            } from '../modules/local/batch_correction_rpca.nf'
include { BATCH_CORRECT_CCA                             } from '../modules/local/batch_correction_cca.nf'
include { BATCH_CORRECT_SCVI                            } from '../modules/local/batch_correction_scvi.nf'
include { BATCH_CORRECT_LIGER                           } from '../modules/local/batch_correction_liger.nf'
include { BATCH_CORRECT_INTEGRATION                     } from '../modules/local/batch_correction_integration.nf'
/*
=======================================================================================================
RUN MAIN WORKFLOW
=======================================================================================================
*/
workflow GEX_EXQC {
    take:
        ch_fqdir_h5
        group_samplesheet
    main:
        // Set output path to relative, species
        outdir_path = Channel.fromPath(params.outdir,relative:true)
        // Run Seurat for individual samples
        SEURAT_PREPROCESS (
            ch_fqdir_h5,
            params.species,
            params.qc_filtering,
            params.nCount_RNA_max,
            params.nCount_RNA_min,
            params.nFeature_RNA_max,
            params.nFeature_RNA_min,
            params.percent_mt_max,
            params.percent_mt_min,
            params.run_doublet_finder,
            params.npcs,
            params.script_preprocess,
            params.script_functions
        )

        // creates metadata
        // creates tuple of group_id_key: all RDS files from INPUT_CHECK_GEX
        // which matches those samples:
        //     example: INPUT_CHECK_GEX.out.group_samplesheet; SEURAT_SINGLE.out.rds
        //     example: sample1:group1_group2; sample1:sample1.rds
        //     example: sample1:group1_group2_group3; sample1:sample1.rds
        //     example: sample2:group1_group2; sample2:sample2.rds
        // output: group1_group2: [sample1.rds,sample2.rds]
        // output: group1_group2_group3: [sample1.rds]
        ch_groups = group_samplesheet
            .splitCsv( header:true, sep:',', strip:true )
            .map { row ->
                    def key = row["keyid"]
                    def sample = row["sampleid"]
                    return [sample, key]
                }
            .combine(SEURAT_PREPROCESS.out.rds, by: 0)
            .map { sample, key, rds_file -> tuple( key, rds_file ) }
            .groupTuple()

        SEURAT_MERGE (
            ch_groups,
            ch_input,
            params.species,
            params.npcs,
            params.vars_to_regress,
            params.script_merge,
            params.script_functions
        )

        // Run batch corrections
        BATCH_CORRECT_HARMONY (
            SEURAT_MERGE.out.rds,
            params.species,
            params.npcs,
            params.vars_to_regress,
            params.resolution_list,
            params.script_bc_harmony,
            params.script_functions
        )

        // Run batch corrections
        BATCH_CORRECT_RPCA (
            SEURAT_MERGE.out.rds,
            params.species,
            params.npcs,
            params.vars_to_regress,
            params.resolution_list,
            params.script_bc_rpca,
            params.script_functions
        )

        // Run batch corrections
        BATCH_CORRECT_CCA (
            SEURAT_MERGE.out.rds,
            params.species,
            params.npcs,
            params.vars_to_regress,
            params.resolution_list,
            params.script_bc_cca,
            params.script_functions
        )

        // Run batch corrections
        BATCH_CORRECT_SCVI (
            SEURAT_MERGE.out.rds,
            params.species,
            params.npcs,
            params.vars_to_regress,
            params.resolution_list,
            params.conda_path,
            params.python_path,
            params.script_scvi,
            params.script_functions
        )

        // Run batch corrections
        BATCH_CORRECT_LIGER (
            SEURAT_MERGE.out.rds,
            params.species,
            params.npcs,
            params.vars_to_regress,
            params.resolution_list,
            params.script_liger,
            params.script_functions
        )

        // Integrate batch corrections
        BATCH_CORRECT_INTEGRATION (
            SEURAT_MERGE.out.rds,
            BATCH_CORRECT_HARMONY.out.rds,
            BATCH_CORRECT_RPCA.out.rds,
            BATCH_CORRECT_CCA.out.rds,
            BATCH_CORRECT_SCVI.out.rds,
            BATCH_CORRECT_LIGER.out.rds,
            params.species,
            params.npcs,
            params.resolution_list,
            params.script_bc_integration,
            params.script_functions
        )

    emit:
        harmony_rds         = BATCH_CORRECT_HARMONY.out.rds
        rpca_rds            = BATCH_CORRECT_RPCA.out.rds
        cca_rds             = BATCH_CORRECT_CCA.out.rds
}
