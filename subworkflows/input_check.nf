/*
===================================================================
Assign Local Modules
===================================================================
*/
include {SAMPLESHEET_CHECK      } from '../modules/local/samplesheet_check'

/*
=======================================================================================================
RUN MAIN WORKFLOW
=======================================================================================================
*/
workflow INPUT_CHECK_GEX {
    take:
        samplesheet
        contrast_samplesheet
        run_cellranger

    main:
        SAMPLESHEET_CHECK (
            samplesheet,
            contrast_samplesheet,
            run_cellranger
        )

    emit:
        gex_samplesheet = SAMPLESHEET_CHECK.out.gex_samplesheet
        contrast_samplesheet = SAMPLESHEET_CHECK.out.contrast_samplesheet
        group_samplesheet = SAMPLESHEET_CHECK.out.group_samplesheet
}

workflow INPUT_CHECK_ATAC {
    take:
        samplesheet

    main:
        SAMPLESHEET_CHECK (samplesheet)
            .atac_samplesheet
            .splitCsv()
            .map { row ->
                    def id = row[0]
                    def inDir = row[1]
                    return [ id, inDir ]
                }
            .set { atac_input }
    emit:
        atac_input               // channel: [ val(id), path(inDir) ]
        atac_samplesheet = SAMPLESHEET_CHECK.out.atac_samplesheet
}


// Channel.fromPath( file(params.fastq_raw_sheet) )
//         .splitCsv()
//         .map { row ->
//             def id = row[0]
//             def read1 = file(row[1])
//             def read2 = file(row[2])

//             return [ id, read1, read2 ]
//         }
//         .groupTuple()
//         .into { sample_fastq_r1r2; sample_fastq_r1r2_2 }

// def create_mapped_channel(row){
//     // [WB_Lysis_1:atac-'0',  '/data/CCBR_Pipeliner/Pipelines/TechDev_scRNASeq_Dev2023/test_dir/WB_Lysis_Granulocytes_3p_Introns_8kCells_fastqs/WB_Lysis_1'|gex-'0',  '/data/CCBR_Pipeliner/Pipelines/TechDev_scRNASeq_Dev2023/test_dir/WB_Lysis_Granulocytes_3p_Introns_8kCells_fastqs/WB_Lysis_1']
//     // [WB_Lysis_2:gex-'0',  '/data/CCBR_Pipeliner/Pipelines/TechDev_scRNASeq_Dev2023/test_dir/WB_Lysis_Granulocytes_3p_Introns_8kCells_fastqs/WB_Lysis_2']
//     array=[]
//     sample=row.split(":")[0]
//     print(sample)

//     //     array = [meta]
//     //     if dt=="gex":
//     //         print("Y")
//     //         array = [meta.id,[dt]]
//     //     else:
//     //         print("N")
//     //     // } else {
//     //     //     if(!file(row.fastq_2).exists()){
//     //     //         exit1, "ERROR: missing FASTQ file R2\n${row.fastq_2}"
//     //     //     }
//     //     //     array=[meta,[file(row.fastq_1),file(row.fastq_2)]]
//     //     // }
//     // return array
// }
