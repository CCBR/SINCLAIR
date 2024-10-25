process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_low'

    container "${params.containers.base}"

    input:
    path (samplesheet)
    path (contrast_samplesheet)
    val (run_cellranger)

    output:
    path '*gex_samplesheet.csv'                      , optional:true, emit:gex_samplesheet
    path '*atac_samplesheet.csv'                     , optional:true, emit:atac_samplesheet
    path '*contrast_samplesheet.csv'                 , emit:contrast_samplesheet
    path '*groups_samplesheet.csv'                   , emit:group_samplesheet

    script:
    """
    if [[ $run_cellranger == N ]]; then
        check_samplesheet_input_cellranger.py \\
        $samplesheet \\
        $contrast_samplesheet \\
        project
    else
        check_samplesheet.py \\
        $samplesheet \\
        $contrast_samplesheet \\
        project
    fi
    """

    stub:
    """
    cp $samplesheet STUB_gex_samplesheet.csv
    cp $contrast_samplesheet STUB_contrast_samplesheet.csv
    """
}
