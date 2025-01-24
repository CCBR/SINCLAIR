process CELLRANGER_COUNT {
    tag "${id}"

    container 'litd/docker-cellranger:v7.2.0'

    input:
    tuple val(id), val(inDir)
    val(genome_dir)
    val(run_cellranger)

    output:
    tuple val(id), path('*/outs/filtered_feature_bc_matrix.h5'), emit:h5

    script:
    def args = task.ext.args ?: ''
    def localmem = task.memory.toGiga()
    """
    if [[ "$run_cellranger" == "N" ]]; then
        mkdir -p $id/outs/
        cp $inDir/filtered_feature_bc_matrix.h5 $id/outs/filtered_feature_bc_matrix.h5
    else
        cellranger count \
        --id $id \
        --fastqs $inDir \
        --transcriptome=${genome_dir} \
        --localcores=$task.cpus \
        --localmem=$localmem
    fi
    """

    stub:
    """
    mkdir -p $id/outs
    touch $id/outs/filtered_feature_bc_matrix.h5
    """
}
