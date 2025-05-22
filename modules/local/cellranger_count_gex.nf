process CELLRANGER_COUNT {
    tag "${id}"

    input:
    tuple val(id), val(inDir)
    val(genome_dir)
    val(run_cellranger)

    output:
    tuple val(id), path('*/outs/filtered_feature_bc_matrix.h5'), emit: h5
    //tuple val(id), path('*/outs/*.bam*'),                        emit: bam_bai

    script:
    def args = task.ext.args ?: ''
    def localmem = task.memory.toGiga()
    """
    cellranger count \
        --id $id \
        --fastqs $inDir \
        --transcriptome=${genome_dir} \
        --localcores=$task.cpus \
        --localmem=$localmem
    """

    stub:
    """
    mkdir -p $id/outs
    touch $id/outs/filtered_feature_bc_matrix.h5
    """
}
