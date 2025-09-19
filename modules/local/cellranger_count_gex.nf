process CELLRANGER_COUNT {
    tag "${id}"
    label 'process_high'

    input:
    tuple val(id), val(inDir)
    val(genome_dir)
    val(save_cellranger_extra_files)

    output:
    tuple val(id), path('*/outs/filtered_feature_bc_matrix.h5'), emit: h5
    tuple val(id), path('*/outs/*.bam*'),                        emit: bam_bai, optional: true
    tuple val(id), path('*/outs/*.cloupe'),                      emit: cloupe, optional: true

    script:
    def args = task.ext.args ?: ''
    def localmem = task.memory.toGiga()
    // if not saving extra files, clean up bam, bai, & cloupe after run
    def cleanup = save_cellranger_extra_files ? "" : "rm -rf $id/outs/*.bam* $id/outs/*.cloupe"
    """
    cellranger count \
        --id $id \
        --fastqs $inDir \
        --transcriptome=${genome_dir} \
        --localcores=$task.cpus \
        --localmem=$localmem
    ${cleanup}
    """

    stub:
    """
    mkdir -p $id/outs
    touch $id/outs/filtered_feature_bc_matrix.h5
    """
}
