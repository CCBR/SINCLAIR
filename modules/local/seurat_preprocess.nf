process SEURAT_PREPROCESS {
    tag "${id}"
    label 'process_medium'

    container "${params.containers.seurat_preproc}"

    input:
    tuple val(id), val(inDir), path(h5)
    val(species)
    val(qc_filtering)
    val(nCount_RNA_max)
    val(nCount_RNA_min)
    val(nFeature_RNA_max)
    val(nFeature_RNA_min)
    val(percent_mt_max)
    val(percent_mt_min)
    val(run_doublet_finder)
    val(npcs)
    val(Rlib_dir)
    path(Rpkg_config)
    path(rmd)
    path(scRNA_functions)

    output:
    tuple val(id), path ("*.rds")                 , emit:rds
    tuple val(id), path ("*.pdf")                 , emit:logs

    script:
    def args = task.ext.args ?: ''
    """
    Rscript -e 'rmarkdown::render("${rmd}",
        params=list(
            species="$species",
            sampleid="$id",
            h5="$h5",
            qc_filtering="$qc_filtering",
            nCount_RNA_max=$nCount_RNA_max,
            nCount_RNA_min=$nCount_RNA_min,
            nFeature_RNA_max=$nFeature_RNA_max,
            nFeature_RNA_min=$nFeature_RNA_min,
            percent_mt_max=$percent_mt_max,
            percent_mt_min=$percent_mt_min,
            run_doublet_finder="$run_doublet_finder",
            npcs=$npcs,
            Rlib_dir="$Rlib_dir",
            Rpkg_config="$Rpkg_config",
            scRNA_functions="$scRNA_functions"),
        output_file = "${id}_seurat_preprocess.pdf")'
    """

    stub:
    """
    touch ${id}_seurat_preprocess.rds
    touch ${id}_seurat_preprocess.pdf
    """
}
