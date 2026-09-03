process SEURAT_PREPROCESS {
    tag "${id}"
    label 'process_medium'

    container "${params.container_seurat}"

    input:
    tuple val(id), val(inDir), path(h5), path(celldex_path)
    val(species)
    val(qc_filtering)
    val(nCount_RNA_max)
    val(nCount_RNA_min)
    val(nFeature_RNA_max)
    val(nFeature_RNA_min)
    val(percent_mt_max)
    val(percent_mt_min)
    val(doublet_finder)
    val(npcs)
    path(rmd)

    output:
    tuple val(id), path ("*seurat_preprocess.rds"), emit:rds
    tuple val(id), path ("*seurat_prefilter.rds"), emit:prefilter
    tuple val(id), path ("*seurat_preprocess.html"), emit:logs
 //   tuple val(id), path ("*seurat_preprocess_report.html"), emit:report

    script:
    def args = task.ext.args ?: ''
    """
    #!/usr/bin/env Rscript --vanilla
    rmarkdown::render("${rmd}",
        params=list(
            species="$species",
            sampleid="$id",
            h5="$h5",
            qc_filtering="$qc_filtering",

            nCount_RNA_max="$nCount_RNA_max",
            nCount_RNA_min="$nCount_RNA_min",
            nFeature_RNA_max="$nFeature_RNA_max",
            nFeature_RNA_min="$nFeature_RNA_min",
            percent_mt_max="$percent_mt_max",
            percent_mt_min="$percent_mt_min",
            doublet_finder="$doublet_finder",
            npcs="$npcs",
            celldex_cache="$celldex_path"
            ),
        output_file = "${id}_seurat_preprocess.html")

    """

    stub:
    """
    touch ${id}_seurat_prefilter.rds
    touch ${id}_seurat_preprocess.rds
    touch ${id}_seurat_preprocess.html
    """
}
