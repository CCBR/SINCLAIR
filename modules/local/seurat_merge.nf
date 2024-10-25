process SEURAT_MERGE {
    tag "${gid}"
    label 'process_high'

    container "${params.containers.seurat_merge}"

    input:
    tuple val(gid), path(rdsFiles)
    path(samplesheet)
    val(species)
    val(npcs)
    val(vars_to_regress)
    path(rmd)
    path(scRNA_functions)

    output:
    tuple val(gid), path ("*_seurat_merged.rds")                 , emit:rds
    tuple val(gid), path ("*_seurat_merged.pdf")                 , emit:logs

    script:
    def args = task.ext.args ?: ''
    """
    Rscript -e 'rmarkdown::render("${rmd}",
        params=list(species="$species",
            npcs="$npcs",
            vars_to_regress="$vars_to_regress",
            rdsFiles="$rdsFiles",
            gid="$gid",
            samplesheet="$samplesheet",
            scRNA_functions="$scRNA_functions",
            testing="N"),
        output_file = "${gid}_seurat_merged.pdf")'
    """

    stub:
    """
    echo $rdsFiles > ${gid}_seurat_merged.rds
    touch ${gid}_seurat_merged.pdf
    """
}
