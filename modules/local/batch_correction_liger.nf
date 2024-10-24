process BATCH_CORRECT_LIGER {
    tag "${gid}"
    label 'process_high'

    container "${params.containers.scvi}"

    input:
    tuple val(gid), path(mergedObj)
    val(species)
    val(npcs)
    val(vars_to_regress)
    val(resolution_list)
    val(Rlib_dir)
    path(Rpkg_config)
    path(rmd)
    path(scRNA_functions)

    output:
    tuple val(gid), path ("*.rds")                 , emit:rds
    tuple val(gid), path ("*.html")                 , emit:logs

    script:
    def args = task.ext.args ?: ''
    """
    Rscript -e 'rmarkdown::render("${rmd}",
        params=list(gid="$gid",
            mergedObj="$mergedObj",
            species="$species",
            npcs="$npcs",
            vars_to_regress="$vars_to_regress",
            resolution_list="$resolution_list",
            Rlib_dir="$Rlib_dir",
            Rpkg_config="$Rpkg_config",
            scRNA_functions="$scRNA_functions",
            testing="N"),
        output_file = "${gid}_batch_correction_liger.html")'
    """

    stub:
    """
    touch ${gid}_batch_correction_liger.rds
    touch ${gid}_batch_correction_liger.html
    """
}
