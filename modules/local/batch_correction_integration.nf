process BATCH_CORRECT_INTEGRATION {
    tag "${gid}"
    label 'process_high'

    input:
    tuple val(gid), path(rds_m)
    tuple val(gid), path(rds_h)
    tuple val(gid), path(rds_r)
    tuple val(gid), path(rds_c)
    tuple val(gid), path(rds_s)
    tuple val(gid), path(rds_l)
    val(species)
    val(npcs)
    val(resolution_list)
    path(Rlib_dir)
    path(Rpkg_config)
    path(rmd)
    path(scRNA_functions)

    output:
    tuple val(gid), path ("*.html")                 , emit:logs

    script:
    def args = task.ext.args ?: ''
    """
    Rscript -e 'rmarkdown::render("${rmd}",
        params=list(gid="$gid",
            mergedObj="$rds_m",
            harmonyObj="$rds_h",
            ccaObj="$rds_c",
            rpcaObj="$rds_r",
            scviObj="$rds_s",
            ligerObj="$rds_l",
            npcs="$npcs",
            resolution_list="$resolution_list",
            citeseq="",
            annot="",
            Rlib_dir="$Rlib_dir",
            Rpkg_config="$Rpkg_config",
            scRNA_functions="$scRNA_functions",
            testing="N"),
        output_file = "${gid}_batch_correction_integration.html")'
    """

    stub:
    """
    touch ${gid}_batch_correction_integration.html
    """
}
