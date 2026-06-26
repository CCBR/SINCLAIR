process COPY_DIR {
    // copy directory
    // -r recursively
    // -L, --dereference follow symbolic links
    input:
    path(input_dir)

    output:
    path("dir/")

    script:
    """
    cp -Lr ${input_dir} dir/
    """
}
