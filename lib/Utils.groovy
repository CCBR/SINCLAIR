class Utils {
    // run spooker for the workflow
    public static String spooker(workflow) {
        def pipeline_name = "${workflow.manifest.name.tokenize('/')[-1]}"
        def command_string = "spooker ${workflow.launchDir} ${pipeline_name}"
        def out = new StringBuilder()
        def err = new StringBuilder()
        try {
            def command = command_string.execute()
            command.consumeProcessOutput(out, err)
            command.waitFor()
        } catch(IOException e) {
            err = e
        }
        new FileWriter("${workflow.launchDir}/log/spooker.log").with {
            write("${out}\n${err}")
            flush()
        }
        return err
    }
}
