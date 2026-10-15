// Helpers for the nf-test pipeline scenarios (tests/pipeline/*.nf.test).
class FomoTest {

    // Copy a repository samplesheet into dir with absolute input paths (nf-test runs elsewhere).
    static String sheet(String baseDir, String relSheet, String dir) {
        def out = new File(dir, 'samplesheet.csv')
        out.parentFile.mkdirs()
        out.text = new File(baseDir, relSheet).text.replace(',assets/', ",${baseDir}/assets/")
        return out.toString()
    }

    // Number of tasks per process (last component of the fully qualified name).
    static Map<String, Integer> taskCounts(List tasks) {
        return tasks.collect { t -> t.name.tokenize(' ')[0].tokenize(':')[-1] }.countBy { n -> n }.sort() as Map
    }

    // Published files, relative to outdir, without reports and pipeline_info.
    static List<String> published(String outdir) {
        def root = new File(outdir).toPath()
        def files = []
        root.toFile().eachFileRecurse(groovy.io.FileType.FILES) { f ->
            def rel = root.relativize(f.toPath()).toString()
            if (!rel.startsWith('pipeline_info/') && !rel.split('/').contains('multiqc')) {
                files << rel
            }
        }
        return files.sort()
    }

    static List<String> expected(String baseDir, String name) {
        return new File(baseDir, "tests/pipeline/expected/${name}.paths").readLines().findAll { l -> l }
    }
}
