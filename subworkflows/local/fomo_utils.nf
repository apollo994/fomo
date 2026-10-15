// Software versions as MultiQC's `*_mqc_versions.yml` ({process: {tool: version}}), from the
// [process, tool, version] items of the versions topic plus the workflow's own versions.
def softwareVersionsYaml(ch_versions) {
    def workflow_versions = channel.of(
        ['Workflow', 'Nextflow', "${nextflow.version}".toString()],
        ['Workflow', 'fomo', "${workflow.manifest.version}".toString()]
    )
    return ch_versions
        .map { proc, tool, version -> [proc.tokenize(':').last(), tool, "${version}".trim()] }
        .mix(workflow_versions)
        .unique()
        .toSortedList { a, b -> a[0] <=> b[0] ?: a[1] <=> b[1] }
        .map { rows ->
            def tree = [:]
            rows.each { proc, tool, version -> tree.get(proc, [:]).put(tool, version) }
            new org.yaml.snakeyaml.Yaml().dumpAsMap(tree)
        }
        .collectFile(name: 'fomo_software_mqc_versions.yml')
}
