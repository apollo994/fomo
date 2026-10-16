// A track is a feature type plus the decoy flag. Its keys, stamped on meta once: `track` (the
// file-name segment), `label` (report label and FASTA header type), `order` (row order in the
// reports) and `td2` (whether TD2 filters it).
def trackOf(String featureType, boolean decoy) {
    def labels = [lnc_RNA: 'lncRNA', mRNA: 'mRNA']
    if (!labels.containsKey(featureType)) {
        error("Unknown feature type '${featureType}': expected one of ${labels.keySet()}")
    }
    return [
        track: decoy ? "${featureType}.decoy".toString() : featureType,
        label: decoy ? "decoy_${labels[featureType]}".toString() : labels[featureType],
        order: (featureType == 'lnc_RNA' ? 1 : 3) + (decoy ? 1 : 0),
        td2:   featureType == 'lnc_RNA' && !decoy,
    ]
}

def withTrack(meta) {
    return meta + trackOf(meta.feature_type, meta.decoy as boolean)
}

// Pairs each [meta, ...] item of `left` with every [key, value] item of `right` whose key is
// key(meta): [meta, ..., value]. A left item without a partner is dropped.
def pairBy(left, right, key) {
    return left
        .map { item -> [key.call(item[0]), item] }
        .combine(right, by: 0)
        .map { _key, item, value -> item.toList() + [value] }
}

def pairById(left, right) {
    return pairBy(left, right, { meta -> meta.id })
}

// [target, [files]] for every target: an item whose meta has a target_id goes to that target, the
// others to every target. Items are flattened first: combine spreads a list across the tuple.
def routeByTarget(ch_items, ch_target_ids) {
    def routed = ch_items.branch { meta, _p ->
        scoped: meta.target_id != null
        shared: true
    }
    return routed.shared
        .flatMap { _meta, p -> [p].flatten() }
        .combine(ch_target_ids)
        .map { p, tid -> tuple(tid, p) }
        .mix(routed.scoped.flatMap { meta, p -> [p].flatten().collect { f -> tuple(meta.target_id, f) } })
        .groupTuple()
}

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
