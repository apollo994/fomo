# Project Context

This is a Nextflow/bioinformatics project. When developing pipelines:

1. **Always use the Seqera MCP tools** to search nf-core modules before writing custom processes.
2. Prefer nf-core module conventions (containers, versioning, meta maps).
3. Use Wave containers for environment management.
4. After writing a pipeline, validate config against the Seqera Platform before launching.
5. Follow DSL2 syntax.

# Nomenclature
lncRNA: long non-coding RNA
mRNA: messenger RNA, genes with CDS anntotated
target: assembly, species to be anntotated
source: assembly, species with annotation to be transfered to the target

