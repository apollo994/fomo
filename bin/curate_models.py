#!/usr/bin/env python3
"""Curate fomo projected models: report intron chains shared exactly by >= N species.

Input:  a fomo GFF3 (<target>.allModels.<track>.raw.gff3 or any per-source subset; .gz ok)
        with gene / transcript / exon rows and source=<species> on the transcript, and optionally
        the target's reference annotation (GFF3 or GTF, .gz ok).

Method
  1. Drop the target's own projections (never support), unstranded and single-exon models.
  2. Reference filter (with --ref): drop every model with an exon overlapping any exon entry of the
     reference on the same strand (--strand both: either strand), whatever its type (mRNA,
     lnc_RNA, antisense_RNA, tRNA, misc_RNA, ...). Exon-less transcript-like features count by
     their CDS/UTR rows or, with no children at all, by their own span.
  3. Exact intron chains: one GTF per species -> `gffcompare --no-merge` -> .tracking cliques,
     checked against the chains derived here (exit 1 on any disagreement), expanded back to every
     model with that chain (gffcompare drops same-species duplicates). A chain is supported when
     its models (after the filter) come from >= --min-species species.
  4. Genes: supported chains sharing >= 1 splice junction (identical intron), connected components.
  5. One transcript per supported chain of the gene. Its exons are those of the chain's
     representative: the longest exonic length among the chain's models (after the filter) -> lowest
     de -> ID. Transcripts are numbered t1, t2, ... by support: most species -> most SJs -> longest
     -> leftmost start -> rightmost end -> representative ID. Each lists the IDs of every model with
     its chain (supporting_transcripts) and their species.
  6. Location relative to the reference (ref_location on the gene; first match wins), computed on the
     gene's footprint (union of its transcripts' exons):
     antisense_exonic (exon overlaps an opposite-strand reference exon) > intronic (span inside one
     same-strand reference transcript span) > sense_span_overlap (overlaps a same-strand reference
     transcript span, not contained) > antisense_intronic (overlaps an opposite-strand reference
     transcript span, no exon overlap) > divergent (no overlap; head-to-head with an opposite-strand
     reference transcript, 5' ends <= --divergent-dist apart) > intergenic; no_reference without
     --ref. A reference transcript span is its first-to-last exon (CDS/UTR or own span if exon-less);
     unstranded reference rows count on both strands.
  7. Reference gene tag (with --ref), on the gene: ref_gene_id / ref_gene_name / ref_gene_biotype /
     ref_gene_orientation (sense|antisense) / ref_distance (0 = overlapping) — the gene the
     ref_location class refers to (overlapping or divergent partner), else the nearest reference gene
     on either strand. Keys without a value are omitted (e.g. no ref_* key at all when the seqid has
     no reference transcript; no ref_gene_name for an unnamed gene).

IDs: H = first 12 hex of sha1('<track>|<seqid>|<strand>|<sorted intron chains>') of the gene; gene
<P>_g_H, transcript <P>_t<n>_H, exon <P>_t<n>e<k>_H with k numbered 5'->3', where P = --id-prefix
(default fomo_<tag>, tag = lncRNA | mRNA | decoy_lncRNA | decoy_mRNA, also the biotype). The same
chains on the same assembly always give the same IDs, whatever the source set. Transcripts are
typed by the track (lnc_RNA | mRNA; decoys as their real track); alignment metrics are not copied
(they stay in the per-source projected GFF3s, keyed by supporting_transcripts).

Output (--prefix P)
  P.gff3.gz           curated genes (gene + one transcript per supported chain + exons)
  P.merged.gff3.gz    the reference verbatim + curated genes, each after its seqid's last
                      reference line (with --ref only)
  P.curation.tsv      report: table / metric / value (observed, discarded_by_reference, novel)
  P.curation.json     the same numbers
  P.curation.data.json  the same numbers as a `curation` data record (with --identity)

Rationale: plans/21, plans/25.
"""

import argparse
import bisect
import collections
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from array import array
from dataclasses import dataclass, field
from pathlib import Path
from typing import NamedTuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fomo_io  # noqa: E402
import fomo_stats  # noqa: E402

GENERIC_TYPES = {'transcript', 'ncRNA', 'primary_transcript', 'RNA', 'gene', 'pseudogene', 'ncRNA_gene'}
BIOTYPE_KEYS = (
    'transcript_biotype',
    'transcript_type',
    'biotype',
    'ncrna_class',
    'gene_biotype',
    'gene_type',
    'gbkey',
)
SMALL_RNA_MARKERS = (
    'trna',
    'rrna',
    'snorna',
    'snrna',
    'mirna',
    'scrna',
    'y_rna',
    'srp',
    'rnase',
    'guide',
    'ncrna',
    'misc',
)
ESCAPE = re.compile(r'[;=&,\t\n\r]')
MERGED_ID_RE = re.compile(r'(?:^|;)\s*ID=([^;\t\n]+)')
HASH_LENGTH = 12  # hex digits of the gene content hash in every curated ID
MAX_PARENT_DEPTH = 10  # Parent hops followed to find a reference feature's gene
ORPHAN_EXON = '<orphan exon>'  # interval key of a reference exon without Parent
LOCATIONS = fomo_stats.CURATION_LOCATIONS
OVERLAPPING = fomo_stats.CURATION_OVERLAPPING
REPORT_TABLES = ('observed', 'discarded_by_reference', 'novel')
log = fomo_io.get_logger('curate_models')


def esc(value) -> str:
    return ESCAPE.sub(lambda c: f'%{ord(c.group(0)):02X}', str(value))


def gtf_attrs(col9: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for key, value in re.findall(r'(\S+)\s+"([^"]*)"', col9):
        attrs.setdefault(key, value)
    return attrs


def is_gtf(path: str) -> bool:
    with fomo_io.open_text(path) as fh:
        for line in fh:
            if line[0] != '#':
                cols = line.rstrip('\n').split('\t')
                if len(cols) >= 9:
                    return '"' in cols[8] and '=' not in cols[8].split('"')[0]
    return False


def track_info(track: str | None) -> tuple[str, str]:
    """-> (tag, transcript feature type): ('lncRNA' | 'mRNA' | 'decoy_lncRNA' | 'decoy_mRNA', 'lnc_RNA' | 'mRNA').
    The tag is the biotype and the default ID prefix token; decoys keep their real track's feature type."""
    name = track or ''
    mrna = bool(re.search(r'(^|[._])mRNA($|[._])', name))
    tag = ('decoy_' if 'decoy' in name else '') + ('mRNA' if mrna else 'lncRNA')
    return tag, 'mRNA' if mrna else 'lnc_RNA'


# ------------------------------------------------------------------ models
class ChainKey(NamedTuple):
    chr: str
    strand: str
    introns: tuple[tuple[int, int], ...]


@dataclass(slots=True)
class Model:
    """One projected model. Its attributes are not kept: the curated output carries none of them."""

    id: str
    chr: str
    strand: str
    exons: tuple[tuple[int, int], ...]
    species: str
    de: float
    chain: tuple[tuple[int, int], ...]
    ref_hits: set = field(default_factory=set)

    @property
    def key(self) -> ChainKey:
        return ChainKey(self.chr, self.strand, self.chain)

    @property
    def exonic_len(self) -> int:
        return sum(e - s + 1 for s, e in self.exons)


def read_models(path: str) -> list[Model]:
    """fomo GFF3 -> models in transcript order. Exits 1 on structural problems."""
    transcripts: dict[str, tuple[str, str, dict[str, str]]] = {}  # ID -> (seqid, strand, attributes)
    exons: dict[str, list] = {}  # transcript ID -> [seqid, strand, [(start, end)]]
    for lineno, cols in fomo_io.iter_gff(path, strict=True):
        chrom = sys.intern(cols[0])
        attrs = fomo_io.parse_attrs(cols[8])
        if cols[2].lower() == 'exon':
            if 'Parent' not in attrs:
                fomo_io.die(f'{path}:{lineno}: exon without Parent')
            for parent in fomo_io.split_ids(attrs['Parent']):
                entry = exons.get(parent)
                if entry is None:
                    exons[parent] = [chrom, cols[6], [(int(cols[3]), int(cols[4]))]]
                elif entry[0] is not chrom or entry[1] != cols[6]:
                    fomo_io.die(f'transcript {parent}: exons on several seqids/strands')
                else:
                    entry[2].append((int(cols[3]), int(cols[4])))
        elif fomo_io.is_transcript(cols[2]) and attrs.get('ID'):
            if attrs['ID'] in transcripts:
                fomo_io.die(f'{path}:{lineno}: duplicate transcript ID {attrs["ID"]}')
            transcripts[attrs['ID']] = (chrom, cols[6], attrs)
    missing = set(exons) - set(transcripts)
    if missing:
        fomo_io.die(f'{len(missing)} exon parent(s) not defined as transcripts, e.g. {sorted(missing)[0]}')
    return [model_of(tid, row, exons[tid]) for tid, row in transcripts.items() if exons.get(tid)]


def model_of(tid: str, row: tuple, exon_entry: list) -> Model:
    chrom, strand, attrs = row
    if (exon_entry[0], exon_entry[1]) != (chrom, strand):
        fomo_io.die(f'transcript {tid}: exons on {exon_entry[0]}{exon_entry[1]}, transcript row on {chrom}{strand}')
    exs = tuple(sorted(exon_entry[2]))
    for (s1, e1), (s2, e2) in zip(exs, exs[1:], strict=False):
        if s2 <= e1 + 1:
            fomo_io.die(f'transcript {tid}: overlapping or abutting exons {s1}-{e1} / {s2}-{e2}')
    species = attrs['source'] if 'source' in attrs else (tid.rsplit('|', 1)[1] if '|' in tid else '')
    if not species:
        fomo_io.die(f'transcript {tid}: no source= attribute and no |<source> ID suffix')
    try:
        de = float(attrs['de']) if 'de' in attrs else math.nan
    except ValueError:
        de = math.nan
    chain = tuple((exs[i][1], exs[i + 1][0]) for i in range(len(exs) - 1))
    return Model(tid, chrom, strand, exs, sys.intern(species), math.inf if math.isnan(de) else de, chain)


# ------------------------------------------------------------------ reference
class RefIndex:
    """(chr, strand) -> intervals sorted by start, in typed arrays, each with a label (any hashable,
    stored once); hits() returns every label overlapping [a, b] (>= 1 bp).

    Built with add(…, key) where key is a small int chosen by the caller (e.g. the index of the
    parent transcript); freeze(labels) maps every key to its label in one pass, so the labels —
    which need the whole reference read first (gene of a transcript) — cost no per-interval memory."""

    def __init__(self):
        self.labels = []
        self.build = collections.defaultdict(lambda: (array('l'), array('l'), array('l')))
        self.idx = {}

    def add(self, chrom, strand, start, end, key):
        starts, ends, keys = self.build[(chrom, strand)]
        starts.append(start)
        ends.append(end)
        keys.append(key)

    def freeze(self, key_labels):
        """key_labels: list, key -> label."""
        label_index, key_label = {}, []
        for label in key_labels:
            li = label_index.get(label)
            if li is None:
                li = label_index[label] = len(self.labels)
                self.labels.append(label)
            key_label.append(li)
        for chrom_strand, (starts, ends, keys) in self.build.items():
            order = sorted(range(len(starts)), key=starts.__getitem__)
            st = array('l', (starts[k] for k in order))
            en = array('l', (ends[k] for k in order))
            lb = array('l', (key_label[keys[k]] for k in order))
            self.idx[chrom_strand] = (st, en, lb, max(e - s for s, e in zip(st, en, strict=True)) + 1)
        self.build = None
        return self

    def intervals(self, chrom, strand, a, b):
        """-> [(start, end, label)] overlapping [a, b]."""
        x = self.idx.get((chrom, strand))
        out = []
        if x:
            st, en, lb, max_len = x
            for i in range(bisect.bisect_left(st, a - max_len), bisect.bisect_right(st, b)):
                if en[i] >= a:
                    out.append((st[i], en[i], self.labels[lb[i]]))
        return out

    def hits(self, chrom, strand, a, b):
        return {label for _, _, label in self.intervals(chrom, strand, a, b)}

    def nearest(self, chrom, strand, a, b):
        """-> (distance, label) of the closest interval to [a, b] (0 if overlapping), or None."""
        x = self.idx.get((chrom, strand))
        if not x:
            return None
        st, en, lb, max_len = x
        overlapping = self.intervals(chrom, strand, a, b)
        if overlapping:
            return (0, min(label for _, _, label in overlapping))
        best = None
        j = bisect.bisect_right(st, b)  # first interval starting right of b
        if j < len(st):
            best = (st[j] - b, self.labels[lb[j]])
        i = bisect.bisect_left(st, a) - 1  # intervals starting left of a: max end
        max_end, label = None, None
        while i >= 0:
            if max_end is not None and st[i] + max_len <= max_end:
                break
            if max_end is None or en[i] > max_end:
                max_end, label = en[i], self.labels[lb[i]]
            i -= 1
        if max_end is not None and (best is None or a - max_end < best[0]):
            best = (a - max_end, label)
        return best


class _RefBuilder:
    """Collects the reference while it is read: exon intervals keyed by their parent, per-parent
    spans, the CDS/UTR and childless fallbacks, and what resolves a parent to its gene."""

    def __init__(self):
        self.feature_label = {}  # feature id -> type label
        self.parent_of = {}  # feature id -> first Parent (GFF3) / gene_id (GTF)
        self.name = {}  # id -> Name / gene_name / gene
        self.alt_name = {}  # transcript id -> gene_name / gene
        self.biotype = {}  # id -> gene biotype
        self.keys, self.key_ids = {}, []  # parent id -> interval key, and back
        self.cds = collections.defaultdict(list)  # parent -> CDS/UTR intervals (used if it has no exon)
        self.with_exons, self.has_child, self.parented = set(), set(), []
        self.spans = {}  # (key, chrom, strand) -> (min, max)
        self.ex_idx = RefIndex()
        self.origin = collections.Counter()

    @staticmethod
    def label(ftype: str, attrs: dict[str, str]) -> str:
        """The type, refined by the biotype for generic types (e.g. transcript:misc_RNA)."""
        if ftype in GENERIC_TYPES:
            biotype = next((attrs[k] for k in BIOTYPE_KEYS if attrs.get(k)), None)
            if biotype and biotype != ftype:
                return sys.intern(f'{ftype}:{biotype}')
        return sys.intern(ftype)

    def _key(self, pid: str) -> int:
        key = self.keys.get(pid)
        if key is None:
            key = self.keys[pid] = len(self.key_ids)
            self.key_ids.append(pid)
        return key

    def add_interval(self, pid, chrom, strands, start, end):
        key = self._key(pid)
        for strand in strands:
            self.ex_idx.add(chrom, strand, start, end, key)
            span_key = (key, chrom, strand)
            span = self.spans.get(span_key)
            self.spans[span_key] = (start, end) if span is None else (min(span[0], start), max(span[1], end))

    def add_gtf_row(self, cols, chrom, strands, start, end):
        attrs = gtf_attrs(cols[8])
        tid = attrs.get('transcript_id')
        if not tid:
            return
        ftype = cols[2].lower()
        if ftype == 'transcript' or (ftype == 'exon' and tid not in self.feature_label):
            self.feature_label[tid] = self.label('transcript', attrs)
            gid = attrs.get('gene_id') or tid
            self.parent_of[tid] = gid
            if attrs.get('gene_name'):
                self.name.setdefault(gid, attrs['gene_name'])
            biotype = attrs.get('gene_biotype') or attrs.get('gene_type')
            if biotype:
                self.biotype.setdefault(gid, biotype)
        if ftype == 'exon':
            self.with_exons.add(tid)
            self.add_interval(tid, chrom, strands, start, end)
            self.origin['exon'] += 1
        elif ftype in fomo_io.CDS_LIKE_TYPES:
            self.cds[tid].append((chrom, strands, start, end))

    def add_gff3_row(self, cols, chrom, strands, start, end):
        attrs = fomo_io.parse_attrs(cols[8])
        parents = fomo_io.split_ids(attrs.get('Parent', ''))
        ftype = cols[2].lower()
        if ftype == 'exon':
            self.origin['exon'] += 1
            for parent in parents or [ORPHAN_EXON]:
                self.with_exons.add(parent)
                self.add_interval(parent, chrom, strands, start, end)
            return
        if ftype in fomo_io.CDS_LIKE_TYPES:
            for parent in parents:
                self.cds[parent].append((chrom, strands, start, end))
            return
        fid = attrs.get('ID')
        if fid:
            self.describe(fid, cols[2], attrs, parents)
        self.has_child.update(parents)
        if parents and fid and fomo_io.is_transcript(cols[2]):
            self.parented.append((fid, chrom, strands, start, end))

    def describe(self, fid, ftype, attrs, parents):
        self.feature_label[fid] = self.label(ftype, attrs)
        if parents:
            self.parent_of[fid] = parents[0]
        name = attrs.get('Name') or attrs.get('gene_name') or attrs.get('gene')
        if name:
            self.name[fid] = name
        alt = attrs.get('gene_name') or attrs.get('gene')
        if alt:
            self.alt_name[fid] = alt
        biotype = attrs.get('gene_biotype') or attrs.get('biotype') or attrs.get('gene_type')
        if biotype:
            self.biotype[fid] = biotype

    def add_fallback_intervals(self):
        """A transcript-like feature without exons counts by its CDS/UTR rows or, with no children
        at all, by its own span."""
        for parent, rows in self.cds.items():
            if parent in self.with_exons:
                continue
            for chrom, strands, start, end in rows:
                self.add_interval(parent, chrom, strands, start, end)
                self.origin['cds_utr_without_exon'] += 1
        for fid, chrom, strands, start, end in self.parented:
            if fid in self.has_child or fid in self.with_exons or fid in self.cds:
                continue
            self.add_interval(fid, chrom, strands, start, end)
            self.origin['childless_feature_span'] += 1

    def gene_of(self, pid: str) -> str:
        """The top of pid's Parent chain (GTF: its gene_id)."""
        gene, depth = pid, 0
        while gene in self.parent_of and depth < MAX_PARENT_DEPTH:
            gene = self.parent_of[gene]
            depth += 1
        return gene

    def resolve_genes(self):
        """-> (exon RefIndex, transcript-span RefIndex, gene_meta), every interval labelled
        (type, gene); gene_meta: gene -> (name, biotype)."""
        gene_meta, key_labels = {}, []
        for pid in self.key_ids:
            gene = self.gene_of(pid)
            ftype = self.feature_label.get(pid, 'unknown_parent')
            if gene not in gene_meta:
                name = self.name.get(gene) if gene != pid else None  # a transcript's own Name is not a gene name
                gene_meta[gene] = (
                    name or self.alt_name.get(pid) or self.name.get(pid) or '.',
                    self.biotype.get(gene) or ftype,
                )
            key_labels.append((ftype, gene))
        tx_idx = RefIndex()
        for (key, chrom, strand), (lo, hi) in self.spans.items():
            tx_idx.add(chrom, strand, lo, hi, key)
        return self.ex_idx.freeze(key_labels), tx_idx.freeze(key_labels), gene_meta


def read_reference(path: str):
    """-> (exon RefIndex, transcript-span RefIndex, Counter of interval origins, is GTF, gene_meta).

    Every exon row counts. A transcript-like feature without exons counts by its CDS/UTR rows or,
    with no children at all, by its own span. Unstranded rows ('.', '?') count on both strands.
    Both indexes carry the label (type, gene): `type` is the exon's parent type (generic types
    refined by biotype, e.g. transcript:misc_RNA) and `gene` the top of its Parent chain (GFF3) or
    its gene_id (GTF). gene_meta: gene -> (name, biotype); name = the gene's Name / gene_name /
    gene, else the transcript's gene_name / gene, else '.'; biotype = the gene's gene_biotype /
    biotype / gene_type, else the transcript's type label."""
    gtf = is_gtf(path)
    builder = _RefBuilder()
    add_row = builder.add_gtf_row if gtf else builder.add_gff3_row
    for _, cols in fomo_io.iter_gff(path):
        strands = (cols[6],) if cols[6] in ('+', '-') else ('+', '-')  # '.' / '?' -> both
        add_row(cols, sys.intern(cols[0]), strands, int(cols[3]), int(cols[4]))
    builder.add_fallback_intervals()
    ex_idx, tx_idx, gene_meta = builder.resolve_genes()
    return ex_idx, tx_idx, builder.origin, gtf, gene_meta


# ------------------------------------------------------------------ gffcompare
def run_gffcompare(models: list[Model], tmp: str, exe: str):
    """-> (cliques as lists of model indices, number of input files); writes one GTF per species
    with integer IDs."""
    by_species = collections.defaultdict(list)
    for i, m in enumerate(models):
        by_species[m.species].append(i)
    files = []
    for k, species in enumerate(sorted(by_species)):
        path = os.path.join(tmp, f'q{k}.gtf')
        with open(path, 'w', encoding='utf-8') as out:
            for i in by_species[species]:
                m = models[i]
                tag = f'gene_id "g{i}"; transcript_id "t{i}";'
                out.write(f'{m.chr}\tfomo\ttranscript\t{m.exons[0][0]}\t{m.exons[-1][1]}\t.\t{m.strand}\t.\t{tag}\n')
                for s, e in m.exons:
                    out.write(f'{m.chr}\tfomo\texon\t{s}\t{e}\t.\t{m.strand}\t.\t{tag}\n')
        files.append(path)
    with open(os.path.join(tmp, 'inputs.txt'), 'w', encoding='utf-8') as out:
        out.write('\n'.join(files) + '\n')
    r = subprocess.run([exe, '--no-merge', '-o', 'cmp', '-i', 'inputs.txt'], cwd=tmp, capture_output=True, check=False)
    if r.returncode != 0:
        fomo_io.die(f'gffcompare failed (exit {r.returncode}):\n{r.stderr.decode("utf-8", "replace")}')
    cliques = []
    with open(os.path.join(tmp, 'cmp.tracking'), encoding='utf-8') as fh:
        for line in fh:
            ids = []
            for query in line.rstrip('\n').split('\t')[4:]:
                if query == '-':
                    continue
                for item in query.split(','):  # qN:g<i>|t<i>|exons|fpkm|tpm|cov|len
                    ids.append(int(item.split(':', 1)[1].split('|')[1][1:]))
            cliques.append(ids)
    return cliques, len(files)


def exact_chains(models: list[Model], tmp: str, exe: str):
    """-> ({chain key: [models]} as decided by gffcompare, verified against the chains derived here,
    number of gffcompare inputs, models restored that gffcompare dropped as duplicates)."""
    if not models:
        return {}, 0, 0
    cliques, n_files = run_gffcompare(models, tmp, exe)
    by_key = collections.defaultdict(list)
    for m in models:
        by_key[m.key].append(m)
    seen = set()
    for ids in cliques:
        keys = {models[i].key for i in ids}
        if len(keys) != 1:
            fomo_io.die(
                f'gffcompare clique with {len(keys)} different intron chains (ids {ids[:5]}); '
                f'check gffcompare version/options'
            )
        key = keys.pop()
        if key in seen:
            fomo_io.die(f'intron chain {key[0]}:{key[1]} reported in two gffcompare cliques')
        seen.add(key)
    if seen != set(by_key):
        fomo_io.die(f'{len(set(by_key) - seen)} intron chain(s) missing from gffcompare tracking')
    restored = sum(len(by_key[k]) for k in seen) - sum(len(c) for c in cliques)
    return by_key, n_files, restored


# ------------------------------------------------------------------ genes / representative
def link_genes(chains) -> list[list[ChainKey]]:
    """chain keys -> groups of chain keys sharing >= 1 intron (connected components)."""
    chains = sorted(chains)
    parent = list(range(len(chains)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    first = {}
    for i, (chrom, strand, introns) in enumerate(chains):
        for intron in introns:
            k = (chrom, strand, intron)
            if k in first:
                parent[find(i)] = find(first[k])
            else:
                first[k] = i
    groups = collections.defaultdict(list)
    for i, key in enumerate(chains):
        groups[find(i)].append(key)
    return list(groups.values())


def qstats(vals) -> dict:
    if not vals:
        return {'n': 0}
    v = sorted(vals)

    def quantile(p):
        return v[min(len(v) - 1, int(round(p * (len(v) - 1))))]

    return {
        'n': len(v),
        'min': v[0],
        'p25': quantile(0.25),
        'median': statistics.median(v),
        'mean': round(statistics.fmean(v), 2),
        'p75': quantile(0.75),
        'max': v[-1],
    }


def type_rank(ftype: str) -> int:
    lower = ftype.lower()
    if ftype == 'mRNA' or 'protein_coding' in lower:
        return 0
    if ftype in ('lnc_RNA', 'lncRNA') or 'lncrna' in lower or 'lincrna' in lower:
        return 1
    if 'antisense' in lower:
        return 2
    if 'pseudo' in lower:
        return 5
    if ftype.startswith('transcript'):
        return 3
    if any(x in lower for x in SMALL_RNA_MARKERS):
        return 4
    return 6


def primary_type(types) -> str:
    return min(types, key=lambda t: (type_rank(t), t))


def ref_location(gene, ex_idx, tx_idx, gene_meta, divergent_dist):
    """-> (class, ref) of a gene footprint (chr, strand, sorted disjoint exons) vs the reference.
    ref = dict(gene, name, biotype, orientation, distance) for the reference gene the class refers
    to: the gene of the overlapping opposite-strand exon (antisense_exonic), of the innermost
    containing same-strand transcript (intronic), of an overlapping transcript (sense_span_overlap /
    antisense_intronic), the head-to-head partner (divergent), or the nearest gene on either strand
    (intergenic; ties: sense, then gene id). None when the seqid has no reference transcript (or
    without a reference). Distance: 0 when overlapping, else start/end coordinate difference
    (adjacent = 1)."""
    if ex_idx is None:
        return 'no_reference', None
    chrom, strand, exons = gene['chr'], gene['strand'], gene['exons']
    start, end = exons[0][0], exons[-1][1]
    opp = '-' if strand == '+' else '+'

    def tag(g, orientation, distance):
        name, biotype = gene_meta[g]
        return {'gene': g, 'name': name, 'biotype': biotype, 'orientation': orientation, 'distance': distance}

    anti = {g for s, e in exons for _, g in ex_idx.hits(chrom, opp, s, e)}
    if anti:
        return 'antisense_exonic', tag(min(anti), 'antisense', 0)
    same = tx_idx.intervals(chrom, strand, start, end)
    containing = [(e2 - s2, g) for s2, e2, (_, g) in same if s2 <= start and e2 >= end]
    if containing:
        return 'intronic', tag(min(containing)[1], 'sense', 0)
    if same:
        return 'sense_span_overlap', tag(min(g for _, _, (_, g) in same), 'sense', 0)
    opposite = tx_idx.intervals(chrom, opp, start, end)
    if opposite:
        return 'antisense_intronic', tag(min(g for _, _, (_, g) in opposite), 'antisense', 0)
    # head-to-head: an opposite-strand transcript whose 5' end lies within divergent_dist upstream of ours
    if divergent_dist > 0:
        window = (start - divergent_dist, start - 1) if strand == '+' else (end + 1, end + divergent_dist)
        partners = tx_idx.intervals(chrom, opp, *window)
        if partners:
            d, g = min(((start - e2) if strand == '+' else (s2 - end), g) for s2, e2, (_, g) in partners)
            return 'divergent', tag(g, 'antisense', d)
    near = []
    for st, orientation in ((strand, 0), (opp, 1)):
        n = tx_idx.nearest(chrom, st, start, end)
        if n:
            near.append((n[0], orientation, n[1][1]))
    if not near:
        return 'intergenic', None
    d, orientation, g = min(near)
    return 'intergenic', tag(g, 'sense' if orientation == 0 else 'antisense', d)


def check_location(loc, ref, divergent_dist, rid):
    """Fail if the tag contradicts the class (it cannot by construction; this guards changes)."""
    ok = True
    if loc in OVERLAPPING:
        expected = 'antisense' if loc.startswith('antisense') else 'sense'
        ok = ref is not None and ref['distance'] == 0 and ref['orientation'] == expected
    elif loc == 'divergent':
        ok = ref is not None and ref['orientation'] == 'antisense' and 0 < ref['distance'] <= divergent_dist
    elif loc == 'intergenic':
        ok = ref is None or ref['distance'] > 0
    if not ok:
        fomo_io.die(f'inconsistent reference tag for {rid}: ref_location={loc}, {ref}')


# ------------------------------------------------------------------ curated genes
@dataclass(slots=True)
class CuratedTranscript:
    rep: Model
    key: ChainKey
    species: list[str]
    support: list[Model]
    id: str = ''
    exon_prefix: str = ''


@dataclass(slots=True)
class CuratedGene:
    gene_id: str
    hash: str
    chr: str
    strand: str
    start: int
    end: int
    tx: list[CuratedTranscript]
    chains: list[ChainKey]
    species: list[str]
    location: str
    ref: dict | None

    @property
    def n_chains(self) -> int:
        return len(self.chains)


def merge_exons(exons):
    """-> sorted, disjoint union of (start, end) intervals (overlapping or abutting ones merged)."""
    out = []
    for s, e in sorted(exons):
        if out and s <= out[-1][1] + 1:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return tuple(out)


def gene_hash(track, chrom, strand, chains, n=HASH_LENGTH):
    """Content ID of a gene: its track, seqid, strand and sorted intron chains."""
    chain_text = ';'.join(sorted(','.join(f'{a}-{b}' for a, b in ch) for ch in chains))
    return hashlib.sha1(f'{track}|{chrom}|{strand}|{chain_text}'.encode()).hexdigest()[:n]


def ref_attrs(loc, ref):
    """-> [(key, value)] of the gene's reference tags; keys without a value are omitted."""
    attrs = [('ref_location', loc)]
    if loc == 'no_reference' or not ref:
        return attrs
    for key, value in (
        ('ref_gene_id', esc(ref['gene'])),
        ('ref_gene_name', esc(ref['name'])),
        ('ref_gene_biotype', esc(ref['biotype'])),
        ('ref_gene_orientation', ref['orientation']),
        ('ref_distance', str(ref['distance'])),
    ):
        if value not in ('', '.'):
            attrs.append((key, value))
    return attrs


def gff_line(chrom, ftype, start, end, strand, attrs):
    col9 = ';'.join(f'{k}={v}' for k, v in attrs)
    return '\t'.join([chrom, 'fomo', ftype, str(start), str(end), '.', strand, '.', col9])


def gene_blocks(genes: list[CuratedGene], tag: str, ftype: str):
    """-> list of (seqid, start, [lines]) per curated gene: gene, then t1 + its exons, t2 + ..."""
    blocks = []
    for g in genes:
        attrs = [('ID', g.gene_id), ('biotype', tag), ('n_chains', str(g.n_chains)), ('n_species', str(len(g.species)))]
        lines = [gff_line(g.chr, 'gene', g.start, g.end, g.strand, attrs + ref_attrs(g.location, g.ref))]
        for t in g.tx:
            r = t.rep
            attrs = [
                ('ID', t.id),
                ('Parent', g.gene_id),
                ('biotype', tag),
                ('representative_transcript', esc(r.id)),
                ('supporting_transcripts', ','.join(esc(m.id) for m in t.support)),
                ('supporting_species', ','.join(esc(x) for x in t.species)),
                ('n_sj', str(len(r.chain))),
            ]
            lines.append(gff_line(r.chr, ftype, r.exons[0][0], r.exons[-1][1], r.strand, attrs))
            n = len(r.exons)
            for j, (xs, xe) in enumerate(r.exons):
                k = j + 1 if r.strand == '+' else n - j  # exon number 5'->3'
                exon_attrs = [('ID', f'{t.exon_prefix}e{k}_{g.hash}'), ('Parent', t.id)]
                lines.append(gff_line(r.chr, 'exon', xs, xe, r.strand, exon_attrs))
        blocks.append((g.chr, g.start, lines))
    return blocks


def all_ids(genes: list[CuratedGene]) -> set[str]:
    """-> every gene / transcript / exon ID written for these genes."""
    ids = set()
    for g in genes:
        ids.add(g.gene_id)
        for t in g.tx:
            ids.add(t.id)
            ids.update(f'{t.exon_prefix}e{k}_{g.hash}' for k in range(1, len(t.rep.exons) + 1))
    return ids


def header(args, n_genes, n_tx):
    strand = 'same-strand' if args.strand == 'same' else 'either-strand'
    return [
        '##gff-version 3',
        f'# curate_models: intron chains shared exactly by >={args.min_species} species (gffcompare --no-merge), '
        f'linked into genes by shared splice junctions; one transcript per supported chain (t1 = most species); '
        f'reference filter: {strand} overlap with any reference exon; '
        f'{n_genes} genes, {n_tx} transcripts',
    ]


def write_curated(path, args, blocks, n_tx):
    with open(path, 'w', encoding='utf-8') as out:
        out.write('\n'.join(header(args, len(blocks), n_tx)) + '\n')
        for _, _, lines in blocks:
            out.write('\n'.join(lines) + '\n')


def reference_layout(ref_path, new_ids):
    """-> {seqid: index of its last reference line}. Fails if a curated ID is already used by the
    reference."""
    last = {}
    with fomo_io.open_text(ref_path) as fh:
        for n, line in enumerate(fh):
            if line.startswith('##FASTA'):
                break
            if line[0] != '#' and line.strip():
                cols = line.split('\t', 9)
                last[cols[0]] = n
                match = MERGED_ID_RE.search(cols[8]) if len(cols) >= 9 else None
                if match and match.group(1).strip() in new_ids:
                    fomo_io.die(f'curated ID {match.group(1).strip()} already used in the reference; set --id-prefix')
    return last


def write_merged(path, ref_path, args, blocks, new_ids, n_tx):
    """Reference verbatim; each seqid's curated blocks right after its last reference feature line."""
    last = reference_layout(ref_path, new_ids)
    by_seq = collections.defaultdict(list)
    for chrom, _, lines in blocks:
        by_seq[chrom].append(lines)
    rest = [c for c in dict.fromkeys(c for c, _, _ in blocks) if c not in last]
    at = collections.defaultdict(list)
    for chrom, n in last.items():
        at[n].append(chrom)
    hdr = header(args, len(blocks), n_tx)

    def write_seqids(out, seqids):
        for chrom in seqids:
            for lines in by_seq[chrom]:
                out.write('\n'.join(lines) + '\n')

    with fomo_io.open_text(ref_path) as fh, open(path, 'w', encoding='utf-8') as out:
        for n, line in enumerate(fh):
            if line.startswith('##FASTA'):
                write_seqids(out, rest)
                rest = []
                out.write(line)
                out.writelines(fh)
                break
            if not line.endswith('\n'):
                line += '\n'
            if n == 0 and not line.startswith('##gff-version'):
                out.write('\n'.join(hdr) + '\n')  # reference has no header: add one
            out.write(line)
            if n == 0 and line.startswith('##gff-version'):
                out.write(hdr[1] + '\n')
            write_seqids(out, at.get(n, ()))
        write_seqids(out, rest)


def gffread_ok(exe, path):
    r = subprocess.run([exe, path, '-o', os.devnull], capture_output=True, check=False)
    return r.returncode == 0, r.stderr.decode('utf-8', 'replace').strip()[-2000:]


def gzip_file(src, dst):
    with open(src, 'rb') as i, gzip.open(dst, 'wb', compresslevel=6) as o:
        shutil.copyfileobj(i, o)


# ------------------------------------------------------------------ steps of main
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--gff', required=True, help='fomo GFF3 (.gz ok)')
    p.add_argument(
        '--ref', help='reference annotation of the target (GFF3/GTF, .gz ok); omit for an unannotated target'
    )
    p.add_argument('--prefix', required=True, help='output path prefix, e.g. outdir/<target>.curated.lnc_RNA')
    p.add_argument('--target', help='target species as written in source= (default: from the --gff file name)')
    p.add_argument(
        '--track', help='track of the input, e.g. lnc_RNA, mRNA, lnc_RNA.decoy (default: from the file name)'
    )
    p.add_argument(
        '--id-prefix',
        help='prefix of curated gene / transcript / exon IDs (default: fomo_<tag>, e.g. fomo_lncRNA); '
        'change it only if the IDs clash with the reference',
    )
    p.add_argument('--min-species', type=int, default=2)
    p.add_argument('--strand', choices=['same', 'both'], default='same', help='reference overlap strand (default same)')
    p.add_argument('--keep-self', action='store_true', help="count the target's own projections as support")
    p.add_argument(
        '--divergent-dist',
        type=int,
        default=1000,
        help="max distance (bp) between the 5' ends of a curated gene and an opposite-strand reference "
        'transcript pointing away from it for ref_location=divergent (default 1000; 0 disables)',
    )
    p.add_argument('--gffcompare', default='gffcompare')
    p.add_argument('--gffread', default='gffread')
    p.add_argument('--skip-gffread', action='store_true', help='do not validate outputs with gffread')
    p.add_argument('--identity', help='identity of the curation data record (JSON object); without it no record')
    p.add_argument('--keep-tmp', action='store_true', help='keep the gffcompare work directory (<prefix>.gffcompare/)')
    p.add_argument(
        '--chains-tsv',
        action='store_true',
        help='also write <prefix>.chains.tsv.gz: every supported chain (before the reference filter) with its '
        'species, status, overlapped reference types and gene (audit)',
    )
    return p.parse_args()


def resolve_track(args):
    """-> (target, track, tag, transcript feature type, ID prefix)."""
    base = os.path.basename(args.gff)
    target = args.target or re.sub(r'\.(allModels|from_).*$', '', base)
    match = re.search(r'\.((?:lnc_RNA|mRNA)(?:\.decoy)?)\.', base)
    track = args.track or (match.group(1) if match else 'lnc_RNA')
    tag, ftype = track_info(track)
    return target, track, tag, ftype, esc(args.id_prefix or f'fomo_{tag}')


def filter_models(models: list[Model], target: str, keep_self: bool, observed: dict) -> list[Model]:
    """Drop self projections (unless kept), unstranded and single-exon models; count each step."""
    observed['n_models'] = len(models)
    own = [m for m in models if m.species == target]
    observed['self_excluded'] = 0 if keep_self else len(own)
    if own and not keep_self:
        # The pipeline never aligns a species to itself; such models come from an older run.
        log.warning(f'{len(own):,} self-projected models (source={target}) excluded')
        models = [m for m in models if m.species != target]
    observed['unstranded'] = sum(m.strand not in '+-' for m in models)
    models = [m for m in models if m.strand in '+-']
    observed['single_exon'] = sum(not m.chain for m in models)
    models = [m for m in models if m.chain]
    observed['multi_exon_models'] = len(models)
    observed['species'] = len({m.species for m in models})
    return models


def tag_ref_hits(models: list[Model], ref: str | None, strand_mode: str):
    """Set each model's ref_hits (types of the reference exons it overlaps);
    -> (exon index, transcript-span index, gene_meta, reference is GTF), all None/False without --ref."""
    if not ref:
        return None, None, None, False
    ex_idx, tx_idx, origin, ref_gtf, gene_meta = read_reference(ref)
    log.info(f'reference {ref}: ' + ', '.join(f'{v:,} {k}' for k, v in origin.items()))
    for m in models:
        strands = (m.strand,) if strand_mode == 'same' else ('+', '-')
        m.ref_hits = {t for st in strands for s, e in m.exons for t, _ in ex_idx.hits(m.chr, st, s, e)}
    return ex_idx, tx_idx, gene_meta, ref_gtf


@dataclass(slots=True)
class ChainSupport:
    by_chain: dict
    species_pre: dict  # chain -> species of its models, before the reference filter
    species_post: dict  # ... after it
    supported_pre: set
    supported_post: set


def species_buckets(prefix: str, counts, min_species: int) -> dict:
    """{<prefix>_<n>sp: count == n} for n up to the last bucket, then <prefix>_ge<n+1>sp."""
    top = fomo_stats.CURATION_MAX_SPECIES_BUCKET
    out = {f'{prefix}_{n}sp': sum(c == n for c in counts) for n in range(min_species, top + 1)}
    out[f'{prefix}_ge{top + 1}sp'] = sum(c >= top + 1 for c in counts)
    return out


def chain_support(models: list[Model], args, observed: dict) -> ChainSupport:
    """Exact chains (gffcompare) and their species support before and after the reference filter."""
    out_dir = os.path.dirname(os.path.abspath(args.prefix))
    keep = os.path.abspath(args.prefix) + '.gffcompare'
    tmp = keep if args.keep_tmp else tempfile.mkdtemp(prefix='curate_', dir=out_dir)
    os.makedirs(tmp, exist_ok=True)
    try:
        by_chain, n_files, restored = exact_chains(models, tmp, args.gffcompare)
    finally:
        if not args.keep_tmp:
            shutil.rmtree(tmp, ignore_errors=True)
    observed['gffcompare_inputs'] = n_files
    observed['dup_restored'] = restored
    observed['chains'] = len(by_chain)
    species_pre = {k: {m.species for m in v} for k, v in by_chain.items()}
    species_post = {k: {m.species for m in v if not m.ref_hits} for k, v in by_chain.items()}
    supported_pre = {k for k, s in species_pre.items() if len(s) >= args.min_species}
    supported_post = {k for k, s in species_post.items() if len(s) >= args.min_species}
    if not supported_post <= supported_pre:
        fomo_io.die('a chain is supported after the reference filter but not before it')
    observed['supported_chains'] = len(supported_pre)
    observed.update(species_buckets('supported_chains', [len(species_pre[k]) for k in supported_pre], args.min_species))
    observed['models_in_supported_chains'] = sum(len(by_chain[k]) for k in supported_pre)
    return ChainSupport(by_chain, species_pre, species_post, supported_pre, supported_post)


def discarded_report(models: list[Model], support: ChainSupport, genes_pre) -> dict:
    """What the reference filter removed: models, supported chains and genes, by reference type."""
    by_chain = support.by_chain
    removed = [m for m in models if m.ref_hits]
    lost = support.supported_pre - support.supported_post
    report = {
        'models_removed': len(removed),
        'models_removed_in_supported_chains': sum(1 for k in support.supported_pre for m in by_chain[k] if m.ref_hits),
        'supported_chains_lost': len(lost),
        'supported_chains_kept_after_removal': sum(
            1 for k in support.supported_post if any(m.ref_hits for m in by_chain[k])
        ),
        'genes_lost': sum(1 for g in genes_pre if not any(k in support.supported_post for k in g)),
    }
    models_by_type, models_by_primary, chains_by_type, chains_by_primary = (collections.Counter() for _ in range(4))
    for m in removed:
        models_by_type.update(m.ref_hits)
        models_by_primary[primary_type(m.ref_hits)] += 1
    for k in lost:
        types = set().union(*(m.ref_hits for m in by_chain[k]))
        if not types:
            fomo_io.die(f'supported chain lost without a removed model: {k[0]}:{k[1]}')
        chains_by_type.update(types)
        chains_by_primary[primary_type(types)] += 1
    for name, counts in (
        ('models_by_type', models_by_type),
        ('models_by_primary_type', models_by_primary),
        ('chains_lost_by_type', chains_by_type),
        ('chains_lost_by_primary_type', chains_by_primary),
    ):
        for ftype, n in sorted(counts.items(), key=lambda x: (-x[1], x[0])):
            report[f'{name}:{ftype}'] = n
    return report


def curated_transcripts(chains, support: ChainSupport) -> list[CuratedTranscript]:
    """One transcript per supported chain of a gene, ranked t1, t2, ... (IDs set by the caller)."""
    tx = []
    for k in chains:
        candidates = [m for m in support.by_chain[k] if not m.ref_hits]
        rep = min(candidates, key=lambda m: (-m.exonic_len, m.de, m.id))
        species = sorted(support.species_post[k])
        tx.append(CuratedTranscript(rep, k, species, sorted(candidates, key=lambda m: (m.species, m.id))))
    tx.sort(
        key=lambda t: (
            -len(t.species),
            -len(t.key.introns),
            -t.rep.exonic_len,
            t.rep.exons[0][0],
            t.rep.exons[-1][1],
            t.rep.id,
        )
    )
    return tx


def build_genes(support: ChainSupport, track, id_prefix, reference, args, seq_rank) -> list[CuratedGene]:
    """Curated genes from the chains supported after the filter, with IDs and reference location."""
    ex_idx, tx_idx, gene_meta = reference
    genes, by_hash = [], {}
    for chains in link_genes(support.supported_post):
        tx = curated_transcripts(chains, support)
        chrom, strand = chains[0][0], chains[0][1]
        footprint = merge_exons([x for t in tx for x in t.rep.exons])
        h = gene_hash(track, chrom, strand, [k.introns for k in chains])
        if h in by_hash:
            other = by_hash[h]
            fomo_io.die(
                f'gene hash collision {h}: {chrom}:{strand}:{footprint[0][0]}-{footprint[-1][1]} and '
                f'{other.chr}:{other.strand}:{other.start}-{other.end}'
            )
        gid = f'{id_prefix}_g_{h}'
        location = {'chr': chrom, 'strand': strand, 'exons': footprint}
        loc, ref = ref_location(location, ex_idx, tx_idx, gene_meta, args.divergent_dist)
        check_location(loc, ref, args.divergent_dist, gid)
        for i, t in enumerate(tx, 1):
            t.id, t.exon_prefix = f'{id_prefix}_t{i}_{h}', f'{id_prefix}_t{i}'
        species = sorted(set().union(*(support.species_post[k] for k in chains)))
        gene = CuratedGene(gid, h, chrom, strand, footprint[0][0], footprint[-1][1], tx, chains, species, loc, ref)
        by_hash[h] = gene
        genes.append(gene)
    genes.sort(key=lambda g: (seq_rank[g.chr], g.start, g.strand, g.end, g.hash))
    return genes


def novel_report(genes: list[CuratedGene], min_species: int, has_ref: bool) -> dict:
    """Per gene: species, chains, span, location; per transcript: structure and support."""
    txs = [t for g in genes for t in g.tx]
    report = {
        'genes': len(genes),
        'transcripts': len(txs),
        'transcripts_per_gene_mean': round(len(txs) / len(genes), 3) if genes else 0,
    }
    report.update(species_buckets('genes', [len(g.species) for g in genes], min_species))
    for n in (1, 2):
        report[f'genes_{n}chain'] = sum(g.n_chains == n for g in genes)
    report['genes_ge3chain'] = sum(g.n_chains >= 3 for g in genes)
    reps = [t.rep for t in txs]
    for name, vals in (
        ('n_exons', [len(r.exons) for r in reps]),
        ('transcript_length', [r.exonic_len for r in reps]),
        ('gene_span', [g.end - g.start + 1 for g in genes]),
        ('intron_length', [b - a - 1 for r in reps for a, b in r.chain]),
        ('supporting_models', [len(t.support) for t in txs]),
        ('de', [r.de for r in reps if r.de != math.inf]),
    ):
        for k, v in qstats(vals).items():
            report[f'{name}_{k}'] = v
    locations = collections.Counter(g.location for g in genes)
    for loc in LOCATIONS:
        report[f'genes_location:{loc}'] = locations[loc]
    if has_ref:
        distances = [g.ref['distance'] for g in genes if g.location == 'intergenic' and g.ref]
        for k, v in qstats(distances).items():
            report[f'intergenic_distance_{k}'] = v
    return report


def write_outputs(args, genes: list[CuratedGene], tag: str, ftype: str, ref_gtf: bool) -> None:
    """<prefix>.gff3.gz (validated with gffread) and, with --ref, <prefix>.merged.gff3.gz."""
    blocks = gene_blocks(genes, tag, ftype)
    n_tx = sum(len(g.tx) for g in genes)
    out_plain = args.prefix + '.gff3'
    write_curated(out_plain, args, blocks, n_tx)
    if not args.skip_gffread:
        ok, err = gffread_ok(args.gffread, out_plain)
        if not ok:
            fomo_io.die(f'gffread rejected {out_plain}:\n{err}')
    gzip_file(out_plain, out_plain + '.gz')
    os.remove(out_plain)
    if not args.ref:
        return
    ref_for_merge = args.ref
    if ref_gtf:
        ref_for_merge = args.prefix + '.ref_from_gtf.gff3'
        cmd = [args.gffread, args.ref, '--keep-genes', '-F', '-o', ref_for_merge]
        r = subprocess.run(cmd, capture_output=True, check=False)
        if r.returncode != 0:
            fomo_io.die(f'gffread could not convert the GTF reference:\n{r.stderr.decode("utf-8", "replace")}')
        log.info('reference is GTF: converted with gffread --keep-genes -F for the merged output')
    merged_plain = args.prefix + '.merged.gff3'
    # Not re-validated with gffread: the reference part is copied verbatim (gffread hard-errors on
    # some real references, e.g. trans-spliced strand '?'), and the curated part is validated above.
    write_merged(merged_plain, ref_for_merge, args, blocks, all_ids(genes), n_tx)
    gzip_file(merged_plain, merged_plain + '.gz')
    os.remove(merged_plain)
    if ref_gtf:
        os.remove(ref_for_merge)


def write_chains_tsv(path: str, support: ChainSupport, genes: list[CuratedGene], seq_rank: dict) -> None:
    gene_of = {k: g.gene_id for g in genes for k in g.chains}
    tx_of = {t.key: t for g in genes for t in g.tx}
    columns = 'chr strand introns n_species_pre n_species_post status ref_types gene_id transcript_id representative'
    with gzip.open(path, 'wt', encoding='utf-8') as out:
        out.write('\t'.join(columns.split() + ['species_pre', 'models']) + '\n')
        for k in sorted(support.supported_pre, key=lambda k: (seq_rank.get(k[0], 0), k[2][0], k[1])):
            types = set().union(*(m.ref_hits for m in support.by_chain[k]))
            row = [
                k[0],
                k[1],
                ','.join(f'{a}-{b}' for a, b in k[2]),
                str(len(support.species_pre[k])),
                str(len(support.species_post[k])),
                'kept' if k in support.supported_post else 'lost_ref',
                ','.join(sorted(types)) or '.',
                gene_of.get(k, '.'),
                tx_of[k].id if k in tx_of else '.',
                tx_of[k].rep.id if k in tx_of else '.',
                ','.join(sorted(support.species_pre[k])),
                ','.join(m.id for m in support.by_chain[k]),
            ]
            out.write('\t'.join(row) + '\n')


def write_report(prefix: str, report: dict) -> None:
    with open(prefix + '.curation.tsv', 'w', encoding='utf-8') as out:
        out.write('table\tmetric\tvalue\n')
        for table in REPORT_TABLES:
            for k, v in report[table].items():
                out.write(f'{table}\t{k}\t{v}\n')
    with open(prefix + '.curation.json', 'w', encoding='utf-8') as out:
        json.dump(report, out, indent=1)


def main() -> int:
    args = parse_args()
    t0 = time.time()
    target, track, tag, ftype, id_prefix = resolve_track(args)
    os.makedirs(os.path.dirname(os.path.abspath(args.prefix)), exist_ok=True)

    observed: dict = {}
    models = filter_models(read_models(args.gff), target, args.keep_self, observed)
    log.info(
        f'{observed["n_models"]:,} models in {args.gff}; {observed["multi_exon_models"]:,} multi-exon non-self '
        f'from {observed["species"]} species'
    )
    ex_idx, tx_idx, gene_meta, ref_gtf = tag_ref_hits(models, args.ref, args.strand)
    support = chain_support(models, args, observed)
    genes_pre = link_genes(support.supported_pre)
    observed['genes'] = len(genes_pre)
    discarded = discarded_report(models, support, genes_pre) if args.ref else {}

    seq_rank = {c: i for i, c in enumerate(dict.fromkeys(m.chr for m in models))}
    genes = build_genes(support, track, id_prefix, (ex_idx, tx_idx, gene_meta), args, seq_rank)
    novel = novel_report(genes, args.min_species, bool(args.ref))

    write_outputs(args, genes, tag, ftype, ref_gtf)
    if args.chains_tsv:
        write_chains_tsv(args.prefix + '.chains.tsv.gz', support, genes, seq_rank)
    params = {
        'gff': args.gff,
        'ref': args.ref,
        'target': target,
        'track': track,
        'biotype': tag,
        'id_prefix': id_prefix,
        'min_species': args.min_species,
        'strand': args.strand,
        'keep_self': args.keep_self,
        'divergent_dist': args.divergent_dist,
    }
    report = {'observed': observed, 'discarded_by_reference': discarded, 'novel': novel, 'params': params}
    write_report(args.prefix, report)
    if args.identity:
        data = {table: report[table] for table in REPORT_TABLES} | {'has_reference': bool(args.ref)}
        record = fomo_stats.make_record('curation', fomo_stats.identity_arg(args.identity), data)
        fomo_stats.write_records(args.prefix + '.curation.data.json', record)
    n_tx = sum(len(g.tx) for g in genes)
    log.info(
        f'supported chains {observed["supported_chains"]:,} -> {len(support.supported_post):,} after reference '
        f'filter; {len(genes):,} curated genes, {n_tx:,} transcripts -> {args.prefix}.gff3.gz '
        f'in {time.time() - t0:.0f}s'
    )
    return 0


if __name__ == '__main__':
    fomo_io.run(main)
