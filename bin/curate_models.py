#!/usr/bin/env python3
"""Curate fomo projected models: report intron chains shared exactly by >= N species (plans/21).

Input:  a fomo GFF3 (<target>.allModels.<gtype>.raw.gff3 or any per-source / top3 subset; .gz ok)
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
  5. One representative per gene: most SJs -> the chain shared by the most species (after the
     filter) -> longest exonic length -> lowest de -> ID. It lists its supporting species and the IDs of
     every model with its chain.
  6. Location relative to the reference (ref_location on gene and transcript; first match wins):
     antisense_exonic (exon overlaps an opposite-strand reference exon) > intronic (span inside one
     same-strand reference transcript span) > sense_span_overlap (overlaps a same-strand reference
     transcript span, not contained) > antisense_intronic (overlaps an opposite-strand reference
     transcript span, no exon overlap) > divergent (no overlap; head-to-head with an opposite-strand
     reference transcript, 5' ends <= --divergent-dist apart) > intergenic; no_reference without
     --ref. A reference transcript span is its first-to-last exon (CDS/UTR or own span if exon-less);
     unstranded reference rows count on both strands. The transcript also gets ref_nearest (type of
     the nearest reference transcript, either strand) and ref_distance (bp, 0 when overlapping).

Output (--prefix P)
  P.gff3.gz           curated genes (gene + representative transcript + exons)
  P.merged.gff3.gz    the reference verbatim + curated genes, each after its seqid's last
                      reference line (with --ref only)
  P.curation.tsv      report: table / metric / value (observed, discarded_by_reference, novel)
  P.curation.json     the same numbers
"""
import argparse
import bisect
import collections
from array import array
import gzip
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

GENE_LEVEL = {'gene', 'pseudogene', 'ncRNA_gene', 'region', 'chromosome', 'scaffold', 'contig'}
CDS_LIKE = {'CDS', 'five_prime_UTR', 'three_prime_UTR', 'UTR', 'start_codon', 'stop_codon'}
NOT_TRANSCRIPT = {'exon', 'intron', 'polyA_site', 'polyA_signal_sequence', 'TSS', 'transcription_start_site',
                  'sequence_feature', 'Selenocysteine', 'stop_codon_redefined_as_selenocysteine'} | CDS_LIKE
GENERIC_TYPES = {'transcript', 'ncRNA', 'primary_transcript', 'RNA', 'gene', 'pseudogene', 'ncRNA_gene'}
BIOTYPE_KEYS = ('transcript_biotype', 'transcript_type', 'biotype', 'ncrna_class', 'gene_biotype', 'gene_type', 'gbkey')
ESCAPE = re.compile(r'[;=&,\t\n\r]')


def log(msg):
    print(f'[curate_models] {msg}', file=sys.stderr, flush=True)


def die(msg):
    log(f'ERROR: {msg}')
    sys.exit(1)


def opener(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def esc(v):
    return ESCAPE.sub(lambda c: '%{:02X}'.format(ord(c.group(0))), str(v))


def gff3_attrs(s):
    """-> list of (key, raw value), order kept."""
    out = []
    for kv in s.strip().rstrip(';').split(';'):
        if '=' in kv:
            k, v = kv.split('=', 1)
            out.append((k.strip(), v.strip()))
    return out


def gtf_attrs(s):
    d = {}
    for k, v in re.findall(r'(\S+)\s+"([^"]*)"', s):
        d.setdefault(k, v)
    return d


def is_gtf(path):
    with opener(path) as fh:
        for line in fh:
            if line[0] != '#':
                f = line.rstrip('\n').split('\t')
                if len(f) >= 9:
                    return '"' in f[8] and '=' not in f[8].split('"')[0]
    return False


def track_info(track):
    t = track or ''
    decoy = 'decoy' in t
    kind = 'mRNA' if re.search(r'(^|[._])mRNA($|[._])', t) else 'lncRNA'
    return 'fomo_' + ('decoy_' if decoy else '') + kind


# ------------------------------------------------------------------ models
class Model:
    """One projected model. Compact (slots, interned seqid/species, raw transcript line kept as one
    string); m['f'] / m['a'] parse that line on demand, only needed for representatives."""
    __slots__ = ('id', 'chr', 'strand', 'exons', 'species', 'de', 'line', 'chain', 'ref_hits')

    def __getitem__(self, k):
        if k == 'f':
            return self.line.rstrip('\n').split('\t')
        if k == 'a':
            return gff3_attrs(self['f'][8])
        return getattr(self, k)

    def __setitem__(self, k, v):
        setattr(self, k, v)


DE_RE = re.compile(r'(?:^|;)de=([^;]*)')
SOURCE_RE = re.compile(r'(?:^|;)source=([^;]*)')
ID_RE = re.compile(r'(?:^|;)ID=([^;]*)')
PARENT_RE = re.compile(r'(?:^|;)Parent=([^;]*)')


def read_models(path):
    """fomo GFF3 -> (list of Model, seqid order). Exits 1 on structural problems."""
    tx, exons, order = {}, {}, {}
    intern = sys.intern
    with opener(path) as fh:
        for n, line in enumerate(fh, 1):
            if line[0] == '#' or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                die(f'{path}:{n}: not a 9-column GFF3 line')
            chrom = intern(f[0])
            order.setdefault(chrom, len(order))
            if f[2] == 'exon':
                m = PARENT_RE.search(f[8])
                if not m:
                    die(f'{path}:{n}: exon without Parent')
                for p in m.group(1).split(','):
                    e = exons.get(p)
                    if e is None:
                        exons[p] = [chrom, f[6], [(int(f[3]), int(f[4]))]]
                    elif e[0] is not chrom or e[1] != f[6]:
                        die(f'transcript {p}: exons on several seqids/strands')
                    else:
                        e[2].append((int(f[3]), int(f[4])))
            elif f[2] not in GENE_LEVEL and f[2] not in NOT_TRANSCRIPT:
                m = ID_RE.search(f[8])
                if not m:
                    continue
                tid = m.group(1).strip()
                if tid in tx:
                    die(f'{path}:{n}: duplicate transcript ID {tid}')
                tx[tid] = line
    missing = set(exons) - set(tx)
    if missing:
        die(f'{len(missing)} exon parent(s) not defined as transcripts, e.g. {sorted(missing)[0]}')
    models = []
    for tid, line in tx.items():
        e = exons.pop(tid, None)
        if not e:
            continue                                              # transcript-like rows without exons
        f = line.split('\t', 9)
        chrom, strand, ex = e
        if (chrom, strand) != (f[0], f[6]):
            die(f'transcript {tid}: exons on {chrom}{strand}, transcript row on {f[0]}{f[6]}')
        ex = tuple(sorted(ex))
        for (s1, e1), (s2, e2) in zip(ex, ex[1:]):
            if s2 <= e1 + 1:
                die(f'transcript {tid}: overlapping or abutting exons {s1}-{e1} / {s2}-{e2}')
        ms = SOURCE_RE.search(f[8])
        sp = ms.group(1).strip() if ms else (tid.rsplit('|', 1)[1] if '|' in tid else '')
        if not sp:
            die(f'transcript {tid}: no source= attribute and no |<source> ID suffix')
        md = DE_RE.search(f[8])
        try:
            de = float(md.group(1)) if md else float('nan')
        except ValueError:
            de = float('nan')
        m = Model()
        m.id, m.chr, m.strand, m.exons, m.species, m.line = tid, chrom, strand, ex, intern(sp), line
        m.de = de if not math.isnan(de) else float('inf')
        m.chain = tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))
        m.ref_hits = set()
        models.append(m)
    return models, list(order)


# ------------------------------------------------------------------ reference
class RefIndex:
    """(chr, strand) -> intervals sorted by start, in typed arrays, each with a type label;
    hits() returns every label overlapping [a, b] (>= 1 bp)."""

    def __init__(self):
        self.labels, self._lab = [], {}
        self.build = collections.defaultdict(lambda: (array('l'), array('l'), array('l')))
        self.idx = {}

    def add(self, chrom, strand, s, e, label):
        li = self._lab.get(label)
        if li is None:
            li = self._lab[label] = len(self.labels)
            self.labels.append(label)
        st, en, lb = self.build[(chrom, strand)]
        st.append(s)
        en.append(e)
        lb.append(li)

    def freeze(self):
        for key, (st, en, lb) in self.build.items():
            order = sorted(range(len(st)), key=st.__getitem__)
            st2 = array('l', (st[k] for k in order))
            en2 = array('l', (en[k] for k in order))
            lb2 = array('l', (lb[k] for k in order))
            self.idx[key] = (st2, en2, lb2, max(e - s for s, e in zip(st2, en2)) + 1)
        self.build = None
        return self

    def hits(self, chrom, strand, a, b):
        x = self.idx.get((chrom, strand))
        out = set()
        if x:
            st, en, lb, ml = x
            for i in range(bisect.bisect_left(st, a - ml), bisect.bisect_right(st, b)):
                if en[i] >= a:
                    out.add(self.labels[lb[i]])
        return out

    def intervals(self, chrom, strand, a, b):
        """-> [(start, end, label)] overlapping [a, b]."""
        x = self.idx.get((chrom, strand))
        out = []
        if x:
            st, en, lb, ml = x
            for i in range(bisect.bisect_left(st, a - ml), bisect.bisect_right(st, b)):
                if en[i] >= a:
                    out.append((st[i], en[i], self.labels[lb[i]]))
        return out

    def nearest(self, chrom, strand, a, b):
        """-> (distance, label) of the closest interval to [a, b] (0 if overlapping), or None."""
        x = self.idx.get((chrom, strand))
        if not x:
            return None
        st, en, lb, ml = x
        best = None
        ov = self.intervals(chrom, strand, a, b)
        if ov:
            return (0, min(l for _, _, l in ov))
        j = bisect.bisect_right(st, b)                          # first interval starting right of b
        if j < len(st):
            best = (st[j] - b, self.labels[lb[j]])
        i = bisect.bisect_left(st, a) - 1                       # intervals starting left of a: max end
        max_end, lab = None, None
        while i >= 0:
            if max_end is not None and st[i] + ml <= max_end:
                break
            if max_end is None or en[i] > max_end:
                max_end, lab = en[i], self.labels[lb[i]]
            i -= 1
        if max_end is not None and (best is None or a - max_end < best[0]):
            best = (a - max_end, lab)
        return best


def read_reference(path):
    """-> (exon RefIndex, transcript-span RefIndex, Counter of interval origins, is GTF).

    Every exon row counts, labelled by its parent's type (generic types refined by biotype). A
    transcript-like feature without exons counts by its CDS/UTR rows or, with no children at all,
    by its own span. Unstranded rows ('.', '?') count on both strands."""
    gtf = is_gtf(path)
    intern = sys.intern
    feat = {}                                  # transcript/gene id -> type label
    pending = []                               # exons / CDS seen before their parent: (pid, kind, chrom, strands, s, e)
    cds = collections.defaultdict(list)        # parent -> CDS/UTR intervals (used only if it has no exon)
    with_exons, has_child, parented = set(), set(), []
    spans = {}                                 # (parent, chrom, strand) -> (min, max)
    ex_idx = RefIndex()
    origin = collections.Counter()

    def lab_of(t, a):
        if t in GENERIC_TYPES:
            bt = next((a[k] for k in BIOTYPE_KEYS if a.get(k)), None)
            if bt and bt != t:
                return intern(f'{t}:{bt}')
        return intern(t)

    def add_iv(pid, chrom, strands, s, e):
        lab = feat.get(pid, 'unknown_parent')
        for st in strands:
            ex_idx.add(chrom, st, s, e, lab)
            k = (pid, chrom, st)
            lo_hi = spans.get(k)
            spans[k] = (s, e) if lo_hi is None else (min(lo_hi[0], s), max(lo_hi[1], e))

    with opener(path) as fh:
        for line in fh:
            if line[0] == '#':
                if line.startswith('##FASTA'):
                    break
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            chrom = intern(f[0])
            strands = (f[6],) if f[6] in ('+', '-') else ('+', '-')      # '.' / '?' -> both
            s, e = int(f[3]), int(f[4])
            if gtf:
                a = gtf_attrs(f[8])
                tid = a.get('transcript_id')
                if not tid:
                    continue
                if f[2] == 'transcript' or (f[2] == 'exon' and tid not in feat):
                    feat[tid] = lab_of('transcript', a)
                if f[2] == 'exon':
                    with_exons.add(tid)
                    add_iv(tid, chrom, strands, s, e)
                    origin['exon'] += 1
                elif f[2] in CDS_LIKE:
                    cds[tid].append((chrom, strands, s, e))
                continue
            t = f[2]
            if t == 'exon':
                m = PARENT_RE.search(f[8])
                parents = m.group(1).split(',') if m else ['<orphan exon>']
                origin['exon'] += 1
                for p in parents:
                    with_exons.add(p)
                    if p in feat:
                        add_iv(p, chrom, strands, s, e)
                    else:
                        pending.append((p, chrom, strands, s, e))
                continue
            if t in CDS_LIKE:
                m = PARENT_RE.search(f[8])
                for p in (m.group(1).split(',') if m else ()):
                    cds[p].append((chrom, strands, s, e))
                continue
            a = dict(gff3_attrs(f[8]))
            parents = [p for p in a.get('Parent', '').split(',') if p]
            if 'ID' in a:
                feat[a['ID']] = lab_of(t, a)
            has_child.update(parents)
            if parents and 'ID' in a and t not in NOT_TRANSCRIPT and t not in GENE_LEVEL:
                parented.append((a['ID'], chrom, strands, s, e))
    for p, chrom, strands, s, e in pending:                     # exons listed before their parent
        add_iv(p, chrom, strands, s, e)
    for p, rows in cds.items():                                 # CDS/UTR of a model without exons
        if p in with_exons:
            continue
        for chrom, strands, s, e in rows:
            add_iv(p, chrom, strands, s, e)
            origin['cds_utr_without_exon'] += 1
    for fid, chrom, strands, s, e in parented:                  # childless transcript-like features
        if fid in has_child or fid in with_exons or fid in cds:
            continue
        add_iv(fid, chrom, strands, s, e)
        origin['childless_feature_span'] += 1
    tx_idx = RefIndex()
    for (p, chrom, st), (lo, hi) in spans.items():
        tx_idx.add(chrom, st, lo, hi, feat.get(p, 'unknown_parent'))
    return ex_idx.freeze(), tx_idx.freeze(), origin, gtf


# ------------------------------------------------------------------ gffcompare
def run_gffcompare(models, tmp, exe):
    """-> list of cliques (lists of model indices); writes one GTF per species with integer IDs."""
    by_sp = collections.defaultdict(list)
    for i, m in enumerate(models):
        by_sp[m['species']].append(i)
    files = []
    for k, sp in enumerate(sorted(by_sp)):
        path = os.path.join(tmp, f'q{k}.gtf')
        with open(path, 'w') as o:
            for i in by_sp[sp]:
                m = models[i]
                tag = f'gene_id "g{i}"; transcript_id "t{i}";'
                o.write(f"{m['chr']}\tfomo\ttranscript\t{m['exons'][0][0]}\t{m['exons'][-1][1]}\t.\t{m['strand']}\t.\t{tag}\n")
                for s, e in m['exons']:
                    o.write(f"{m['chr']}\tfomo\texon\t{s}\t{e}\t.\t{m['strand']}\t.\t{tag}\n")
        files.append(path)
    with open(os.path.join(tmp, 'inputs.txt'), 'w') as o:
        o.write('\n'.join(files) + '\n')
    r = subprocess.run([exe, '--no-merge', '-o', 'cmp', '-i', 'inputs.txt'], cwd=tmp, capture_output=True)
    if r.returncode != 0:
        die(f"gffcompare failed (exit {r.returncode}):\n{r.stderr.decode('utf-8', 'replace')}")
    cliques = []
    with open(os.path.join(tmp, 'cmp.tracking')) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            ids = []
            for q in f[4:]:
                if q == '-':
                    continue
                for item in q.split(','):                             # qN:g<i>|t<i>|exons|fpkm|tpm|cov|len
                    ids.append(int(item.split(':', 1)[1].split('|')[1][1:]))
            cliques.append(ids)
    return cliques, len(files)


def exact_chains(models, tmp, exe):
    """-> {chain key: [models]} as decided by gffcompare, verified against the chains derived here."""
    if not models:
        return {}, 0, 0
    cliques, n_files = run_gffcompare(models, tmp, exe)
    key = lambda m: (m['chr'], m['strand'], m['chain'])
    by_key = collections.defaultdict(list)
    for m in models:
        by_key[key(m)].append(m)
    seen = set()
    for ids in cliques:
        ks = {key(models[i]) for i in ids}
        if len(ks) != 1:
            die(f'gffcompare clique with {len(ks)} different intron chains (ids {ids[:5]}); check gffcompare version/options')
        k = ks.pop()
        if k in seen:
            die(f'intron chain {k[0]}:{k[1]} reported in two gffcompare cliques')
        seen.add(k)
    if seen != set(by_key):
        die(f'{len(set(by_key) - seen)} intron chain(s) missing from gffcompare tracking')
    dup = sum(len(by_key[k]) for k in seen) - sum(len(c) for c in cliques)
    return by_key, n_files, dup


# ------------------------------------------------------------------ genes / representative
def link_genes(chains):
    """chains: iterable of chain keys -> list of lists of chain keys sharing >= 1 intron (connected)."""
    chains = sorted(chains)
    par = list(range(len(chains)))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    first = {}
    for i, (c, s, ch) in enumerate(chains):
        for j in ch:
            k = (c, s, j)
            if k in first:
                par[find(i)] = find(first[k])
            else:
                first[k] = i
    g = collections.defaultdict(list)
    for i, k in enumerate(chains):
        g[find(i)].append(k)
    return list(g.values())


def exonic_len(m):
    return sum(e - s + 1 for s, e in m['exons'])


def qstats(vals):
    if not vals:
        return dict(n=0)
    v = sorted(vals)
    q = lambda p: v[min(len(v) - 1, int(round(p * (len(v) - 1))))]
    return dict(n=len(v), min=v[0], p25=q(.25), median=statistics.median(v), mean=round(statistics.fmean(v), 2),
                p75=q(.75), max=v[-1])


def primary_type(types):
    def rank(t):
        tl = t.lower()
        if t == 'mRNA' or 'protein_coding' in tl:
            return 0
        if t in ('lnc_RNA', 'lncRNA') or 'lncrna' in tl or 'lincrna' in tl:
            return 1
        if 'antisense' in tl:
            return 2
        if 'pseudo' in tl:
            return 5
        if t.startswith('transcript'):
            return 3
        if any(x in tl for x in ('trna', 'rrna', 'snorna', 'snrna', 'mirna', 'scrna', 'y_rna', 'srp', 'rnase', 'guide', 'ncrna', 'misc')):
            return 4
        return 6
    return min(types, key=lambda t: (rank(t), t))


LOCATIONS = ('antisense_exonic', 'intronic', 'sense_span_overlap', 'antisense_intronic', 'divergent',
             'intergenic', 'no_reference')


def ref_location(r, ex_idx, tx_idx, divergent_dist):
    """-> (class, nearest type label or '.', distance or '.') of representative r vs the reference."""
    if ex_idx is None:
        return 'no_reference', None, None
    c, strand, ex = r['chr'], r['strand'], r['exons']
    S, E = ex[0][0], ex[-1][1]
    opp = '-' if strand == '+' else '+'
    near = [n for n in (tx_idx.nearest(c, strand, S, E), tx_idx.nearest(c, opp, S, E)) if n]
    nb = min(near) if near else None
    nlab, ndist = (nb[1], nb[0]) if nb else ('.', '.')
    if any(ex_idx.hits(c, opp, s, e) for s, e in ex):
        return 'antisense_exonic', nlab, ndist
    same = tx_idx.intervals(c, strand, S, E)
    if any(s <= S and e >= E for s, e, _ in same):
        return 'intronic', nlab, ndist
    if same:
        return 'sense_span_overlap', nlab, ndist
    if tx_idx.intervals(c, opp, S, E):
        return 'antisense_intronic', nlab, ndist
    # head-to-head: an opposite-strand transcript whose 5' end lies within divergent_dist upstream of ours
    win = (S - divergent_dist, S - 1) if strand == '+' else (E + 1, E + divergent_dist)
    if divergent_dist > 0 and tx_idx.intervals(c, opp, *win):
        return 'divergent', nlab, ndist
    return 'intergenic', nlab, ndist


# ------------------------------------------------------------------ output
def gene_blocks(genes, gene_prefix, biotype, track):
    """-> list of (seqid, start, [lines]) per curated gene, numbered in output order."""
    blocks = []
    for i, g in enumerate(genes, 1):
        r = g['rep']
        gid = f'{gene_prefix}_G{i:06d}'
        g['gene_id'] = gid
        s, e = r['exons'][0][0], r['exons'][-1][1]
        lines = ['\t'.join([r['chr'], 'fomo', 'gene', str(s), str(e), '.', r['strand'], '.',
                            f"ID={esc(gid)};gene_biotype={biotype};n_chains={g['n_chains']};"
                            f"n_species_gene={g['n_species_gene']};ref_location={g['location']};curated_track={esc(track)}"])]
        extra = [('transcript_biotype', biotype), ('n_supporting_species', str(len(g['species']))),
                 ('supporting_species', ','.join(esc(x) for x in g['species'])),
                 ('supporting_models', ','.join(esc(m['id']) for m in g['support'])),
                 ('n_sj', str(len(r['chain']))), ('ref_location', g['location'])]
        if g['nearest'] is not None:
            extra += [('ref_nearest', esc(g['nearest'])), ('ref_distance', str(g['distance']))]
        keys = {k for k, _ in extra}
        rid = r['id']
        a = [('ID', rid), ('Parent', esc(gid))]
        a += [(k, v) for k, v in r['a'] if k not in ('ID', 'Parent') and k not in keys] + extra
        lines.append('\t'.join([r['chr'], 'fomo', r['f'][2], str(s), str(e), r['f'][5], r['strand'], '.',
                                ';'.join(f'{k}={v}' for k, v in a)]))
        for j, (xs, xe) in enumerate(r['exons'], 1):
            lines.append('\t'.join([r['chr'], 'fomo', 'exon', str(xs), str(xe), '.', r['strand'], '.',
                                    f'ID={rid}.exon{j};Parent={rid}']))
        blocks.append((r['chr'], s, lines))
    return blocks


def header(args, n_genes):
    return ['##gff-version 3',
            f'# curate_models: intron chains shared exactly by >={args.min_species} species (gffcompare --no-merge), '
            f'linked into genes by shared splice junctions; one representative per gene; reference filter: '
            f"{'same-strand' if args.strand == 'same' else 'either-strand'} overlap with any reference exon; "
            f'{n_genes} genes']


def write_curated(path, args, blocks):
    with open(path, 'w') as o:
        o.write('\n'.join(header(args, len(blocks))) + '\n')
        for _, _, lines in blocks:
            o.write('\n'.join(lines) + '\n')


def write_merged(path, ref_path, args, blocks, new_ids):
    """Reference verbatim; each seqid's curated blocks right after its last reference feature line.
    Exits 1 if a curated ID is already used by the reference."""
    last = {}
    id_re = re.compile(r'(?:^|;)\s*ID=([^;\t\n]+)')
    with opener(ref_path) as fh:
        for n, line in enumerate(fh):
            if line.startswith('##FASTA'):
                break
            if line[0] != '#' and line.strip():
                f = line.split('\t', 9)
                last[f[0]] = n
                if len(f) >= 9:
                    m = id_re.search(f[8])
                    if m and m.group(1).strip() in new_ids:
                        die(f'curated ID {m.group(1).strip()} already used in the reference; set --gene-prefix')
    by_seq = collections.defaultdict(list)
    for c, _, lines in blocks:
        by_seq[c].append(lines)
    rest = [c for c in dict.fromkeys(c for c, _, _ in blocks) if c not in last]
    at = collections.defaultdict(list)
    for c, n in last.items():
        at[n].append(c)
    with opener(ref_path) as fh, open(path, 'w') as o:
        wrote_hdr = False
        for n, line in enumerate(fh):
            if line.startswith('##FASTA'):
                for c in rest:
                    for lines in by_seq[c]:
                        o.write('\n'.join(lines) + '\n')
                rest = []
                o.write(line)
                o.writelines(fh)
                break
            if not line.endswith('\n'):
                line += '\n'
            if n == 0 and not line.startswith('##gff-version'):
                o.write('\n'.join(header(args, len(blocks))) + '\n')   # reference has no header: add one
                wrote_hdr = True
            o.write(line)
            if not wrote_hdr and n == 0:
                o.write(header(args, len(blocks))[1] + '\n')
                wrote_hdr = True
            for c in at.get(n, ()):
                for lines in by_seq[c]:
                    o.write('\n'.join(lines) + '\n')
        for c in rest:
            for lines in by_seq[c]:
                o.write('\n'.join(lines) + '\n')


def gffread_ok(exe, path):
    r = subprocess.run([exe, path, '-o', os.devnull], capture_output=True)
    return r.returncode == 0, r.stderr.decode('utf-8', 'replace').strip()[-2000:]


def gzip_file(src, dst):
    with open(src, 'rb') as i, gzip.open(dst, 'wb', compresslevel=6) as o:
        shutil.copyfileobj(i, o)


# ------------------------------------------------------------------ main
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--gff', required=True, help='fomo GFF3 (.gz ok)')
    p.add_argument('--ref', help='reference annotation of the target (GFF3/GTF, .gz ok); omit for an unannotated target')
    p.add_argument('--prefix', required=True, help='output path prefix, e.g. outdir/<target>.curated.lnc_RNA')
    p.add_argument('--target', help='target species as written in source= (default: from the --gff file name)')
    p.add_argument('--track', help='gene type of the input, e.g. lnc_RNA, mRNA, lnc_RNA.decoy (default: from the file name)')
    p.add_argument('--gene-prefix', help='prefix of curated gene IDs (default: basename of --prefix)')
    p.add_argument('--min-species', type=int, default=2)
    p.add_argument('--strand', choices=['same', 'both'], default='same', help='reference overlap strand (default same)')
    p.add_argument('--keep-self', action='store_true', help="count the target's own projections as support")
    p.add_argument('--divergent-dist', type=int, default=1000,
                   help="max distance (bp) between the 5' ends of a curated gene and an opposite-strand reference "
                        "transcript pointing away from it for ref_location=divergent (default 1000; 0 disables)")
    p.add_argument('--gffcompare', default='gffcompare')
    p.add_argument('--gffread', default='gffread')
    p.add_argument('--skip-gffread', action='store_true', help='do not validate outputs with gffread')
    p.add_argument('--keep-tmp', action='store_true', help='keep the gffcompare work directory (<prefix>.gffcompare/)')
    p.add_argument('--chains-tsv', action='store_true', help='also write <prefix>.chains.tsv.gz: every supported chain (before the '
                   'reference filter) with its species, status, overlapped reference types and gene (audit)')
    args = p.parse_args()
    t0 = time.time()

    base = os.path.basename(args.gff)
    target = args.target or re.sub(r'\.(allModels|top\d+|from_).*$', '', base)
    track = args.track or next((m.group(1) for m in [re.search(r'\.((?:lnc_RNA|mRNA)(?:\.decoy)?)\.', base)] if m), 'lnc_RNA')
    biotype = track_info(track)
    gene_prefix = args.gene_prefix or os.path.basename(args.prefix)
    out_dir = os.path.dirname(os.path.abspath(args.prefix))
    os.makedirs(out_dir, exist_ok=True)

    models, _ = read_models(args.gff)
    R = collections.OrderedDict()
    obs, dis, nov = collections.OrderedDict(), collections.OrderedDict(), collections.OrderedDict()
    obs['n_models'] = len(models)
    is_self = [m for m in models if m['species'] == target]
    obs['self_excluded'] = 0 if args.keep_self else len(is_self)
    if not is_self and not args.keep_self:
        log(f'note: no models with source={target} (pure target, or set --target)')
    if not args.keep_self:
        models = [m for m in models if m['species'] != target]
    obs['unstranded'] = sum(m['strand'] not in '+-' for m in models)
    models = [m for m in models if m['strand'] in '+-']
    obs['single_exon'] = sum(not m['chain'] for m in models)
    models = [m for m in models if m['chain']]
    obs['multi_exon_models'] = len(models)
    obs['species'] = len({m['species'] for m in models})
    log(f"{obs['n_models']:,} models in {args.gff}; {obs['multi_exon_models']:,} multi-exon non-self from {obs['species']} species")

    # reference filter (per model; the overlapped types are kept for the report)
    ref_gtf = False
    ex_idx = tx_idx = None
    if args.ref:
        ex_idx, tx_idx, origin, ref_gtf = read_reference(args.ref)
        log(f'reference {args.ref}: ' + ', '.join(f'{v:,} {k}' for k, v in origin.items()))
        for m in models:
            hits = set()
            strands = (m['strand'],) if args.strand == 'same' else ('+', '-')
            for st in strands:
                for s, e in m['exons']:
                    hits |= ex_idx.hits(m['chr'], st, s, e)
            m['ref_hits'] = hits
    else:
        for m in models:
            m['ref_hits'] = set()

    # exact chains (gffcompare), support before / after the filter
    tmp = (os.path.abspath(args.prefix) + '.gffcompare') if args.keep_tmp else tempfile.mkdtemp(prefix='curate_', dir=out_dir)
    os.makedirs(tmp, exist_ok=True)
    by_chain, n_files, dup = exact_chains(models, tmp, args.gffcompare)
    if not args.keep_tmp:
        shutil.rmtree(tmp, ignore_errors=True)
    obs['gffcompare_inputs'] = n_files
    obs['dup_restored'] = dup
    obs['chains'] = len(by_chain)
    sp_pre = {k: {m['species'] for m in v} for k, v in by_chain.items()}
    sp_post = {k: {m['species'] for m in v if not m['ref_hits']} for k, v in by_chain.items()}
    sup_pre = {k for k, s in sp_pre.items() if len(s) >= args.min_species}
    sup_post = {k for k, s in sp_post.items() if len(s) >= args.min_species}
    assert sup_post <= sup_pre
    obs['supported_chains'] = len(sup_pre)
    for n in range(args.min_species, 5):
        obs[f'supported_chains_{n}sp'] = sum(len(sp_pre[k]) == n for k in sup_pre)
    obs['supported_chains_ge5sp'] = sum(len(sp_pre[k]) >= 5 for k in sup_pre)
    obs['models_in_supported_chains'] = sum(len(by_chain[k]) for k in sup_pre)
    genes_pre = link_genes(sup_pre)
    obs['genes'] = len(genes_pre)

    # discarded by the reference
    if args.ref:
        removed = [m for m in models if m['ref_hits']]
        dis['models_removed'] = len(removed)
        dis['models_removed_in_supported_chains'] = sum(1 for k in sup_pre for m in by_chain[k] if m['ref_hits'])
        lost = sup_pre - sup_post
        dis['supported_chains_lost'] = len(lost)
        dis['supported_chains_kept_after_removal'] = sum(1 for k in sup_post if any(m['ref_hits'] for m in by_chain[k]))
        dis['genes_lost'] = sum(1 for g in genes_pre if not any(k in sup_post for k in g))
        mt, mp, ct, cp = (collections.Counter() for _ in range(4))
        for m in removed:
            mt.update(m['ref_hits'])
            mp[primary_type(m['ref_hits'])] += 1
        for k in lost:
            types = set().union(*(m['ref_hits'] for m in by_chain[k]))
            if not types:
                die(f'supported chain lost without a removed model: {k[0]}:{k[1]}')
            ct.update(types)
            cp[primary_type(types)] += 1
        for name, cnt in (('models_by_type', mt), ('models_by_primary_type', mp),
                          ('chains_lost_by_type', ct), ('chains_lost_by_primary_type', cp)):
            for t, v in sorted(cnt.items(), key=lambda x: (-x[1], x[0])):
                dis[f'{name}:{t}'] = v

    # genes, representative
    genes = []
    for g in link_genes(sup_post):
        cand = [m for k in g for m in by_chain[k] if not m['ref_hits']]
        rep = min(cand, key=lambda m: (-len(m['chain']), -len(sp_post[(m['chr'], m['strand'], m['chain'])]),
                                       -exonic_len(m), m['de'], m['id']))
        rk = (rep['chr'], rep['strand'], rep['chain'])
        support = sorted((m for m in by_chain[rk] if not m['ref_hits']), key=lambda m: (m['species'], m['id']))
        loc, nlab, ndist = ref_location(rep, ex_idx, tx_idx, args.divergent_dist)
        genes.append(dict(rep=rep, chains=g, n_chains=len(g), species=sorted(sp_post[rk]), support=support,
                          n_species_gene=len(set().union(*(sp_post[k] for k in g))),
                          location=loc, nearest=nlab, distance=ndist))
    seq_rank = {c: i for i, c in enumerate(dict.fromkeys(m['chr'] for m in models))}
    genes.sort(key=lambda g: (seq_rank[g['rep']['chr']], g['rep']['exons'][0][0], g['rep']['strand'],
                              g['rep']['exons'][-1][1], g['rep']['id']))
    blocks = gene_blocks(genes, gene_prefix, biotype, track)

    # novel model stats
    nov['genes'] = len(genes)
    for n in range(args.min_species, 5):
        nov[f'genes_{n}sp'] = sum(len(g['species']) == n for g in genes)
    nov['genes_ge5sp'] = sum(len(g['species']) >= 5 for g in genes)
    for n in (1, 2):
        nov[f'genes_{n}chain'] = sum(g['n_chains'] == n for g in genes)
    nov['genes_ge3chain'] = sum(g['n_chains'] >= 3 for g in genes)
    reps = [g['rep'] for g in genes]
    for name, vals in (('n_exons', [len(r['exons']) for r in reps]),
                       ('transcript_length', [exonic_len(r) for r in reps]),
                       ('gene_span', [r['exons'][-1][1] - r['exons'][0][0] + 1 for r in reps]),
                       ('intron_length', [b - a - 1 for r in reps for a, b in r['chain']]),
                       ('supporting_models', [len(g['support']) for g in genes]),
                       ('de', [r['de'] for r in reps if r['de'] != float('inf')])):
        for k, v in qstats(vals).items():
            nov[f'{name}_{k}'] = v
    loc_n = collections.Counter(g['location'] for g in genes)
    for c in LOCATIONS:
        nov[f'genes_location:{c}'] = loc_n[c]

    # write + validate
    out_plain = args.prefix + '.gff3'
    write_curated(out_plain, args, blocks)
    if not args.skip_gffread:
        ok, err = gffread_ok(args.gffread, out_plain)
        if not ok:
            die(f'gffread rejected {out_plain}:\n{err}')
    gzip_file(out_plain, out_plain + '.gz')
    os.remove(out_plain)
    if args.ref:
        new_ids = {g['gene_id'] for g in genes} | {g['rep']['id'] for g in genes} | \
                  {f"{g['rep']['id']}.exon{j}" for g in genes for j in range(1, len(g['rep']['exons']) + 1)}
        ref_for_merge = args.ref
        if ref_gtf:
            ref_for_merge = args.prefix + '.ref_from_gtf.gff3'
            r = subprocess.run([args.gffread, args.ref, '--keep-genes', '-F', '-o', ref_for_merge], capture_output=True)
            if r.returncode != 0:
                die(f"gffread could not convert the GTF reference:\n{r.stderr.decode('utf-8', 'replace')}")
            log('reference is GTF: converted with gffread --keep-genes -F for the merged output')
        merged_plain = args.prefix + '.merged.gff3'
        # not re-validated with gffread: the reference part is copied verbatim (gffread hard-errors on some
        # real references, e.g. RefSeq trans-spliced strand '?'), and the curated part was validated above
        write_merged(merged_plain, ref_for_merge, args, blocks, new_ids)
        gzip_file(merged_plain, merged_plain + '.gz')
        os.remove(merged_plain)
        if ref_gtf:
            os.remove(ref_for_merge)

    if args.chains_tsv:
        gene_of = {k: g['gene_id'] for g in genes for k in g['chains']}
        rep_of = {(g['rep']['chr'], g['rep']['strand'], g['rep']['chain']): g['rep']['id'] for g in genes}
        with gzip.open(args.prefix + '.chains.tsv.gz', 'wt') as o:
            o.write('chr\tstrand\tintrons\tn_species_pre\tn_species_post\tstatus\tref_types\tgene_id\trepresentative\tspecies_pre\tmodels\n')
            for k in sorted(sup_pre, key=lambda k: (seq_rank.get(k[0], 0), k[2][0], k[1])):
                types = set().union(*(m['ref_hits'] for m in by_chain[k]))
                o.write('\t'.join([k[0], k[1], ','.join(f'{a}-{b}' for a, b in k[2]), str(len(sp_pre[k])), str(len(sp_post[k])),
                                   'kept' if k in sup_post else 'lost_ref', ','.join(sorted(types)) or '.',
                                   gene_of.get(k, '.'), rep_of.get(k, '.'), ','.join(sorted(sp_pre[k])),
                                   ','.join(m['id'] for m in by_chain[k])]) + '\n')

    R['observed'], R['discarded_by_reference'], R['novel'] = obs, dis, nov
    R['params'] = dict(gff=args.gff, ref=args.ref, target=target, track=track, biotype=biotype,
                       min_species=args.min_species, strand=args.strand, keep_self=args.keep_self,
                       divergent_dist=args.divergent_dist)
    with open(args.prefix + '.curation.tsv', 'w') as o:
        o.write('table\tmetric\tvalue\n')
        for t in ('observed', 'discarded_by_reference', 'novel'):
            for k, v in R[t].items():
                o.write(f'{t}\t{k}\t{v}\n')
    with open(args.prefix + '.curation.json', 'w') as o:
        json.dump(R, o, indent=1)
    log(f"supported chains {obs['supported_chains']:,} -> {len(sup_post):,} after reference filter; "
        f"{len(genes):,} curated genes -> {args.prefix}.gff3.gz in {time.time() - t0:.0f}s")


if __name__ == '__main__':
    main()
