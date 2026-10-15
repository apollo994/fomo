#!/usr/bin/env python3
"""
Relocate spliced transcript models into intergenic intervals: the decoy track.

Each decoy keeps a source model's shape (exon and intron lengths, strand) and moves to a random
position inside one intergenic interval, chosen with a chromosome probability proportional to its
eligible intergenic length (intervals at least as long as the model). A used stretch leaves its
interval, so decoys never overlap. Models with fewer than two exons are not relocated; with
--decoy-cap N a random sample of N models is. The same --seed gives the same decoys.

IDs: transcript <tid>.reloc, gene <tid>.gene.reloc, exon <exon ID>.reloc (or <tid>.reloc.exon<i>);
every gene and transcript carries decoy=true. The output is grouped by seqid, as gff-feature-stats
requires, with a header that records the relocated and skipped counts.

Usage: relocate_loci.py --gff models.gff3 --bed intergenic.bed --output decoys.gff3 [--seed N]
"""

import argparse
import random
import sys
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fomo_io  # noqa: E402

Interval = tuple[int, int]
log = fomo_io.get_logger('relocate_loci')


@dataclass(slots=True)
class Feature:
    seqid: str
    source: str
    ftype: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attrs: dict[str, str]

    def to_gff3(self) -> str:
        cols = [self.seqid, self.source, self.ftype, str(self.start), str(self.end)]
        return '\t'.join(cols + [self.score, self.strand, self.phase, format_attrs(self.attrs)])


@dataclass(slots=True)
class TranscriptModel:
    tid: str
    transcript: Feature
    exons: list[Feature]
    gene: Feature | None = None


def format_attrs(attrs: dict[str, str]) -> str:
    """ID, Parent and Name first, then the other keys sorted; '.' when empty."""
    if not attrs:
        return '.'
    first = [k for k in ('ID', 'Parent', 'Name') if k in attrs]
    keys = first + sorted(k for k in attrs if k not in first)
    return ';'.join(f'{k}={attrs[k]}' for k in keys)


def interval_len(iv: Interval) -> int:
    return iv[1] - iv[0] + 1


def remove_subinterval(intervals: list[Interval], idx: int, used: Interval) -> None:
    """Replace interval idx by what remains of it around `used` (appended at the end)."""
    s, e = intervals.pop(idx)
    if used[0] > s:
        intervals.append((s, used[0] - 1))
    if used[1] < e:
        intervals.append((used[1] + 1, e))


def choose_region(intervals: list[Interval], length: int, rng: random.Random) -> int | None:
    """Index of a random interval at least `length` long, or None."""
    candidates = [i for i, iv in enumerate(intervals) if interval_len(iv) >= length]
    return rng.choice(candidates) if candidates else None


def transcript_span(exons: list[Feature]) -> Interval:
    return min(e.start for e in exons), max(e.end for e in exons)


def relocate_model(model: TranscriptModel, seqid: str, start: int, source: str = 'relocated') -> TranscriptModel:
    """The model moved to seqid:start, with decoy IDs and attributes."""
    old_start, old_end = transcript_span(model.exons)
    end = start + old_end - old_start
    t = model.transcript
    tid, gid = f'{model.tid}.reloc', f'{model.tid}.gene.reloc'
    gene = Feature(seqid, source, 'gene', start, end, t.score, t.strand, '.', {'ID': gid, 'decoy': 'true'})
    t_attrs = dict(t.attrs) | {'ID': tid, 'Parent': gid, 'decoy': 'true'}
    transcript = Feature(seqid, source, t.ftype, start, end, t.score, t.strand, t.phase, t_attrs)
    exons = []
    for i, ex in enumerate(model.exons, start=1):
        e_attrs = dict(ex.attrs) | {'Parent': tid}
        e_attrs['ID'] = f'{e_attrs["ID"]}.reloc' if 'ID' in e_attrs else f'{tid}.exon{i}'
        shift = start - old_start
        exons.append(
            Feature(seqid, source, ex.ftype, ex.start + shift, ex.end + shift, ex.score, ex.strand, ex.phase, e_attrs)
        )
    exons.sort(key=lambda f: f.start)
    return TranscriptModel(tid=tid, transcript=transcript, exons=exons, gene=gene)


def read_intergenic_bed(path: str, min_length: int = 0) -> dict[str, list[Interval]]:
    """Chromosome -> sorted 1-based closed intervals of a BED file."""
    per_chr: dict[str, list[Interval]] = {}
    with fomo_io.open_text(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 3 or line.startswith('#'):
                continue
            iv = (int(parts[1]) + 1, int(parts[2]))
            if iv[1] < iv[0] or (min_length > 0 and interval_len(iv) < min_length):
                continue
            per_chr.setdefault(parts[0], []).append(iv)
    for ivs in per_chr.values():
        ivs.sort()
    return per_chr


def read_gff3_models(path: str, transcript_type: str, exon_type: str) -> list[TranscriptModel]:
    """The transcripts of `transcript_type` with at least two exons, in file order."""
    transcripts: dict[str, Feature] = {}
    exons_by_parent: dict[str, list[Feature]] = {}
    for _, cols in fomo_io.iter_gff(path):
        seqid, source, ftype, start, end, score, strand, phase, col9 = cols[:9]
        attrs = fomo_io.parse_attrs(col9)
        feat = Feature(seqid, source, ftype, int(start), int(end), score, strand, phase, attrs)
        if ftype == transcript_type and attrs.get('ID'):
            transcripts[attrs['ID']] = feat
        elif ftype == exon_type:
            for parent in fomo_io.split_ids(attrs.get('Parent', '')):
                exons_by_parent.setdefault(parent, []).append(feat)
    models = []
    for tid, tfeat in transcripts.items():
        exons = exons_by_parent.get(tid, [])
        if len(exons) >= 2:
            models.append(TranscriptModel(tid=tid, transcript=tfeat, exons=sorted(exons, key=lambda f: f.start)))
    return models


def relocate_all(
    models: list[TranscriptModel], intergenic: dict[str, list[Interval]], rng: random.Random
) -> tuple[list[TranscriptModel], int]:
    """-> (relocated models grouped by seqid, number skipped for lack of space)."""
    relocated, skipped = [], 0
    for model in models:
        span_start, span_end = transcript_span(model.exons)
        length = span_end - span_start + 1
        weights = []
        for chrom, ivs in intergenic.items():
            eligible = sum(interval_len(iv) for iv in ivs if interval_len(iv) >= length)
            if eligible > 0:
                weights.append((chrom, eligible))
        if not weights:
            skipped += 1
            continue
        chroms, chrom_weights = zip(*weights, strict=True)
        chrom = rng.choices(chroms, weights=chrom_weights, k=1)[0]
        idx = choose_region(intergenic[chrom], length, rng)
        if idx is None:
            skipped += 1
            continue
        region_start, region_end = intergenic[chrom][idx]
        start = rng.randint(region_start, region_end - length + 1)
        remove_subinterval(intergenic[chrom], idx, (start, start + length - 1))
        relocated.append(relocate_model(model, chrom, start))
    # Each model picks its chromosome independently; a stable sort on seqid groups the output
    # (gff-feature-stats rejects interleaved seqids) and keeps every model's lines together.
    relocated.sort(key=lambda m: m.transcript.seqid)
    return relocated, skipped


def write_decoys(path: str, relocated: list[TranscriptModel], skipped: int, feature_type: str) -> None:
    with fomo_io.open_text(path, 'wt') as out:
        out.write('##gff-version 3\n##source=relocate_loci\n')
        out.write(f'##feature_type={feature_type}\n##relocated_models={len(relocated)}\n##skipped_models={skipped}\n')
        for model in relocated:
            out.write(model.gene.to_gff3() + '\n')
            out.write(model.transcript.to_gff3() + '\n')
            out.writelines(ex.to_gff3() + '\n' for ex in model.exons)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--gff', required=True, help='GFF3 of the models to relocate')
    p.add_argument('--bed', required=True, help='intergenic intervals (BED)')
    p.add_argument('--output', required=True, help='decoy GFF3')
    p.add_argument('--seed', type=int, default=None, help='random seed for reproducible placement')
    p.add_argument('--feature-type', default='lnc_RNA', help='transcript feature type to relocate (lnc_RNA)')
    p.add_argument('--exon-type', default='exon', help='exon feature type (exon)')
    p.add_argument('--min-intergenic-length', type=int, default=0, help='skip shorter intervals (bp, 0)')
    p.add_argument('--decoy-cap', type=int, default=1000, help='relocate at most N random models (1000; 0: all)')
    return p.parse_args()


def main() -> int:
    args = parse_args()
    rng = random.Random(args.seed)
    intergenic = read_intergenic_bed(args.bed, min_length=args.min_intergenic_length)
    models = read_gff3_models(args.gff, args.feature_type, args.exon_type)
    rng.shuffle(models)
    if args.decoy_cap > 0:
        models = models[: args.decoy_cap]
    relocated, skipped = relocate_all(models, intergenic, rng)
    write_decoys(args.output, relocated, skipped, args.feature_type)
    log.info(f'relocated {len(relocated)} models, skipped {skipped}')
    if not relocated:
        log.warning('no model was relocated: check the input GFF3 and the intergenic BED')
    return 0


if __name__ == '__main__':
    fomo_io.run(main)
