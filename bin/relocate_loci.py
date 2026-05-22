#!/usr/bin/env python3
import argparse
import random
import sys
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional


# ----------------------------
# Data structures
# ----------------------------
@dataclass
class Feature:
    seqid: str
    source: str
    ftype: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attrs_raw: str
    attrs: Dict[str, str]

    def to_gff3(self) -> str:
        attr_str = format_gff3_attrs(self.attrs)
        return "\t".join([
            self.seqid,
            self.source,
            self.ftype,
            str(self.start),
            str(self.end),
            self.score,
            self.strand,
            self.phase,
            attr_str
        ])


@dataclass
class TranscriptModel:
    tid: str
    transcript: Feature
    exons: List[Feature]
    gene: Optional[Feature] = None


# ----------------------------
# GFF3 attribute parsing
# ----------------------------
def parse_gff3_attrs(attr_str: str) -> Dict[str, str]:
    d: Dict[str, str] = {}
    attr_str = attr_str.strip()
    if not attr_str or attr_str == ".":
        return d
    for part in attr_str.split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            k = k.strip()
            v = v.strip()
            if k and k not in d:
                d[k] = v
        else:
            toks = part.split()
            if len(toks) >= 2 and toks[0] not in d:
                d[toks[0]] = toks[1]
    return d


def format_gff3_attrs(attrs: Dict[str, str]) -> str:
    if not attrs:
        return "."
    keys = list(attrs.keys())
    preferred = []
    for k in ["ID", "Parent", "Name"]:
        if k in attrs:
            preferred.append(k)
    rest = sorted([k for k in keys if k not in set(preferred)])
    out_keys = preferred + rest
    return ";".join([f"{k}={attrs[k]}" for k in out_keys])


# ----------------------------
# Intergenic interval handling
# ----------------------------
def bed_to_gff_interval(bed_start0: int, bed_end: int) -> Tuple[int, int]:
    return bed_start0 + 1, bed_end


def interval_len(iv: Tuple[int, int]) -> int:
    return iv[1] - iv[0] + 1


def remove_subinterval(intervals: List[Tuple[int, int]], idx: int, used: Tuple[int, int]) -> None:
    s, e = intervals[idx]
    us, ue = used
    new_ivs: List[Tuple[int, int]] = []
    if us > s:
        new_ivs.append((s, us - 1))
    if ue < e:
        new_ivs.append((ue + 1, e))
    intervals.pop(idx)
    for iv in new_ivs:
        intervals.append(iv)


def choose_region_for_length(intervals: List[Tuple[int, int]], length: int) -> Optional[int]:
    candidates = [i for i, iv in enumerate(intervals) if interval_len(iv) >= length]
    if not candidates:
        return None
    return random.choice(candidates)


# ----------------------------
# Model shape calculation and relocation
# ----------------------------
def transcript_span(exons: List[Feature]) -> Tuple[int, int]:
    return min(e.start for e in exons), max(e.end for e in exons)


def build_shape(exons: List[Feature]) -> List[Tuple[int, int]]:
    exons_sorted = sorted(exons, key=lambda f: f.start)
    t_start = exons_sorted[0].start
    return [(ex.start - t_start, ex.end - t_start) for ex in exons_sorted]


def relocate_model(model: TranscriptModel, new_seqid: str, new_t_start: int, new_source: str = "relocated") -> TranscriptModel:
    t_old_start, t_old_end = transcript_span(model.exons)
    t_len = t_old_end - t_old_start + 1
    new_t_end = new_t_start + t_len - 1

    t_feat = model.transcript
    new_t_attrs = dict(t_feat.attrs)
    old_tid = model.tid
    new_tid = f"{old_tid}.reloc"
    new_gid = f"{old_tid}.gene.reloc"
    new_t_attrs["ID"] = new_tid
    new_t_attrs["Parent"] = new_gid

    new_gene = Feature(
        seqid=new_seqid, source=new_source, ftype="gene",
        start=new_t_start, end=new_t_end,
        score=t_feat.score, strand=t_feat.strand, phase=".",
        attrs_raw="", attrs={"ID": new_gid, "decoy": "true"}
    )

    new_t_attrs["decoy"] = "true"
    new_transcript = Feature(
        seqid=new_seqid, source=new_source, ftype=t_feat.ftype,
        start=new_t_start, end=new_t_end,
        score=t_feat.score, strand=t_feat.strand, phase=t_feat.phase,
        attrs_raw=t_feat.attrs_raw, attrs=new_t_attrs
    )

    new_exons: List[Feature] = []
    for i, ex in enumerate(model.exons, start=1):
        ex_off_s = ex.start - t_old_start
        ex_off_e = ex.end - t_old_start
        new_ex_s = new_t_start + ex_off_s
        new_ex_e = new_t_start + ex_off_e

        new_e_attrs = dict(ex.attrs)
        new_e_attrs["Parent"] = new_tid
        if "ID" in new_e_attrs:
            new_e_attrs["ID"] = f"{new_e_attrs['ID']}.reloc"
        else:
            new_e_attrs["ID"] = f"{new_tid}.exon{i}"

        new_exons.append(Feature(
            seqid=new_seqid, source=new_source, ftype=ex.ftype,
            start=new_ex_s, end=new_ex_e,
            score=ex.score, strand=ex.strand, phase=ex.phase,
            attrs_raw=ex.attrs_raw, attrs=new_e_attrs
        ))

    new_exons.sort(key=lambda f: f.start)
    return TranscriptModel(tid=new_tid, transcript=new_transcript, exons=new_exons, gene=new_gene)


# ----------------------------
# Parsing input files
# ----------------------------
def read_intergenic_bed(path: str, min_length: int = 0) -> Dict[str, List[Tuple[int, int]]]:
    per_chr: Dict[str, List[Tuple[int, int]]] = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start0 = int(parts[1])
            end = int(parts[2])
            s1, e1 = bed_to_gff_interval(start0, end)
            if e1 < s1:
                continue
            if min_length > 0 and interval_len((s1, e1)) < min_length:
                continue
            per_chr.setdefault(chrom, []).append((s1, e1))
    for chrom in per_chr:
        per_chr[chrom].sort()
    return per_chr


def read_gff3_models(path: str, transcript_type: str = "lnc_RNA", exon_type: str = "exon") -> List[TranscriptModel]:
    transcripts: Dict[str, Feature] = {}
    exons_by_parent: Dict[str, List[Feature]] = {}

    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            seqid, source, ftype, start, end, score, strand, phase, attrs_raw = fields
            start_i = int(start)
            end_i = int(end)
            attrs = parse_gff3_attrs(attrs_raw)

            feat = Feature(
                seqid=seqid, source=source, ftype=ftype,
                start=start_i, end=end_i,
                score=score, strand=strand, phase=phase,
                attrs_raw=attrs_raw, attrs=attrs
            )

            if ftype == transcript_type:
                tid = attrs.get("ID")
                if tid:
                    transcripts[tid] = feat
            elif ftype == exon_type:
                parent = attrs.get("Parent")
                if parent:
                    for p in parent.split(","):
                        p = p.strip()
                        if p:
                            exons_by_parent.setdefault(p, []).append(feat)

    models: List[TranscriptModel] = []
    for tid, tfeat in transcripts.items():
        exs = exons_by_parent.get(tid, [])
        if len(exs) < 2:
            continue
        exs_sorted = sorted(exs, key=lambda f: f.start)
        models.append(TranscriptModel(tid=tid, transcript=tfeat, exons=exs_sorted))

    return models


# ----------------------------
# CLI
# ----------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Relocate transcript models into intergenic intervals of a target genome"
    )
    p.add_argument("--input-gff",            required=True,  metavar="GFF3",
                   help="Source filtered GFF3 (plain, not gzipped)")
    p.add_argument("--intergenic-bed",        required=True,  metavar="BED",
                   help="Target intergenic intervals in BED format")
    p.add_argument("--output-gff",           required=True,  metavar="GFF3",
                   help="Output decoy GFF3")
    p.add_argument("--seed",                 type=int, default=None,
                   help="Random seed for reproducible placement")
    p.add_argument("--feature-type",         default="lnc_RNA",
                   help="Transcript feature type to relocate (default: lnc_RNA)")
    p.add_argument("--exon-type",            default="exon",
                   help="Exon feature type (default: exon)")
    p.add_argument("--min-intergenic-length", type=int, default=0, metavar="BP",
                   help="Skip intergenic intervals shorter than this (default: 0)")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    intergenic = read_intergenic_bed(args.intergenic_bed, min_length=args.min_intergenic_length)
    models = read_gff3_models(args.input_gff, transcript_type=args.feature_type, exon_type=args.exon_type)
    random.shuffle(models)

    relocated: List[TranscriptModel] = []
    skipped = 0

    for m in models:
        old_span_s, old_span_e = transcript_span(m.exons)
        t_len = old_span_e - old_span_s + 1

        # Collect target chromosomes that have at least one interval large enough
        eligible_chroms = [
            chrom for chrom, ivs in intergenic.items()
            if any(interval_len(iv) >= t_len for iv in ivs)
        ]
        if not eligible_chroms:
            skipped += 1
            continue

        tgt_chrom = random.choice(eligible_chroms)
        idx = choose_region_for_length(intergenic[tgt_chrom], t_len)
        if idx is None:
            skipped += 1
            continue

        region_s, region_e = intergenic[tgt_chrom][idx]
        new_start = random.randint(region_s, region_e - t_len + 1)
        used = (new_start, new_start + t_len - 1)

        remove_subinterval(intergenic[tgt_chrom], idx, used)
        relocated.append(relocate_model(m, new_seqid=tgt_chrom, new_t_start=new_start, new_source="relocated"))

    with open(args.output_gff, "w", encoding="utf-8") as out:
        out.write("##gff-version 3\n")
        out.write(f"##source=relocate_loci\n")
        out.write(f"##feature_type={args.feature_type}\n")
        out.write(f"##relocated_models={len(relocated)}\n")
        out.write(f"##skipped_models={skipped}\n")
        for rm in relocated:
            if rm.gene is not None:
                out.write(rm.gene.to_gff3() + "\n")
            out.write(rm.transcript.to_gff3() + "\n")
            for ex in rm.exons:
                out.write(ex.to_gff3() + "\n")

    print(f"Relocated models: {len(relocated)}", file=sys.stderr)
    print(f"Skipped models:   {skipped}", file=sys.stderr)

    if not relocated:
        print("WARNING: zero models were relocated — check input GFF3 and intergenic BED", file=sys.stderr)


if __name__ == "__main__":
    main()
