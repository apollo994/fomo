#!/usr/bin/env python3
import sys
import random
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
    exons: List[Feature]  # sorted by start
    gene: Optional[Feature] = None


# ----------------------------
# GFF3 attribute parsing
# ----------------------------
def parse_gff3_attrs(attr_str: str) -> Dict[str, str]:
    """Parse GFF3 9th column key=value;key2=value2 into dict.
    Keeps first value for repeated keys; keeps raw for writing later."""
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
            # tolerate non-standard "key value"
            toks = part.split()
            if len(toks) >= 2 and toks[0] not in d:
                d[toks[0]] = toks[1]
    return d


def format_gff3_attrs(attrs: Dict[str, str]) -> str:
    if not attrs:
        return "."
    # stable order: ID, Parent, then the rest
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
# We store per-chrom intervals as a list of (start,end), 0-based half-open or 1-based?
# BED is 0-based half-open, GFF3 is 1-based inclusive.
# We'll convert BED to 1-based inclusive for easier use with GFF coordinates.
def bed_to_gff_interval(bed_start0: int, bed_end: int) -> Tuple[int, int]:
    # BED: [start0, end) => GFF: [start1, end1] inclusive
    start1 = bed_start0 + 1
    end1 = bed_end
    return start1, end1


def interval_len(iv: Tuple[int, int]) -> int:
    # 1-based inclusive
    return iv[1] - iv[0] + 1


def remove_subinterval(intervals: List[Tuple[int, int]], idx: int, used: Tuple[int, int]) -> None:
    """Remove 'used' (inclusive) from intervals[idx], splitting if needed."""
    s, e = intervals[idx]
    us, ue = used
    # used must lie within [s,e]
    new_ivs: List[Tuple[int, int]] = []
    if us > s:
        new_ivs.append((s, us - 1))
    if ue < e:
        new_ivs.append((ue + 1, e))
    # replace the interval at idx with 0,1, or 2 new intervals
    intervals.pop(idx)
    # insert back (keep list unsorted is ok; we can sort occasionally)
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
    """Return exon blocks as offsets from transcript start.
    Each block is (offset_start, offset_end) in bp, 0-based within transcript,
    but we will use them to create new 1-based genomic coordinates."""
    exons_sorted = sorted(exons, key=lambda f: f.start)
    t_start = exons_sorted[0].start
    blocks = []
    for ex in exons_sorted:
        blocks.append((ex.start - t_start, ex.end - t_start))
    return blocks


def relocate_model(model: TranscriptModel, new_seqid: str, new_t_start: int, new_source: str = "relocated") -> TranscriptModel:
    """Create a relocated copy of transcript+exons at (new_seqid, new_t_start)."""
    blocks = build_shape(model.exons)
    t_old_start, t_old_end = transcript_span(model.exons)
    t_len = t_old_end - t_old_start + 1

    # transcript feature spans entire exon span
    new_t_end = new_t_start + t_len - 1

    # clone transcript
    t_feat = model.transcript
    new_t_attrs = dict(t_feat.attrs)
    # Keep original ID but tag it; OR create a new ID to avoid collisions.
    # Safer: new ID
    old_tid = model.tid
    new_tid = f"{old_tid}.reloc"
    new_gid = f"{old_tid}.gene.reloc"
    new_t_attrs["ID"] = new_tid
    new_t_attrs["Parent"] = new_gid

    new_gene = Feature(
        seqid=new_seqid,
        source=new_source,
        ftype="gene_decoy",
        start=new_t_start,
        end=new_t_end,
        score=t_feat.score,
        strand=t_feat.strand,
        phase=".",
        attrs_raw="",
        attrs={"ID": new_gid}
    )

    new_transcript = Feature(
        seqid=new_seqid,
        source=new_source,
        ftype=t_feat.ftype + "_decoy",
        start=new_t_start,
        end=new_t_end,
        score=t_feat.score,
        strand=t_feat.strand,
        phase=t_feat.phase,
        attrs_raw=t_feat.attrs_raw,
        attrs=new_t_attrs
    )

    # clone exons
    new_exons: List[Feature] = []
    for i, ex in enumerate(model.exons, start=1):
        ex_off_s = ex.start - t_old_start
        ex_off_e = ex.end - t_old_start
        new_ex_s = new_t_start + ex_off_s
        new_ex_e = new_t_start + ex_off_e

        new_e_attrs = dict(ex.attrs)
        # update Parent to new transcript ID
        new_e_attrs["Parent"] = new_tid
        # update exon ID if present
        if "ID" in new_e_attrs:
            new_e_attrs["ID"] = f"{new_e_attrs['ID']}.reloc"
        else:
            new_e_attrs["ID"] = f"{new_tid}.exon{i}"

        new_exons.append(Feature(
            seqid=new_seqid,
            source=new_source,
            ftype=ex.ftype,
            start=new_ex_s,
            end=new_ex_e,
            score=ex.score,
            strand=ex.strand,
            phase=ex.phase,
            attrs_raw=ex.attrs_raw,
            attrs=new_e_attrs
        ))

    new_exons.sort(key=lambda f: f.start)
    return TranscriptModel(tid=new_tid, transcript=new_transcript, exons=new_exons, gene=new_gene)


# ----------------------------
# Parsing input files
# ----------------------------
def read_intergenic_bed(path: str) -> Dict[str, List[Tuple[int, int]]]:
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
            per_chr.setdefault(chrom, []).append((s1, e1))
    # sort for nicer behavior
    for chrom in per_chr:
        per_chr[chrom].sort()
    return per_chr


def read_gff3_models(path: str, transcript_type: str = "lnc_RNA", exon_type: str = "exon") -> List[TranscriptModel]:
    # First pass: collect transcript features and exon features by transcript Parent
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
                seqid=seqid,
                source=source,
                ftype=ftype,
                start=start_i,
                end=end_i,
                score=score,
                strand=strand,
                phase=phase,
                attrs_raw=attrs_raw,
                attrs=attrs
            )

            if ftype == transcript_type:
                tid = attrs.get("ID")
                if tid:
                    transcripts[tid] = feat
            elif ftype == exon_type:
                parent = attrs.get("Parent")
                if parent:
                    # Parent can be comma-separated in GFF3
                    for p in parent.split(","):
                        p = p.strip()
                        if not p:
                            continue
                        exons_by_parent.setdefault(p, []).append(feat)

    models: List[TranscriptModel] = []
    for tid, tfeat in transcripts.items():
        exs = exons_by_parent.get(tid, [])
        if len(exs) < 2:
            # keep only spliced models (as your original script did)
            continue
        exs_sorted = sorted(exs, key=lambda f: f.start)
        models.append(TranscriptModel(tid=tid, transcript=tfeat, exons=exs_sorted))

    return models


# ----------------------------
# Main
# ----------------------------
def main():
    if len(sys.argv) < 4:
        print(f"Usage: {sys.argv[0]} <intergenic.bed> <genes.gff3> <out.gff3> [seed]", file=sys.stderr)
        sys.exit(1)

    bed_path = sys.argv[1]
    gff_path = sys.argv[2]
    out_path = sys.argv[3]
    seed = int(sys.argv[4]) if len(sys.argv) >= 5 else None
    if seed is not None:
        random.seed(seed)

    intergenic = read_intergenic_bed(bed_path)
    models = read_gff3_models(gff_path, transcript_type="lnc_RNA", exon_type="exon")
    random.shuffle(models)

    relocated: List[TranscriptModel] = []
    skipped = 0

    for m in models:
        chrom = m.transcript.seqid
        if chrom not in intergenic or not intergenic[chrom]:
            skipped += 1
            continue

        # determine transcript length as span of exons (like original)
        old_span_s, old_span_e = transcript_span(m.exons)
        t_len = old_span_e - old_span_s + 1

        idx = choose_region_for_length(intergenic[chrom], t_len)
        if idx is None:
            skipped += 1
            continue

        region_s, region_e = intergenic[chrom][idx]
        # choose new transcript start so [new_start, new_start+t_len-1] fits
        new_start = random.randint(region_s, region_e - t_len + 1)
        used = (new_start, new_start + t_len - 1)

        # update intergenic space to avoid overlaps
        remove_subinterval(intergenic[chrom], idx, used)

        relocated.append(relocate_model(m, new_seqid=chrom, new_t_start=new_start, new_source="relocated"))

    # Write output: only relocated models (transcript + exons)
    with open(out_path, "w", encoding="utf-8") as out:
        out.write("##gff-version 3\n")
        out.write(f"##source=relocate_lncRNA\n")
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
    print(f"Wrote: {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
