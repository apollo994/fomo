#!/usr/bin/env python3
"""
Shared helpers for the FOMO bin/ scripts that read gffcompare `.stats` files:
the model-identity filename grammar, the aggregate pseudo-source ids, and the
Sensitivity/Precision parser.

This is a MODULE, not a script — it is imported by `select_top_sources.py`,
`gffcompare_accuracy_mqc.py` and `run_summary_tables.py`, each of which does

    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import fomo_stats

Nextflow only puts `bin/` on PATH, not on Python's import path — but it stages
(or bind-mounts) the whole directory, so a `__file__`-relative insert finds the
siblings on every executor. Do not "simplify" that preamble to a bare import.

Why this file exists: the identity regex, AGGREGATE_IDS and the stats parser were
duplicated across two scripts that MUST agree with each other and with the channel
filter in `subworkflows/local/consensus_top.nf`. One copy here leaves exactly two
places to keep in sync — this module and that filter.
"""
import os
import re
from typing import Dict, NamedTuple, Optional, Tuple

# Feature levels, in the order gffcompare reports them.
LEVELS = ["Base", "Exon", "Intron", "Intron_chain", "Transcript", "Locus"]

# Pseudo-source ids standing for an aggregate model rather than a real species. They
# use the same `from_<id>` filename convention as per-source stats, so a ranking pool
# must skip them or the top-N would be picked from models built out of the top-N.
# Kept in sync with the channel filter in subworkflows/local/consensus_top.nf, which
# filters before this code ever sees a file — this is the second layer.
AGGREGATE_IDS = ("allModels_raw", "allModels_collapsed", "top3_raw", "top3_collapsed")

# <target>.from_<source>.<feature_type>[.decoy], after the .gffcompare.stats suffix
# is stripped. Species names contain '_' and digits, hence the non-greedy captures
# anchored on the literal '.from_' and the feature-type alternation.
FNAME_RE = re.compile(
    r"^(?P<target>.+?)\.from_(?P<source>.+?)\.(?P<ft>lnc_RNA|mRNA)(?P<decoy>\.decoy)?$"
)


class Identity(NamedTuple):
    """Which model a gffcompare stats file describes."""
    target: str
    source: str          # a real species, or one of AGGREGATE_IDS
    feature_type: str    # 'lnc_RNA' | 'mRNA'
    decoy: bool

    @property
    def is_aggregate(self) -> bool:
        return self.source in AGGREGATE_IDS

    @property
    def is_self(self) -> bool:
        """A species projected onto itself — the Sn/Pr ceiling control, excluded
        from every ranking pool (see consensus_top.nf's meta.id != meta.target_id)."""
        return self.source == self.target


def parse_identity(stats_path: str) -> Optional[Identity]:
    """Derive the Identity from a gffcompare stats filename, or None if it does
    not match the grammar."""
    base = os.path.basename(stats_path)
    base = re.sub(r"\.gffcompare\.stats$", "", base)
    base = re.sub(r"\.stats$", "", base)            # tolerate either suffix
    base = re.sub(r"\.gffcompare$", "", base)
    m = FNAME_RE.match(base)
    if not m:
        return None
    return Identity(
        target=m.group("target"),
        source=m.group("source"),
        feature_type=m.group("ft"),
        decoy=bool(m.group("decoy")),
    )


def feature_class(ft: str, decoy: bool) -> str:
    """Presentation label for a (feature_type, decoy) pair."""
    if decoy:
        return "decoy"
    return "lncRNA" if ft == "lnc_RNA" else "mRNA"


def parse_accuracy(stats_path: str) -> Dict[str, Tuple[float, float]]:
    """Return {level: (sensitivity, precision)} read from a stats file.
    Mirrors MultiQC's native parser: drop '|', normalise 'Intron chain',
    then split — token[0]=level, token[2]=sensitivity, token[3]=precision."""
    out: Dict[str, Tuple[float, float]] = {}
    with open(stats_path, encoding="utf-8") as fh:
        for line in fh:
            if "level:" not in line:
                continue
            toks = line.replace("|", "").replace("Intron chain", "Intron_chain").split()
            if len(toks) < 4:
                continue
            level = toks[0]
            try:
                sens = float(toks[2])
                prec = float(toks[3])
            except ValueError:
                continue
            out[level] = (sens, prec)
    return out


def transcript_sn_pr(stats_path: str) -> Optional[Tuple[float, float]]:
    """(sensitivity, precision) from the Transcript-level row, or None if the file
    has no parseable one. Transcript level is the ranking level throughout FOMO."""
    return parse_accuracy(stats_path).get("Transcript")


def f1(sn: float, pr: float) -> float:
    return 0.0 if (sn + pr) == 0 else 2.0 * sn * pr / (sn + pr)
