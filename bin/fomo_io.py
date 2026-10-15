#!/usr/bin/env python3
"""
Shared I/O for the FOMO bin/ scripts: plain or gzipped text, GFF3 and FASTA reading, the GFF3
attribute parser, the feature-type vocabulary, and one error and logging convention.

A module, not a script. Importers put bin/ on the import path first, as for fomo_stats.py:

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import fomo_io
"""

import gzip
import logging
import os
import sys
from collections.abc import Callable, Iterator

# Feature types, compared in lower case (Sequence Ontology names, any capitalisation). A
# transcript is any type that is none of the others, so a reference's rarer RNA types (scRNA,
# misc_RNA, Y_RNA, ...) are transcripts too.
GENE_TYPES = frozenset({'gene', 'pseudogene', 'ncrna_gene'})
STRUCTURAL_TYPES = frozenset({'region', 'chromosome', 'scaffold', 'contig'})
CDS_LIKE_TYPES = frozenset({'cds', 'five_prime_utr', 'three_prime_utr', 'utr', 'start_codon', 'stop_codon'})
CHILD_TYPES = CDS_LIKE_TYPES | frozenset(
    {
        'exon',
        'intron',
        'polya_site',
        'polya_signal_sequence',
        'tss',
        'transcription_start_site',
        'sequence_feature',
        'selenocysteine',
        'stop_codon_redefined_as_selenocysteine',
    }
)


def feature_kind(ftype: str) -> str:
    """'gene', 'structural', 'child' or 'transcript'."""
    t = ftype.lower()
    if t in GENE_TYPES:
        return 'gene'
    if t in STRUCTURAL_TYPES:
        return 'structural'
    if t in CHILD_TYPES:
        return 'child'
    return 'transcript'


def is_transcript(ftype: str) -> bool:
    return feature_kind(ftype) == 'transcript'


class FomoError(Exception):
    """An input or consistency error: run() reports it and exits 1."""


def die(msg: str):
    raise FomoError(msg)


def script_name() -> str:
    return os.path.splitext(os.path.basename(sys.argv[0]))[0] or 'fomo'


class _StderrHandler(logging.StreamHandler):
    """Writes to whatever sys.stderr is when a record is emitted (tests may swap it)."""

    @property
    def stream(self):
        return sys.stderr

    @stream.setter
    def stream(self, _value):
        pass


def get_logger(script: str | None = None) -> logging.Logger:
    """A stderr logger, '<script>: <LEVEL>: <message>'."""
    logger = logging.getLogger(script or script_name())
    if not logger.handlers:
        handler = _StderrHandler()
        handler.setFormatter(logging.Formatter('%(name)s: %(levelname)s: %(message)s'))
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        logger.propagate = False
    return logger


def run(main: Callable[[], int | None]):
    """Run main() and exit with its status; a FomoError exits 1 with one logged message."""
    try:
        status = main()
    except FomoError as exc:
        get_logger().error(str(exc))
        sys.exit(1)
    sys.exit(status or 0)


def open_text(path: str, mode: str = 'rt'):
    """Open a text file for reading ('rt') or writing ('wt'), gzipped when the name ends in .gz."""
    if str(path).endswith('.gz'):
        return gzip.open(path, mode, encoding='utf-8')
    return open(path, mode.replace('t', ''), encoding='utf-8')


def iter_gff(path: str, strict: bool = False) -> Iterator[tuple[int, list[str]]]:
    """(line number, columns) of every feature line. Comments and blank lines are skipped, and
    reading stops at ##FASTA. A line with fewer than 9 tab-separated columns is skipped, or with
    strict=True is an error."""
    with open_text(path) as fh:
        for lineno, line in enumerate(fh, 1):
            if line.startswith('#'):
                if line.startswith('##FASTA'):
                    return
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) >= 9:
                yield lineno, cols
            elif strict and line.strip():
                die(f'{path}:{lineno}: not a 9-column GFF3 line')


def parse_attrs(col9: str) -> dict[str, str]:
    """GFF3 column 9 as a dict. Fields and their key and value are stripped; a field without '='
    is ignored; a repeated key keeps its first value. Values stay URL-escaped."""
    attrs: dict[str, str] = {}
    for field in col9.split(';'):
        key, sep, value = field.strip().partition('=')
        key = key.strip()
        if sep and key and key not in attrs:
            attrs[key] = value.strip()
    return attrs


def split_ids(value: str) -> list[str]:
    """The non-empty ids of a comma-separated attribute value (Parent=a,b)."""
    return [x for x in (v.strip() for v in value.split(',')) if x]


def iter_fasta(path: str) -> Iterator[tuple[str, str, list[str]]]:
    """(name, header line, sequence lines) of every record; the name is the header's first word.
    Lines keep their newlines, so a record is written back verbatim."""
    name = header = None
    seq: list[str] = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if header is not None:
                    yield name, header, seq
                header = line
                name = line[1:].split()[0] if line[1:].split() else ''
                seq = []
            elif header is not None:
                seq.append(line)
    if header is not None:
        yield name, header, seq


def fasta_names(path: str) -> list[str]:
    """Record names in file order."""
    return [name for name, _, _ in iter_fasta(path)]
