#!/usr/bin/env python3
"""curate_models.py report (<prefix>.curation.json) -> three MultiQC custom-content tables (plans/21).

One row per curated track; Sample = <target>.<track> (track = <feature_type>[.decoy]):
  <prefix>_curation_observed_mqc.tsv   (1) multi-species chains observed, before the reference filter
  <prefix>_curation_discarded_mqc.tsv  (2) discarded by overlap with the reference, chains lost per
                                           primary reference type (with a reference only; otherwise NA)
  <prefix>_curation_novel_mqc.tsv      (3) the novel genes and transcripts (one per supported chain, plans/25):
                                           support, chains, ref_location, medians
                                           (incl. distance of intergenic genes to the nearest reference gene)

No `# id:` header: the sections are routed by `sp` filename patterns only (assets/multiqc/sections.yml),
since a file found by both renders twice. Stdlib only.
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fomo_stats  # noqa: E402

OBSERVED = fomo_stats.CURATION_OBSERVED
DISCARDED = fomo_stats.CURATION_DISCARDED
NOVEL = fomo_stats.CURATION_NOVEL
LOCATIONS = fomo_stats.CURATION_LOCATION_COLUMNS
MEDIANS = fomo_stats.CURATION_MEDIANS


def write(path, sample, cols):
    with open(path, 'w') as o:
        o.write('Sample\t' + '\t'.join(k for k, _ in cols) + '\n')
        o.write(sample + '\t' + '\t'.join('NA' if v is None else str(v) for _, v in cols) + '\n')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--json', required=True, help='<prefix>.curation.json from curate_models.py')
    p.add_argument('--prefix', required=True, help='output prefix (e.g. <target>.curated.lnc_RNA)')
    p.add_argument('--sample', help='MultiQC sample name (default: <target>.<track> from the report params)')
    args = p.parse_args()

    with open(args.json, encoding='utf-8') as fh:
        r = json.load(fh)
    obs, dis, nov, par = r['observed'], r['discarded_by_reference'], r['novel'], r['params']
    sample = args.sample or f"{par['target']}.{par['track']}"
    has_ref = bool(par.get('ref'))

    write(f'{args.prefix}_curation_observed_mqc.tsv', sample, [(k, obs.get(k, 0)) for k in OBSERVED])

    cols = [('reference', 'yes' if has_ref else 'no')]
    cols += [(k, dis.get(k, 0) if has_ref else None) for k in DISCARDED]
    pref = 'chains_lost_by_primary_type:'
    cols += [(f'chains_lost:{k[len(pref) :]}', v) for k, v in dis.items() if k.startswith(pref)]
    write(f'{args.prefix}_curation_discarded_mqc.tsv', sample, cols)

    cols = [(k, nov.get(k, 0)) for k in NOVEL]
    cols += [(f'location:{c}', nov.get(f'genes_location:{c}', 0)) for c in LOCATIONS]
    cols += [(f'median_{m}', nov.get(f'{m}_median')) for m in MEDIANS]
    cols += [('median_intergenic_distance', nov.get('intergenic_distance_median') if has_ref else None)]
    write(f'{args.prefix}_curation_novel_mqc.tsv', sample, cols)


if __name__ == '__main__':
    main()
