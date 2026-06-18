#!/usr/bin/env python3
"""
plot_combined_trio_mendel.py

Combined, all-species supplemental figure for the per-trio Mendelian validation
(parental opposite-homozygote test; see trio_mendel_check.py). Reads one
``{SP}_trio_mendel.tsv`` per species and writes a single multi-row PDF.

Layout: one row per species, two columns (the per-species histogram /
"distribution" panel of plot_trio_mendel.py is intentionally dropped here):

  Left column : per-offspring consistency grouped along the x-axis by family and
                coloured by verdict, with a dashed consistency_min guide line.
                As in plot_trio_mendel.py the y value is the empirical het
                fraction recomputed from n_het/n_informative, so
                insufficient-sites trios (consistency = NA in the TSV) are still
                shown (grey); a trio with zero informative sites is drawn at y=0
                with an 'x' marker.
  Right column: n_informative (opposite-homozygote sites, log x) vs the same het
                fraction, coloured by verdict -- shows whether a flagged trio is
                merely low-powered rather than genuinely inconsistent.

Species are ordered alphabetically by Latin binomial and titled in italics.
Rows/trios can be excluded with --exclude SP:family (repeatable); the supplement
excludes BIC family6.

Usage:
    python plot_combined_trio_mendel.py \
        --input_tsvs steps/AFR/AFR_trio_mendel.tsv steps/BIC/BIC_trio_mendel.tsv ... \
        --output_pdf steps/combined/all_species_trio_mendel.pdf \
        --consistency_min 0.9 --exclude BIC:family6
"""
import argparse
import csv

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

VERDICT_COLORS = {
    'consistent': '#2ca02c',
    'INCONSISTENT': '#d62728',
    'insufficient-sites': '#9e9e9e',
}

# Species code -> Latin binomial (rendered in italics as the row title).
SPECIES_LATIN = {
    'AFR': 'S. africanus',
    'BIC': 'S. bicolor',
    'DUM': 'S. dumicola',
    'LIN': 'S. lineatus',
    'MIM': 'S. mimosarum',
    'SAR': 'S. sarasinorum',
    'TEN': 'S. tentoriicola',
}


def load(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float('nan')


def het_fraction(r):
    """Empirical offspring het fraction, recomputed from the count columns.

    Derived from n_het / n_informative rather than the 'consistency' field, so
    insufficient-sites trios (consistency = NA in the TSV) are still plotted.
    Returns nan only when n_informative == 0 (no opposite-homozygote site).
    """
    n = _f(r['n_informative'])
    h = _f(r['n_het'])
    if n != n or h != h or n <= 0:
        return float('nan')
    return h / n


def family_num(fam):
    """Numeric suffix of a 'family{N}' label for sorting/x-ticks (else hash)."""
    digits = ''.join(c for c in fam if c.isdigit())
    return int(digits) if digits else fam


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input_tsvs', nargs='+', required=True,
                    help='one {SP}_trio_mendel.tsv per species (species read from '
                         'the TSV "species" column)')
    ap.add_argument('--output_pdf', required=True)
    ap.add_argument('--consistency_min', type=float, default=0.9)
    ap.add_argument('--exclude', action='append', default=[],
                    help='SP:family rows to drop, e.g. BIC:family6 (repeatable)')
    args = ap.parse_args()

    excluded = set()
    for token in args.exclude:
        sp, _, fam = token.partition(':')
        if sp and fam:
            excluded.add((sp, fam))

    # Collect rows per species (keyed by the TSV species column), dropping excludes.
    by_species = {}
    for path in args.input_tsvs:
        for r in load(path):
            sp = r['species']
            if (sp, r['family']) in excluded:
                continue
            by_species.setdefault(sp, []).append(r)

    # Order rows alphabetically by Latin binomial (fall back to the code).
    species = sorted(by_species, key=lambda s: SPECIES_LATIN.get(s, s))
    n = len(species)

    rng = np.random.default_rng(0)
    fig, axes = plt.subplots(n, 2, figsize=(11, 2.6 * n), squeeze=False)

    for i, sp in enumerate(species):
        rows = by_species[sp]
        axL, axR = axes[i][0], axes[i][1]
        latin = SPECIES_LATIN.get(sp, sp)

        families = sorted({r['family'] for r in rows}, key=family_num)
        fam_x = {fam: j for j, fam in enumerate(families)}

        # ---- Left: consistency grouped by family ----
        for r in rows:
            x = fam_x[r['family']] + (rng.random() - 0.5) * 0.5
            y = het_fraction(r)
            color = VERDICT_COLORS.get(r['verdict'], '#000000')
            if y == y:
                axL.scatter(x, y, c=color, s=42, edgecolors='k', linewidths=0.3,
                            alpha=0.9)
            else:
                axL.scatter(x, 0.0, c=color, s=42, marker='x', linewidths=1.2,
                            alpha=0.9)
        axL.axhline(args.consistency_min, color='grey', ls='--', lw=0.9)
        axL.set_xticks(range(len(families)))
        axL.set_xticklabels([str(family_num(f)) for f in families], fontsize=8)
        axL.set_xlim(-0.6, len(families) - 0.4)
        axL.set_ylim(-0.02, 1.05)
        axL.set_ylabel('het fraction', fontsize=9)
        axL.set_title('$%s$' % latin.replace(' ', r'\ '), fontsize=11, loc='left')
        if i == n - 1:
            axL.set_xlabel('family', fontsize=9)

        # ---- Right: power (n_informative, log x) vs consistency ----
        for v, color in VERDICT_COLORS.items():
            pts = [(_f(r['n_informative']), het_fraction(r))
                   for r in rows if r['verdict'] == v]
            pts = [(x, y) for x, y in pts if x == x and x > 0 and y == y]
            if pts:
                xs, ys = zip(*pts)
                axR.scatter(xs, ys, c=color, s=40, edgecolors='k', linewidths=0.3,
                            alpha=0.9)
        axR.axhline(args.consistency_min, color='grey', ls='--', lw=0.9)
        axR.set_xscale('log')
        axR.set_ylim(-0.02, 1.05)
        axR.set_ylabel('het fraction', fontsize=9)
        if i == n - 1:
            axR.set_xlabel('n_informative (opposite-homozygote sites, log)', fontsize=9)

    # Shared verdict legend at the figure level.
    handles = [Line2D([0], [0], marker='o', ls='', mec='k', mew=0.3,
                      mfc=VERDICT_COLORS[v], label=v) for v in VERDICT_COLORS]
    handles.append(Line2D([0], [0], color='grey', ls='--', lw=0.9,
                          label='consistency_min=%.2f' % args.consistency_min))
    fig.legend(handles=handles, loc='upper center', ncol=4, fontsize=9,
               frameon=False, bbox_to_anchor=(0.5, 1.0))

    fig.suptitle('Per-trio Mendelian validation across species '
                 '(parental opposite-homozygote test)', fontsize=13, y=1.005)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(args.output_pdf, bbox_inches='tight')
    print('Wrote %s (%d species, %d excluded rule(s))'
          % (args.output_pdf, n, len(excluded)))


if __name__ == '__main__':
    main()
