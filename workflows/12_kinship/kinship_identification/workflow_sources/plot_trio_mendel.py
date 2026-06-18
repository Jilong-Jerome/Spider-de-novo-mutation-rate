#!/usr/bin/env python3
"""
plot_trio_mendel.py

Figure for the per-trio Mendelian validation produced by trio_mendel_check.py.
Reads the trio TSV and writes a 3-panel PDF:

  Panel A: per-offspring consistency (fraction of parental opposite-homozygote
           sites at which the offspring is heterozygous), one point per trio,
           grouped along the x-axis by family and coloured by verdict. A dashed
           guide line marks consistency_min. Valid trios sit near 1.0. The
           fraction is recomputed from the n_het/n_informative columns so that
           insufficient-sites trios (consistency = NA in the TSV) are shown too,
           in grey; a trio with zero informative sites is drawn at y=0 with an
           'x' marker so it is never silently dropped.
  Panel B: n_informative (opposite-homozygote sites, log x) vs consistency,
           coloured by verdict -- shows whether any flagged trio is merely
           low-powered (few sites) rather than genuinely inconsistent.
  Panel C: histogram of consistency across all trios.

Usage:
    python plot_trio_mendel.py --input_tsv {SP}_trio_mendel.tsv \
        --output_pdf {SP}_trio_mendel.pdf --species MIM \
        --consistency_min 0.9
"""
import argparse
import csv

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

VERDICT_COLORS = {
    'consistent': '#2ca02c',
    'INCONSISTENT': '#d62728',
    'insufficient-sites': '#9e9e9e',
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

    Derived from n_het / n_informative rather than read from the 'consistency'
    field, so insufficient-sites trios (which carry consistency = NA in the TSV
    because n_informative < min_sites) are still plotted. Returns nan only when
    n_informative == 0 (no opposite-homozygote site at all -> no fraction).
    """
    n = _f(r['n_informative'])
    h = _f(r['n_het'])
    if n != n or h != h or n <= 0:
        return float('nan')
    return h / n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input_tsv', required=True)
    ap.add_argument('--output_pdf', required=True)
    ap.add_argument('--species', default='')
    ap.add_argument('--consistency_min', type=float, default=0.9)
    args = ap.parse_args()

    rows = load(args.input_tsv)

    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(18, 5.5))

    families = sorted({r['family'] for r in rows})
    fam_x = {fam: i for i, fam in enumerate(families)}

    # ---- Panel A: consistency per trio, grouped by family ----
    # y is the empirical het fraction (n_het/n_informative); insufficient-sites
    # trios are plotted too (grey). Trios with no informative site at all cannot
    # have a fraction -- they are drawn at y=0 with an 'x' marker so they are
    # still visible rather than silently dropped.
    rng = np.random.default_rng(0)
    n_zero_info = 0
    for r in rows:
        x = fam_x[r['family']] + (rng.random() - 0.5) * 0.5  # jitter within family
        y = het_fraction(r)
        v = r['verdict']
        color = VERDICT_COLORS.get(v, '#000000')
        if y == y:
            axA.scatter(x, y, c=color, s=55, edgecolors='k', linewidths=0.3,
                        alpha=0.9, label=v)
        else:
            n_zero_info += 1
            axA.scatter(x, 0.0, c=color, s=55, marker='x', linewidths=1.2,
                        alpha=0.9, label='%s (0 sites)' % v)
    axA.axhline(args.consistency_min, color='grey', ls='--', lw=0.9,
                label='consistency_min=%.2f' % args.consistency_min)
    axA.set_xticks(range(len(families)))
    axA.set_xticklabels(families, rotation=45, ha='right', fontsize=8)
    axA.set_ylabel('offspring het fraction at parental\nopposite-homozygote sites')
    axA.set_ylim(-0.02, 1.05)
    axA.set_title('A. Trio consistency by family')
    # De-duplicate legend labels.
    handles, labels = axA.get_legend_handles_labels()
    seen = {}
    for h, l in zip(handles, labels):
        seen.setdefault(l, h)
    axA.legend(seen.values(), seen.keys(), fontsize=7, loc='lower right')

    # ---- Panel B: n_informative vs consistency ----
    # log x-axis cannot show n_informative == 0; those are annotated as a count.
    for v, color in VERDICT_COLORS.items():
        pts = [(_f(r['n_informative']), het_fraction(r))
               for r in rows if r['verdict'] == v]
        pts = [(x, y) for x, y in pts if x == x and x > 0 and y == y]
        if pts:
            xs, ys = zip(*pts)
            axB.scatter(xs, ys, c=color, s=50, edgecolors='k', linewidths=0.3,
                        alpha=0.9, label=v)
    if n_zero_info:
        axB.text(0.02, 0.02, '%d trio(s) with 0 informative sites (not shown)'
                 % n_zero_info, transform=axB.transAxes, fontsize=7,
                 color='#9e9e9e', va='bottom')
    axB.axhline(args.consistency_min, color='grey', ls='--', lw=0.9)
    axB.set_xscale('log')
    axB.set_xlabel('n_informative (opposite-homozygote sites, log)')
    axB.set_ylabel('consistency')
    axB.set_ylim(-0.02, 1.05)
    axB.set_title('B. Power vs consistency')
    axB.legend(fontsize=7, loc='lower right')

    # ---- Panel C: histogram of consistency ----
    vals = [het_fraction(r) for r in rows]
    vals = [v for v in vals if v == v]
    if vals:
        axC.hist(vals, bins=np.linspace(0, 1, 21), color='#4c72b0',
                 edgecolor='k', alpha=0.85)
    axC.axvline(args.consistency_min, color='red', ls='--', lw=0.9)
    axC.set_xlabel('consistency')
    axC.set_ylabel('number of trios')
    axC.set_title('C. Consistency distribution')

    title = 'Per-trio Mendelian validation (parental opposite-homozygote test)'
    if args.species:
        title += ' — %s' % args.species
    fig.suptitle(title, fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(args.output_pdf)
    print('Wrote %s' % args.output_pdf)


if __name__ == '__main__':
    main()
