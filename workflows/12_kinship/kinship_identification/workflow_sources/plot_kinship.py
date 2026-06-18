#!/usr/bin/env python3
"""
plot_kinship.py

Validation figure for the pairwise kinship classification. Reads the classified
TSV produced by classify_kinship.py and writes a multi-panel PDF:

  Panel A: KING kinship coefficient vs IBS0, coloured by inferred class, marker
           shape by known class. Threshold guide lines drawn. This is the panel
           that visually separates 1st-degree pairs into parent-offspring
           (IBS0 ~ 0) and full-sibling (IBS0 > 0).
  Panel B: PLINK IBD probabilities Z0 vs Z1 (independent cross-check); the
           parent-offspring corner (Z1~1) and full-sib region (Z0~0.25,Z1~0.5)
           are visible.
  Panel C: Confusion matrix (inferred vs known) as an annotated heatmap.

Usage:
    python plot_kinship.py --input_tsv {SP}_kinship_classified.tsv \
        --output_pdf {SP}_kinship_validation.pdf --species MIM \
        --kinship_first 0.177 --kinship_second 0.0884 --kinship_third 0.0442 \
        --ibs0_po_max 0.005
"""
import argparse
import csv

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

CLASS_ORDER = ['unrelated', '3rd-degree', '2nd-degree',
               'full-sibling', 'parent-offspring', 'duplicate/MZ']
CLASS_COLORS = {
    'unrelated': '#9e9e9e',
    '3rd-degree': '#80cbc4',
    '2nd-degree': '#4db6ac',
    'full-sibling': '#1f77b4',
    'parent-offspring': '#d62728',
    'duplicate/MZ': '#000000',
    'NA': '#dddddd',
}
KNOWN_MARKERS = {
    'unrelated': 'o',
    'parent-offspring': '^',
    'full-sibling': 's',
    'unknown': 'x',
}


def load(path):
    rows = []
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for r in reader:
            rows.append(r)
    return rows


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float('nan')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input_tsv', required=True)
    ap.add_argument('--output_pdf', required=True)
    ap.add_argument('--species', default='')
    ap.add_argument('--kinship_first', type=float, default=0.177)
    ap.add_argument('--kinship_second', type=float, default=0.0884)
    ap.add_argument('--kinship_third', type=float, default=0.0442)
    ap.add_argument('--ibs0_po_max', type=float, default=0.005)
    args = ap.parse_args()

    rows = load(args.input_tsv)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    axA, axB, axC = axes

    # ---- Panel A: kinship vs IBS0 ----
    for inf_class in CLASS_ORDER + ['NA']:
        for known, marker in KNOWN_MARKERS.items():
            xs = [_f(r['IBS0']) for r in rows
                  if r['inferred_class'] == inf_class and r['known_class'] == known]
            ys = [_f(r['kinship']) for r in rows
                  if r['inferred_class'] == inf_class and r['known_class'] == known]
            if xs:
                axA.scatter(xs, ys, c=CLASS_COLORS.get(inf_class, '#dddddd'),
                            marker=marker, s=45, edgecolors='k', linewidths=0.3,
                            alpha=0.85)
    for y in (args.kinship_first, args.kinship_second, args.kinship_third):
        axA.axhline(y, color='grey', ls='--', lw=0.7)
    axA.axvline(args.ibs0_po_max, color='red', ls=':', lw=0.8)
    axA.set_xlabel('IBS0 (proportion, KING)')
    axA.set_ylabel('Kinship coefficient (KING)')
    axA.set_title('A. Kinship vs IBS0')

    # ---- Panel B: Z0 vs Z1 ----
    for inf_class in CLASS_ORDER + ['NA']:
        xs = [_f(r['Z0']) for r in rows if r['inferred_class'] == inf_class]
        ys = [_f(r['Z1']) for r in rows if r['inferred_class'] == inf_class]
        if xs:
            axB.scatter(xs, ys, c=CLASS_COLORS.get(inf_class, '#dddddd'),
                        marker='o', s=40, edgecolors='k', linewidths=0.3,
                        alpha=0.85, label=inf_class)
    axB.set_xlabel('Z0 (P(IBD=0), PLINK)')
    axB.set_ylabel('Z1 (P(IBD=1), PLINK)')
    axB.set_title('B. PLINK IBD probabilities')
    axB.set_xlim(-0.05, 1.05)
    axB.set_ylim(-0.05, 1.05)
    axB.legend(fontsize=7, loc='upper right', title='inferred')

    # ---- Panel C: confusion matrix ----
    known_labels = ['unrelated', 'parent-offspring', 'full-sibling']
    inferred_labels = ['unrelated', '3rd-degree', '2nd-degree',
                       'full-sibling', 'parent-offspring', 'duplicate/MZ']
    mat = np.zeros((len(known_labels), len(inferred_labels)), dtype=int)
    ki = {k: i for i, k in enumerate(known_labels)}
    ii = {k: i for i, k in enumerate(inferred_labels)}
    for r in rows:
        if r['known_class'] in ki and r['inferred_class'] in ii:
            mat[ki[r['known_class']], ii[r['inferred_class']]] += 1
    im = axC.imshow(mat, cmap='Blues', aspect='auto')
    axC.set_xticks(range(len(inferred_labels)))
    axC.set_xticklabels(inferred_labels, rotation=45, ha='right', fontsize=8)
    axC.set_yticks(range(len(known_labels)))
    axC.set_yticklabels(known_labels, fontsize=8)
    axC.set_xlabel('inferred class')
    axC.set_ylabel('known class')
    axC.set_title('C. Confusion (counts)')
    for i in range(len(known_labels)):
        for j in range(len(inferred_labels)):
            if mat[i, j]:
                axC.text(j, i, str(mat[i, j]), ha='center', va='center',
                         color='black', fontsize=8)
    fig.colorbar(im, ax=axC, fraction=0.046, pad=0.04)

    title = 'Pairwise kinship validation'
    if args.species:
        title += ' — %s' % args.species
    fig.suptitle(title, fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(args.output_pdf)
    print('Wrote %s' % args.output_pdf)


if __name__ == '__main__':
    main()
