#!/usr/bin/env python3
"""
plot_phylo_aligned_rates.py

Per-species overall germline + somatic mutation rates aligned to the species
phylogeny. Produces:

  * TSV: per-species germline rate point estimate, bootstrap 95% CI (resampling
    trios with replacement, n_bootstrap iterations), somatic point estimate,
    and somatic/germline fold change.
  * PDF: 4-panel figure (1 row, 4 columns), shared species y-axis. Panel A is
    the species cladogram; panels B/C/D are germline rate (point + 95% CI
    segment), somatic rate (point), and fold change (point). Species are
    ordered top-to-bottom by Newick leaf order.

BIC family2/family3 hyper-mutated trios and SAR family3 low-callability trios
are already excluded upstream via the `included` flag in
per_trio_germline_mutations_callable.tsv; we filter on `included == 'yes'`.
"""
import argparse
import math

import yaml
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.size': 13,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 11,
    'legend.title_fontsize': 12,
})

GERMLINE_COLOR = '#B8860B'  # dark goldenrod
SOMATIC_COLOR = '#2D7DD2'   # blue
FC_COLOR = '#444444'
SOCIAL_COLOR = '#8B2E2E'
SUBSOCIAL_COLOR = '#2F4F7F'


# --- Newick ----------------------------------------------------------------
# Minimal parser, copied from dnm_spectrum/workflow_sources/spectrum_phylo_mantel.py
# so this script does not depend on a sibling workflow's import path.

def parse_newick(s):
    s = s.strip().rstrip(';').strip()
    pos = [0]

    def _clade():
        children = []
        if pos[0] < len(s) and s[pos[0]] == '(':
            pos[0] += 1
            children.append(_clade())
            while s[pos[0]] == ',':
                pos[0] += 1
                children.append(_clade())
            if s[pos[0]] != ')':
                raise ValueError(f'Expected ) at pos {pos[0]}: {s}')
            pos[0] += 1
        label_start = pos[0]
        while pos[0] < len(s) and s[pos[0]] not in ',()':
            pos[0] += 1
        label = s[label_start:pos[0]].strip()
        return (children, label)

    return _clade()


def assign_layout(tree):
    """Walk the parsed tree and lay out a cladogram with all leaves aligned
    at x = max_height and each internal node placed at max_height - subtree
    height (so terminal branch length is identical for every leaf, regardless
    of clade depth). Returns (segments, leaf_order, max_height)."""
    leaf_counter = [0]
    leaf_order = []

    def walk(clade):
        children, label = clade
        if not children:
            y = leaf_counter[0]
            leaf_counter[0] += 1
            leaf_order.append(label)
            return {'height': 0, 'y': y, 'label': label, 'children': []}
        kids = [walk(c) for c in children]
        y = sum(k['y'] for k in kids) / len(kids)
        height = 1 + max(k['height'] for k in kids)
        return {'height': height, 'y': y, 'label': label, 'children': kids}

    root = walk(tree)
    max_height = root['height']

    def assign_x(node):
        node['x'] = max_height - node['height']
        for c in node['children']:
            assign_x(c)
    assign_x(root)

    segments = []

    def emit(node):
        if not node['children']:
            return
        ys = [c['y'] for c in node['children']]
        segments.append(((node['x'], min(ys)), (node['x'], max(ys))))
        for c in node['children']:
            segments.append(((node['x'], c['y']), (c['x'], c['y'])))
            emit(c)
    emit(root)

    return segments, leaf_order, max_height


# --- TSV readers -----------------------------------------------------------

def read_tsv(path):
    rows = []
    with open(path) as f:
        header = None
        for line in f:
            if line.startswith('#'):
                continue
            if header is None:
                header = line.rstrip('\n').split('\t')
                continue
            fields = line.rstrip('\n').split('\t')
            rows.append(dict(zip(header, fields)))
    return rows


def per_trio_overall(per_trio_path):
    """Return {species: [(n_mutations, callable_bp), ...]} aggregated per
    (species, offspring) over chromosomes, restricted to mutation_class ==
    'overall' and included == 'yes'."""
    by_trio = {}
    for r in read_tsv(per_trio_path):
        if r['mutation_class'] != 'overall':
            continue
        if r['included'] != 'yes':
            continue
        sp = r['species']
        off = r['offspring_id']
        n = int(r['n_germline_mutations'])
        bp = int(r['callable_denom_bp'])
        key = (sp, off)
        prev = by_trio.get(key, (0, 0))
        by_trio[key] = (prev[0] + n, prev[1] + bp)
    out = {}
    for (sp, _off), (n, bp) in by_trio.items():
        out.setdefault(sp, []).append((n, bp))
    return out


def per_species_somatic_overall(per_species_somatic_path):
    out = {}
    for r in read_tsv(per_species_somatic_path):
        if r['mutation_class'] != 'overall':
            continue
        out[r['species']] = (
            int(r['n_somatic_mutations']),
            int(r['callable_trinuc_count']),
        )
    return out


# --- Bootstrap -------------------------------------------------------------

def germline_rate_with_ci(trios, n_bootstrap, rng):
    """trios: list of (n_mut, callable_bp) per trio. Rate is summed over the
    list with denom = 2 * sum(callable_bp) (diploid haploid base pairs).
    Bootstrap resamples trios with replacement."""
    if not trios:
        return float('nan'), float('nan'), float('nan'), 0
    arr = np.array(trios, dtype=np.int64)
    total_n = int(arr[:, 0].sum())
    total_bp = int(arr[:, 1].sum())
    point = total_n / (2 * total_bp) if total_bp > 0 else float('nan')
    k = arr.shape[0]
    rates = np.empty(n_bootstrap, dtype=float)
    for b in range(n_bootstrap):
        idx = rng.integers(0, k, size=k)
        sub_n = int(arr[idx, 0].sum())
        sub_bp = int(arr[idx, 1].sum())
        rates[b] = sub_n / (2 * sub_bp) if sub_bp > 0 else float('nan')
    lo, hi = np.nanpercentile(rates, [2.5, 97.5])
    return point, float(lo), float(hi), k


# --- Plot ------------------------------------------------------------------

def render_pdf(rows, segments, max_depth, social_set, subsocial_set,
               pairs, output_pdf):
    """rows: list of dicts in y-axis order top->bottom. species y-position is
    the row index (matches the leaf order assigned by assign_layout)."""
    n = len(rows)
    species = [r['species'] for r in rows]
    y_pos = np.arange(n)

    fig, axes = plt.subplots(
        1, 3, figsize=(12, 4),
        gridspec_kw={'width_ratios': [1.3, 1.2, 1.0]},
        sharey=True,
    )
    axTree, axRate, axFC = axes

    # --- Phylogeny with social/subsocial tip indicator ---
    for (x1, y1), (x2, y2) in segments:
        axTree.plot([x1, x2], [y1, y2], color='black', linewidth=1.2)
    tip_x = max_depth + 0.15
    for i, sp in enumerate(species):
        if sp in social_set:
            color = SOCIAL_COLOR
        elif sp in subsocial_set:
            color = SUBSOCIAL_COLOR
        else:
            color = '#888888'
        axTree.plot(tip_x, i, marker='s', markersize=8, color=color,
                    markeredgecolor='black', markeredgewidth=0.4)
        axTree.text(tip_x + 0.18, i, sp, ha='left', va='center')
    axTree.set_xlim(-0.2, max_depth + 1.6)
    for s in axTree.spines.values():
        s.set_visible(False)
    axTree.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    # Lifestyle legend
    legend_handles = [
        plt.Line2D([0], [0], marker='s', linestyle='', color=SOCIAL_COLOR,
                   markeredgecolor='black', markeredgewidth=0.4,
                   markersize=8, label='social'),
        plt.Line2D([0], [0], marker='s', linestyle='', color=SUBSOCIAL_COLOR,
                   markeredgecolor='black', markeredgewidth=0.4,
                   markersize=8, label='subsocial'),
    ]
    axTree.legend(handles=legend_handles,
                  loc='lower center', bbox_to_anchor=(0.5, 1.02),
                  ncol=2, frameon=False)

    # --- Combined rate panel: germline (point + 95% CI) and somatic (point) ---
    g_rate = np.array([r['germline_rate'] for r in rows], dtype=float)
    g_lo = np.array([r['germline_ci_low'] for r in rows], dtype=float)
    g_hi = np.array([r['germline_ci_high'] for r in rows], dtype=float)
    xerr_lo = np.where(np.isfinite(g_rate) & np.isfinite(g_lo),
                       np.maximum(g_rate - g_lo, 0.0), 0.0)
    xerr_hi = np.where(np.isfinite(g_rate) & np.isfinite(g_hi),
                       np.maximum(g_hi - g_rate, 0.0), 0.0)
    s_rate = np.array([r['somatic_rate'] for r in rows], dtype=float)
    axRate.errorbar(g_rate, y_pos, xerr=[xerr_lo, xerr_hi], fmt='o',
                    color=GERMLINE_COLOR, ecolor=GERMLINE_COLOR,
                    capsize=3, markersize=5, linewidth=1.2, label='germline')
    axRate.plot(s_rate, y_pos, 'o', color=SOMATIC_COLOR,
                markersize=5, label='somatic')
    # Arrows from subsocial -> sister social, per series.
    sp_to_idx = {sp: i for i, sp in enumerate(species)}
    for social_sp, subsocial_sp in pairs:
        if social_sp not in sp_to_idx or subsocial_sp not in sp_to_idx:
            continue
        i_soc = sp_to_idx[social_sp]
        i_sub = sp_to_idx[subsocial_sp]
        # germline arrow
        if np.isfinite(g_rate[i_soc]) and np.isfinite(g_rate[i_sub]):
            axRate.annotate(
                '', xy=(g_rate[i_soc], i_soc),
                xytext=(g_rate[i_sub], i_sub),
                arrowprops=dict(arrowstyle='->', color=GERMLINE_COLOR,
                                lw=1.0, alpha=0.7,
                                shrinkA=4, shrinkB=4),
            )
        # somatic arrow
        if np.isfinite(s_rate[i_soc]) and np.isfinite(s_rate[i_sub]):
            axRate.annotate(
                '', xy=(s_rate[i_soc], i_soc),
                xytext=(s_rate[i_sub], i_sub),
                arrowprops=dict(arrowstyle='->', color=SOMATIC_COLOR,
                                lw=1.0, alpha=0.7,
                                shrinkA=4, shrinkB=4),
            )

    finite = np.concatenate([g_rate[np.isfinite(g_rate)],
                             s_rate[np.isfinite(s_rate)]])
    if finite.size and finite.min() > 0 and finite.max() / finite.min() > 100:
        axRate.set_xscale('log')
    axRate.set_xlabel('mutation rate per site')
    axRate.grid(axis='x', linestyle=':', alpha=0.4)
    axRate.legend(loc='lower center', bbox_to_anchor=(0.5, 1.02),
                  ncol=2, frameon=False)

    # --- Fold change panel ---
    fc = np.array([r['fold_change'] for r in rows], dtype=float)
    axFC.plot(fc, y_pos, 'o', color=FC_COLOR, markersize=5)
    axFC.axvline(1.0, color='grey', linestyle='--', linewidth=0.8)
    if np.all(np.isfinite(fc)) and np.nanmax(fc) / max(np.nanmin(fc), 1e-12) > 100:
        axFC.set_xscale('log')
    axFC.set_xlabel('fold change (somatic / germline)')
    axFC.grid(axis='x', linestyle=':', alpha=0.4)

    # Shared y-axis: species labels are drawn inside the phylogeny panel, so
    # suppress the y tick labels on the rate / fold-change panels.
    axTree.set_ylim(n - 0.4, -0.6)  # top = first leaf in Newick
    for ax in axes:
        ax.tick_params(left=False, labelleft=False)

    fig.tight_layout()
    fig.savefig(output_pdf)
    plt.close(fig)


# --- Main ------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--per_trio_tsv', required=True)
    ap.add_argument('--per_species_somatic_tsv', required=True)
    ap.add_argument('--spectrum_config', required=True)
    ap.add_argument('--output_tsv', required=True)
    ap.add_argument('--output_pdf', required=True)
    ap.add_argument('--n_bootstrap', type=int, default=10000)
    ap.add_argument('--seed', type=int, default=42)
    args = ap.parse_args()

    with open(args.spectrum_config) as f:
        cfg = yaml.safe_load(f)
    newick = cfg['species_tree_newick']
    groups = cfg.get('spectrum_groups', {})
    social_set = set(groups.get('social', {}).get('species', []))
    subsocial_set = set(groups.get('subsocial', {}).get('species', []))
    pairs = [tuple(p) for p in cfg.get('social_subsocial_pairs', [])]

    tree = parse_newick(newick)
    segments, leaf_order, max_depth = assign_layout(tree)

    germline_trios = per_trio_overall(args.per_trio_tsv)
    somatic = per_species_somatic_overall(args.per_species_somatic_tsv)

    rng = np.random.default_rng(args.seed)
    rows = []
    for sp in leaf_order:
        trios = germline_trios.get(sp, [])
        g_rate, g_lo, g_hi, n_trios = germline_rate_with_ci(
            trios, args.n_bootstrap, rng,
        )
        s_n, s_d = somatic.get(sp, (0, 0))
        s_rate = (s_n / s_d) if s_d > 0 else float('nan')
        fc = (s_rate / g_rate) if (g_rate and not math.isnan(g_rate)
                                    and g_rate > 0) else float('nan')
        rows.append({
            'species': sp,
            'germline_rate': g_rate,
            'germline_ci_low': g_lo,
            'germline_ci_high': g_hi,
            'germline_n_trios': n_trios,
            'somatic_rate': s_rate,
            'fold_change': fc,
        })

    with open(args.output_tsv, 'w') as out:
        out.write(
            'species\tphylo_order\tgermline_rate\tgermline_ci_low\t'
            'germline_ci_high\tgermline_n_trios\tsomatic_rate\tfold_change\n'
        )
        for i, r in enumerate(rows):
            out.write(
                f"{r['species']}\t{i}\t{r['germline_rate']:g}\t"
                f"{r['germline_ci_low']:g}\t{r['germline_ci_high']:g}\t"
                f"{r['germline_n_trios']}\t{r['somatic_rate']:g}\t"
                f"{r['fold_change']:g}\n"
            )

    render_pdf(rows, segments, max_depth, social_set, subsocial_set,
               pairs, args.output_pdf)


if __name__ == '__main__':
    main()
