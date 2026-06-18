#!/usr/bin/env python3
"""
baseline_somatic_test.py
Use the merged subsocial somatic spectrum as the ancestral baseline and
test each species' somatic spectrum against it, category by category,
to identify mutation types that shift consistently in the social species.

For each species SP and each of the 7 collapsed categories i:
    table = [[sp_i,       sp_total   - sp_i      ],
             [baseline_i, baseline_total - baseline_i]]
    two-sided Fisher's exact -> p, OR; BH across 7 within each species.

Baseline definition:
    baseline(SP) = merge(subsocial species \\ {SP}) if SP is subsocial,
                   else merge(all subsocial species).
(Leave-one-out for subsocial species avoids self-comparison.)

Consistency summary:
    For each category, count how many social species show BH-significant
    shifts in the same direction (log2 OR sign).
"""
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

MUT_TYPES = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
BASES = ['A', 'C', 'G', 'T']
CAT7_LABELS = ['C>A', 'C>G', 'C>T\n(non-CpG)', 'C>T\n(CpG)', 'T>A', 'T>C', 'T>G']
CAT7_SHORT = ['C>A', 'C>G', 'C>T_nonCpG', 'C>T_CpG', 'T>A', 'T>C', 'T>G']


def read_somatic(path):
    data = {}
    with open(path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
            row = dict(zip(header, fields))
            data[(row['mutation_type'], row['context'])] = int(row['count'])
    return data


def counts_to_7category(counts):
    vec = np.zeros(7, dtype=int)
    for mut in MUT_TYPES:
        for five in BASES:
            for three in BASES:
                trinuc = f'{five}{mut[0]}{three}'
                c = counts.get((mut, trinuc), 0)
                if mut == 'C>A':
                    vec[0] += c
                elif mut == 'C>G':
                    vec[1] += c
                elif mut == 'C>T':
                    if three == 'G':
                        vec[3] += c
                    else:
                        vec[2] += c
                elif mut == 'T>A':
                    vec[4] += c
                elif mut == 'T>C':
                    vec[5] += c
                elif mut == 'T>G':
                    vec[6] += c
    return vec


def bh_adjust(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.empty_like(adj)
    out[order] = adj
    return out


def log2_or(a, b, c, d, eps=0.5):
    return np.log2(((a + eps) * (d + eps)) / ((b + eps) * (c + eps)))


def parse_pair_args(pairs):
    """Parse [NAME, PATH] list into OrderedDict-preserving name-list + dict."""
    names = [p[0] for p in pairs]
    paths = {p[0]: p[1] for p in pairs}
    return names, paths


def test_species_vs_baseline(sp_vec, base_vec):
    n_sp = sp_vec.sum()
    n_base = base_vec.sum()
    rows, pvals = [], []
    for i in range(7):
        a = int(sp_vec[i]); b = int(n_sp) - a
        c = int(base_vec[i]); d = int(n_base) - c
        odds, p = stats.fisher_exact([[a, b], [c, d]], alternative='two-sided')
        rows.append({
            'category': CAT7_SHORT[i],
            'sp_count': a, 'baseline_count': c,
            'sp_prop': a / n_sp if n_sp else 0,
            'baseline_prop': c / n_base if n_base else 0,
            'odds_ratio': odds,
            'log2_or': log2_or(a, b, c, d),
            'p_value': p,
        })
        pvals.append(p)
    padj = bh_adjust(pvals)
    for i, r in enumerate(rows):
        r['p_adjusted'] = padj[i]
        r['significant'] = 'yes' if padj[i] < 0.05 else 'no'
    return rows


def sig_stars(p):
    if np.isnan(p):
        return ''
    if p < 0.001:
        return '***'
    if p < 0.01:
        return '**'
    if p < 0.05:
        return '*'
    return 'ns'


def write_tsv(path, baseline_total, per_species_rows,
              social_names, subsocial_names, consistency_rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write(f'# Baseline: merged subsocial somatic spectrum '
                f'(leave-one-out for each subsocial species)\n')
        f.write(f'# Social species: {",".join(social_names)}\n')
        f.write(f'# Subsocial species: {",".join(subsocial_names)}\n')
        f.write(f'# Total subsocial baseline count (full): {baseline_total}\n\n')

        f.write('# Per-species per-category Fisher exact vs subsocial baseline\n')
        f.write('species\tgroup\tcategory\tsp_count\tbaseline_count\t'
                'sp_prop\tbaseline_prop\todds_ratio\tlog2_or\t'
                'p_value\tp_adjusted\tsignificant\n')
        for sp, group, rows in per_species_rows:
            for r in rows:
                f.write(f"{sp}\t{group}\t{r['category']}\t{r['sp_count']}\t{r['baseline_count']}\t"
                        f"{r['sp_prop']:.6f}\t{r['baseline_prop']:.6f}\t"
                        f"{r['odds_ratio']:.4g}\t{r['log2_or']:.4f}\t"
                        f"{r['p_value']:.4g}\t{r['p_adjusted']:.4g}\t{r['significant']}\n")

        f.write('\n# Consistency of significant shifts among SOCIAL species '
                '(BH-significant within species, sign of log2 OR)\n')
        f.write('category\tn_social_sig_up\tn_social_sig_down\tn_social_total\t'
                'consistent_direction\n')
        for r in consistency_rows:
            f.write(f"{r['category']}\t{r['n_up']}\t{r['n_down']}\t{r['n_total']}\t"
                    f"{r['consistent']}\n")


def plot_heatmap(social_names, subsocial_names,
                 per_species_rows, consistency_rows, out_pdf):
    species_order = list(social_names) + list(subsocial_names)
    groups = {sp: 'social' for sp in social_names}
    groups.update({sp: 'subsocial' for sp in subsocial_names})

    rows_by_sp = {sp: rows for sp, grp, rows in per_species_rows}

    n_sp = len(species_order)
    mat = np.zeros((n_sp, 7))
    stars = [[''] * 7 for _ in range(n_sp)]
    for r, sp in enumerate(species_order):
        rows = rows_by_sp[sp]
        for c, row in enumerate(rows):
            mat[r, c] = row['log2_or']
            if row['significant'] == 'yes':
                stars[r][c] = sig_stars(row['p_adjusted'])

    vmax = max(1.0, np.nanmax(np.abs(mat))) * 1.0

    fig, (ax_hm, ax_bar) = plt.subplots(
        2, 1, figsize=(10, 1.0 + 0.55 * n_sp + 2.3),
        gridspec_kw={'height_ratios': [max(3, n_sp), 2.2], 'hspace': 0.35},
    )

    # Heatmap
    im = ax_hm.imshow(mat, cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                      aspect='auto')
    ax_hm.set_xticks(range(7))
    ax_hm.set_xticklabels(CAT7_LABELS, fontsize=10)
    ax_hm.set_yticks(range(n_sp))
    ax_hm.set_yticklabels(
        [f'{sp} ({groups[sp]})' for sp in species_order], fontsize=10)
    for r in range(n_sp):
        for c in range(7):
            if stars[r][c]:
                ax_hm.text(c, r, stars[r][c], ha='center', va='center',
                           fontsize=9, color='black', fontweight='bold')
    # separator between social and subsocial
    if social_names and subsocial_names:
        ax_hm.axhline(len(social_names) - 0.5, color='black', lw=1.2)
    cbar = fig.colorbar(im, ax=ax_hm, shrink=0.75)
    cbar.set_label('log2 OR vs subsocial baseline')
    ax_hm.set_title('Per-species somatic shift relative to merged subsocial baseline\n'
                    '(star = BH-significant within species)', fontsize=11)

    # Consistency bar: number of BH-significant social species per category, colored by direction
    x = np.arange(7)
    up = np.array([r['n_up'] for r in consistency_rows])
    down = np.array([r['n_down'] for r in consistency_rows])
    ax_bar.bar(x, up, color='#B2182B', edgecolor='black', linewidth=0.5,
               label='social: sig. higher than baseline')
    ax_bar.bar(x, -down, color='#2166AC', edgecolor='black', linewidth=0.5,
               label='social: sig. lower than baseline')
    ax_bar.axhline(0, color='black', lw=0.5)
    max_n = max(len(social_names), 1)
    ax_bar.set_ylim(-max_n - 0.5, max_n + 0.5)
    ax_bar.set_yticks(range(-max_n, max_n + 1))
    ax_bar.set_yticklabels([str(abs(v)) for v in range(-max_n, max_n + 1)])
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels(CAT7_LABELS, fontsize=10)
    ax_bar.set_ylabel('# social species\n(significant)')
    ax_bar.set_title(f'Consistency of shifts across {max_n} social species '
                     f'({", ".join(social_names)})', fontsize=11)
    ax_bar.legend(fontsize=8, frameon=False, loc='upper right')
    ax_bar.spines['top'].set_visible(False)
    ax_bar.spines['right'].set_visible(False)

    fig.tight_layout()
    os.makedirs(os.path.dirname(out_pdf), exist_ok=True)
    fig.savefig(out_pdf, dpi=200, bbox_inches='tight')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description='Test each species somatic spectrum vs subsocial baseline')
    parser.add_argument('--social', action='append', nargs=2,
                        metavar=('NAME', 'PATH'), required=True,
                        help='Repeat per social species')
    parser.add_argument('--subsocial', action='append', nargs=2,
                        metavar=('NAME', 'PATH'), required=True,
                        help='Repeat per subsocial species (contributes to baseline)')
    parser.add_argument('--output_tsv', required=True)
    parser.add_argument('--output_pdf', required=True)
    args = parser.parse_args()

    social_names, social_paths = parse_pair_args(args.social)
    subsocial_names, subsocial_paths = parse_pair_args(args.subsocial)

    social_vecs = {sp: counts_to_7category(read_somatic(p))
                   for sp, p in social_paths.items()}
    subsocial_vecs = {sp: counts_to_7category(read_somatic(p))
                      for sp, p in subsocial_paths.items()}

    full_baseline = np.sum(list(subsocial_vecs.values()), axis=0)

    per_species_rows = []
    for sp in social_names:
        per_species_rows.append(
            (sp, 'social', test_species_vs_baseline(social_vecs[sp], full_baseline))
        )
    for sp in subsocial_names:
        # leave-one-out baseline
        loo_vecs = [v for name, v in subsocial_vecs.items() if name != sp]
        loo_baseline = np.sum(loo_vecs, axis=0)
        per_species_rows.append(
            (sp, 'subsocial', test_species_vs_baseline(subsocial_vecs[sp], loo_baseline))
        )

    # Consistency among social species
    social_rows_by_sp = {sp: rows for sp, g, rows in per_species_rows if g == 'social'}
    consistency_rows = []
    for i in range(7):
        n_up = 0; n_down = 0
        for sp in social_names:
            r = social_rows_by_sp[sp][i]
            if r['significant'] == 'yes':
                if r['log2_or'] > 0:
                    n_up += 1
                elif r['log2_or'] < 0:
                    n_down += 1
        n_sig = n_up + n_down
        consistent = 'all' if n_sig == len(social_names) and (n_up == 0 or n_down == 0) else (
            'partial' if n_sig > 0 and (n_up == 0 or n_down == 0) else 'mixed_or_none'
        )
        consistency_rows.append({
            'category': CAT7_SHORT[i],
            'n_up': n_up, 'n_down': n_down,
            'n_total': len(social_names),
            'consistent': consistent,
        })

    write_tsv(args.output_tsv, int(full_baseline.sum()),
              per_species_rows, social_names, subsocial_names, consistency_rows)
    plot_heatmap(social_names, subsocial_names,
                 per_species_rows, consistency_rows, args.output_pdf)

    print(f'TSV -> {args.output_tsv}')
    print(f'PDF -> {args.output_pdf}')
    for r in consistency_rows:
        print(f"{r['category']}: up={r['n_up']}, down={r['n_down']}, "
              f"({r['consistent']})")


if __name__ == '__main__':
    main()
