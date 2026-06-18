#!/usr/bin/env python3
"""
compare_groups_test.py
Test whether the distribution of mutations across the 7 collapsed categories
(C>A, C>G, C>T non-CpG, C>T CpG, T>A, T>C, T>G) differs between the social
and subsocial species groups, separately for germline and somatic mutations.

Global test:    chi-square on 7x2 contingency (scipy.stats.chi2_contingency)
Per-category:   two-sided Fisher's exact on 2x2 (category vs rest)
Correction:     Benjamini-Hochberg (BH) across the 7 per-category p-values
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
CAT7_COLORS = ['#03BCEE', '#010101', '#E32926', '#E32926',
               '#CAC9C9', '#A1CE63', '#EBC6C4']


def read_germline(path):
    data = {}
    with open(path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
            row = dict(zip(header, fields))
            data[(row['mutation_type'], row['trinuc_context'])] = int(row['count'])
    return data


def read_somatic(path):
    data = {}
    with open(path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
            row = dict(zip(header, fields))
            data[(row['mutation_type'], row['context'])] = int(row['count'])
    return data


def merge_counts(dicts):
    merged = {}
    for d in dicts:
        for k, v in d.items():
            merged[k] = merged.get(k, 0) + v
    return merged


def counts_to_7category(counts):
    """Collapse a {(mut, trinuc)->count} dict to a length-7 numpy array.
    Split C>T into non-CpG (3' != G) and CpG (3' == G).
    """
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
                        vec[3] += c  # CpG
                    else:
                        vec[2] += c  # non-CpG
                elif mut == 'T>A':
                    vec[4] += c
                elif mut == 'T>C':
                    vec[5] += c
                elif mut == 'T>G':
                    vec[6] += c
    return vec


def bh_adjust(pvals):
    """Benjamini-Hochberg FDR correction. Returns array of adjusted p-values."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    adj = ranked * n / (np.arange(n) + 1)
    # enforce monotonicity (running min from the end)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.empty_like(adj)
    out[order] = adj
    return out


def multinomial_bootstrap(counts, n_boot, rng):
    total = counts.sum()
    if total == 0:
        return np.zeros((n_boot, len(counts)))
    p = counts / total
    return rng.multinomial(int(total), p, size=n_boot).astype(float) / total


def compare_two_groups(social_counts, subsocial_counts):
    """Return dict of global chi-square result and per-category rows."""
    total_soc = social_counts.sum()
    total_sub = subsocial_counts.sum()

    # Global chi-square on 7x2 contingency (rows=category, cols=group)
    table = np.vstack([social_counts, subsocial_counts]).T  # shape (7, 2)
    chi2, p_global, dof, _ = stats.chi2_contingency(table)

    per_cat = []
    pvals = []
    for i in range(7):
        x_soc = int(social_counts[i])
        x_sub = int(subsocial_counts[i])
        tbl = [[x_soc, int(total_soc) - x_soc],
               [x_sub, int(total_sub) - x_sub]]
        odds, p = stats.fisher_exact(tbl, alternative='two-sided')
        per_cat.append({
            'category': CAT7_SHORT[i],
            'social_count': x_soc,
            'subsocial_count': x_sub,
            'social_prop': x_soc / total_soc if total_soc else 0,
            'subsocial_prop': x_sub / total_sub if total_sub else 0,
            'odds_ratio': odds,
            'p_value': p,
        })
        pvals.append(p)

    padj = bh_adjust(pvals)
    for i, row in enumerate(per_cat):
        row['p_adjusted'] = padj[i]
        row['significant'] = 'yes' if padj[i] < 0.05 else 'no'

    return {
        'chi2': chi2,
        'p_global': p_global,
        'dof': dof,
        'total_social': int(total_soc),
        'total_subsocial': int(total_sub),
        'rows': per_cat,
    }


def sig_stars(p):
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    return 'ns'


def write_tsv(path, germline_result, somatic_result):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write('# Global chi-square tests (7x2 contingency)\n')
        f.write('comparison\tchi2\tdof\tp_global\tn_social\tn_subsocial\n')
        for name, r in [('germline', germline_result), ('somatic', somatic_result)]:
            f.write(f"{name}\t{r['chi2']:.4f}\t{r['dof']}\t{r['p_global']:.4g}\t"
                    f"{r['total_social']}\t{r['total_subsocial']}\n")
        f.write('\n# Per-category Fisher exact tests with BH-adjusted p-values\n')
        f.write('comparison\tcategory\tsocial_count\tsubsocial_count\t'
                'social_prop\tsubsocial_prop\todds_ratio\tp_value\tp_adjusted\tsignificant\n')
        for name, r in [('germline', germline_result), ('somatic', somatic_result)]:
            for row in r['rows']:
                f.write(f"{name}\t{row['category']}\t{row['social_count']}\t{row['subsocial_count']}\t"
                        f"{row['social_prop']:.6f}\t{row['subsocial_prop']:.6f}\t"
                        f"{row['odds_ratio']:.4g}\t{row['p_value']:.4g}\t"
                        f"{row['p_adjusted']:.4g}\t{row['significant']}\n")


def plot_result(social_counts, subsocial_counts, result, ax, title, rng, n_boot=1000):
    soc_draws = multinomial_bootstrap(social_counts, n_boot, rng)
    sub_draws = multinomial_bootstrap(subsocial_counts, n_boot, rng)

    soc_mean = soc_draws.mean(axis=0)
    soc_lo = np.quantile(soc_draws, 0.025, axis=0)
    soc_hi = np.quantile(soc_draws, 0.975, axis=0)
    sub_mean = sub_draws.mean(axis=0)
    sub_lo = np.quantile(sub_draws, 0.025, axis=0)
    sub_hi = np.quantile(sub_draws, 0.975, axis=0)

    x = np.arange(7)
    w = 0.38
    # social bars (left) — solid canonical color per category
    soc_err = np.vstack([soc_mean - soc_lo, soc_hi - soc_mean])
    ax.bar(x - w / 2, soc_mean, width=w,
           color=CAT7_COLORS, edgecolor='black', linewidth=0.5,
           yerr=soc_err, error_kw=dict(ecolor='black', lw=0.8, capsize=3),
           label='social')
    # subsocial bars (right) — hatched
    sub_err = np.vstack([sub_mean - sub_lo, sub_hi - sub_mean])
    ax.bar(x + w / 2, sub_mean, width=w,
           color=CAT7_COLORS, edgecolor='black', linewidth=0.5,
           hatch='///', alpha=0.55,
           yerr=sub_err, error_kw=dict(ecolor='black', lw=0.8, capsize=3),
           label='subsocial')

    # Significance annotation
    ymax = max(soc_hi.max(), sub_hi.max()) * 1.25
    for i, row in enumerate(result['rows']):
        stars = sig_stars(row['p_adjusted'])
        top = max(soc_hi[i], sub_hi[i])
        ax.text(i, top + ymax * 0.02, stars, ha='center', va='bottom',
                fontsize=10, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(CAT7_LABELS, fontsize=10)
    ax.set_ylim(0, ymax)
    ax.set_ylabel('Proportion')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    global_p = result['p_global']
    ax.set_title(
        f'{title}  —  social (n={result["total_social"]}) vs subsocial (n={result["total_subsocial"]})  —  '
        f'global \u03c7\u00b2 p={global_p:.3g}',
        fontsize=11,
    )


def plot_tests(germline_result, germline_counts, somatic_result, somatic_counts, out_pdf, seed=42):
    rng = np.random.default_rng(seed)

    fig, axes = plt.subplots(2, 1, figsize=(11, 8))
    plot_result(germline_counts['social'], germline_counts['subsocial'],
                germline_result, axes[0], 'Germline (autosome)', rng)
    plot_result(somatic_counts['social'], somatic_counts['subsocial'],
                somatic_result, axes[1], 'Somatic (merged)', rng)

    # Legend
    soc_handle = plt.Rectangle((0, 0), 1, 1, facecolor='grey',
                                edgecolor='black', linewidth=0.5)
    sub_handle = plt.Rectangle((0, 0), 1, 1, facecolor='grey',
                                edgecolor='black', linewidth=0.5,
                                hatch='///', alpha=0.55)
    fig.legend([soc_handle, sub_handle], ['Social', 'Subsocial'],
               loc='upper right', fontsize=10, frameon=False)

    fig.suptitle(
        'Social vs Subsocial — 7-category mutation spectrum comparison\n'
        '(95% bootstrap CI; BH-adjusted Fisher exact: *p<0.05, **p<0.01, ***p<0.001)',
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 0.94, 0.94])
    os.makedirs(os.path.dirname(out_pdf), exist_ok=True)
    fig.savefig(out_pdf, dpi=200, bbox_inches='tight')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description='Social vs subsocial 7-category test')
    parser.add_argument('--germline_social', required=True)
    parser.add_argument('--germline_subsocial', required=True)
    parser.add_argument('--somatic_social', nargs='+', required=True)
    parser.add_argument('--somatic_subsocial', nargs='+', required=True)
    parser.add_argument('--output_tsv', required=True)
    parser.add_argument('--output_pdf', required=True)
    args = parser.parse_args()

    # Germline: group-level spectrum TSVs
    g_soc = counts_to_7category(read_germline(args.germline_social))
    g_sub = counts_to_7category(read_germline(args.germline_subsocial))

    # Somatic: merge species TSVs per group
    s_soc = counts_to_7category(merge_counts([read_somatic(p) for p in args.somatic_social]))
    s_sub = counts_to_7category(merge_counts([read_somatic(p) for p in args.somatic_subsocial]))

    germline_result = compare_two_groups(g_soc, g_sub)
    somatic_result = compare_two_groups(s_soc, s_sub)

    write_tsv(args.output_tsv, germline_result, somatic_result)
    plot_tests(germline_result, {'social': g_soc, 'subsocial': g_sub},
               somatic_result, {'social': s_soc, 'subsocial': s_sub},
               args.output_pdf)

    print(f'TSV -> {args.output_tsv}')
    print(f'PDF -> {args.output_pdf}')
    print(f'Germline global chi2 p={germline_result["p_global"]:.4g}')
    print(f'Somatic  global chi2 p={somatic_result["p_global"]:.4g}')


if __name__ == '__main__':
    main()
