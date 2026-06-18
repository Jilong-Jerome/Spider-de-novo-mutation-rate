#!/usr/bin/env python3
"""
pairwise_somatic_test.py
For each social-subsocial species pair, test whether the somatic mutation
spectrum differs between the two species across the 7 collapsed categories,
and check whether the direction of difference is consistent with the merged
social-vs-subsocial somatic comparison.

Per pair, per category: two-sided Fisher's exact on 2x2
    [[social_i,    n_social    - social_i],
     [subsocial_i, n_subsocial - subsocial_i]]
with BH adjustment across 7 categories within each pair.

Across pairs, per category: Cochran-Mantel-Haenszel test on the 3 strata
(one 2x2 per pair), reported with common odds ratio.
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


def cmh_test(tables):
    """Cochran-Mantel-Haenszel test across K 2x2 tables.
    Each table [[a,b],[c,d]] with row1=social, row2=subsocial,
    col1=in-category, col2=rest.
    Returns (chi2, p, or_mh).
    """
    tables = [np.asarray(t, dtype=float) for t in tables]
    A = sum(t[0, 0] for t in tables)
    EA = 0.0
    V = 0.0
    R = 0.0
    S = 0.0
    for t in tables:
        a, b = t[0]
        c, d = t[1]
        n1 = a + b
        n2 = c + d
        m1 = a + c
        m2 = b + d
        N = n1 + n2
        EA += n1 * m1 / N
        if N > 1:
            V += (n1 * n2 * m1 * m2) / (N * N * (N - 1))
        R += a * d / N
        S += b * c / N
    if V <= 0:
        return (np.nan, np.nan, np.nan)
    chi2 = (abs(A - EA) - 0.5) ** 2 / V
    p = 1 - stats.chi2.cdf(chi2, 1)
    or_mh = R / S if S > 0 else np.nan
    return (chi2, p, or_mh)


def per_pair_per_category(soc_counts, sub_counts):
    """Return list of 7 row-dicts with Fisher's exact + OR per category."""
    n_soc = soc_counts.sum()
    n_sub = sub_counts.sum()
    rows = []
    pvals = []
    for i in range(7):
        a = int(soc_counts[i]); b = int(n_soc) - a
        c = int(sub_counts[i]); d = int(n_sub) - c
        odds, p = stats.fisher_exact([[a, b], [c, d]], alternative='two-sided')
        rows.append({
            'category': CAT7_SHORT[i],
            'social_count': a, 'subsocial_count': c,
            'social_prop': a / n_soc if n_soc else 0,
            'subsocial_prop': c / n_sub if n_sub else 0,
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


def write_tsv(path, pair_names, pair_results, cmh_rows, consistency_rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write('# Per-pair per-category Fisher exact (BH within each pair)\n')
        f.write('pair\tcategory\tsocial_count\tsubsocial_count\t'
                'social_prop\tsubsocial_prop\todds_ratio\tlog2_or\t'
                'p_value\tp_adjusted\tsignificant\n')
        for pname, rows in zip(pair_names, pair_results):
            for r in rows:
                f.write(f"{pname}\t{r['category']}\t{r['social_count']}\t{r['subsocial_count']}\t"
                        f"{r['social_prop']:.6f}\t{r['subsocial_prop']:.6f}\t"
                        f"{r['odds_ratio']:.4g}\t{r['log2_or']:.4f}\t"
                        f"{r['p_value']:.4g}\t{r['p_adjusted']:.4g}\t{r['significant']}\n")

        f.write('\n# Per-category Cochran-Mantel-Haenszel test across pairs (BH across 7)\n')
        f.write('category\tcmh_chi2\tcmh_p\tcmh_p_adjusted\tor_mh\tsignificant\n')
        for r in cmh_rows:
            f.write(f"{r['category']}\t{r['cmh_chi2']:.4f}\t{r['cmh_p']:.4g}\t"
                    f"{r['cmh_p_adjusted']:.4g}\t{r['or_mh']:.4g}\t{r['significant']}\n")

        f.write('\n# Direction consistency across pairs (sign of log2 OR)\n')
        f.write('category\tpairs_pos\tpairs_neg\tconsistent\n')
        for r in consistency_rows:
            f.write(f"{r['category']}\t{r['pairs_pos']}\t{r['pairs_neg']}\t{r['consistent']}\n")


def plot_panels(pair_names, pair_soc, pair_sub, pair_results, cmh_rows, out_pdf):
    n_pairs = len(pair_names)
    fig, axes = plt.subplots(1, 7, figsize=(16, 4.2), sharey=True)

    # shared y-range across categories
    all_lor = np.array([[r['log2_or'] for r in rows] for rows in pair_results])
    span = max(1.0, np.nanmax(np.abs(all_lor))) * 1.45

    x = np.arange(n_pairs)
    for i, ax in enumerate(axes):
        vals = [rows[i]['log2_or'] for rows in pair_results]
        stars = [sig_stars(rows[i]['p_adjusted']) for rows in pair_results]
        bars = ax.bar(x, vals, color=CAT7_COLORS[i],
                      edgecolor='black', linewidth=0.5)
        ax.axhline(0, color='black', lw=0.5)
        for xi, (v, s) in enumerate(zip(vals, stars)):
            top = v if v >= 0 else 0
            ax.text(xi, top + span * 0.04, s, ha='center', va='bottom',
                    fontsize=9, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(pair_names, rotation=30, ha='right', fontsize=8)
        cmh_p = cmh_rows[i]['cmh_p']
        cmh_padj = cmh_rows[i]['cmh_p_adjusted']
        ax.set_title(f"{CAT7_LABELS[i]}\nCMH p={cmh_p:.3g}\n(adj {cmh_padj:.3g})",
                     fontsize=9)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    axes[0].set_ylabel('log2 OR (social / subsocial)\nsomatic, per category')
    axes[0].set_ylim(-span, span)

    fig.suptitle(
        'Pairwise somatic spectrum comparison (social vs subsocial)\n'
        f'pairs: {", ".join(pair_names)}  —  BH within each pair; CMH BH across categories',
        fontsize=11,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    os.makedirs(os.path.dirname(out_pdf), exist_ok=True)
    fig.savefig(out_pdf, dpi=200, bbox_inches='tight')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description='Per-pair somatic social vs subsocial comparison')
    parser.add_argument('--pair', action='append', nargs=3,
                        metavar=('PAIR_NAME', 'SOCIAL_TSV', 'SUBSOCIAL_TSV'),
                        required=True,
                        help='Repeat per pair: name, social TSV, subsocial TSV')
    parser.add_argument('--output_tsv', required=True)
    parser.add_argument('--output_pdf', required=True)
    args = parser.parse_args()

    pair_names = [p[0] for p in args.pair]
    pair_soc = [counts_to_7category(read_somatic(p[1])) for p in args.pair]
    pair_sub = [counts_to_7category(read_somatic(p[2])) for p in args.pair]

    pair_results = [per_pair_per_category(s, u) for s, u in zip(pair_soc, pair_sub)]

    cmh_rows = []
    cmh_pvals = []
    for i in range(7):
        tables = []
        for s, u in zip(pair_soc, pair_sub):
            a = int(s[i]); b = int(s.sum()) - a
            c = int(u[i]); d = int(u.sum()) - c
            tables.append([[a, b], [c, d]])
        chi2, p, or_mh = cmh_test(tables)
        cmh_rows.append({
            'category': CAT7_SHORT[i],
            'cmh_chi2': chi2, 'cmh_p': p, 'or_mh': or_mh,
        })
        cmh_pvals.append(p)
    cmh_padj = bh_adjust(cmh_pvals)
    for i, r in enumerate(cmh_rows):
        r['cmh_p_adjusted'] = cmh_padj[i]
        r['significant'] = 'yes' if cmh_padj[i] < 0.05 else 'no'

    consistency_rows = []
    for i in range(7):
        signs = [np.sign(rows[i]['log2_or']) for rows in pair_results]
        pos = sum(1 for s in signs if s > 0)
        neg = sum(1 for s in signs if s < 0)
        consistency_rows.append({
            'category': CAT7_SHORT[i],
            'pairs_pos': pos,
            'pairs_neg': neg,
            'consistent': 'yes' if (pos == len(signs) or neg == len(signs)) else 'no',
        })

    write_tsv(args.output_tsv, pair_names, pair_results, cmh_rows, consistency_rows)
    plot_panels(pair_names, pair_soc, pair_sub, pair_results, cmh_rows, args.output_pdf)

    print(f'TSV -> {args.output_tsv}')
    print(f'PDF -> {args.output_pdf}')
    for i, r in enumerate(cmh_rows):
        print(f"{CAT7_SHORT[i]}: CMH p={r['cmh_p']:.3g}, "
              f"adj p={r['cmh_p_adjusted']:.3g}, OR_MH={r['or_mh']:.3g}")


if __name__ == '__main__':
    main()
