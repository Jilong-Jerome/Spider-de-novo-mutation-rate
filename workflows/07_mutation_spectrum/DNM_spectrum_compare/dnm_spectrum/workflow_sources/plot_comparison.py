#!/usr/bin/env python3
"""
plot_comparison.py
Germline (group, autosome) vs merged-somatic SBS96 comparison.
Top: 6 panels of 96 SBS96 categories, germline vs somatic side-by-side.
Bottom: collapsed 6-mutation-type comparison.
All bars have bootstrapped 95% CIs.
"""
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats

MUT_TYPES = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
MUT_COLORS = {
    'C>A': '#03BCEE',
    'C>G': '#010101',
    'C>T': '#E32926',
    'T>A': '#CAC9C9',
    'T>C': '#A1CE63',
    'T>G': '#EBC6C4',
}
BASES = ['A', 'C', 'G', 'T']


def read_germline(path):
    data = {}
    with open(path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
            row = dict(zip(header, fields))
            key = (row['mutation_type'], row['trinuc_context'])
            data[key] = int(row['count'])
    return data


def read_somatic(path):
    data = {}
    with open(path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
            row = dict(zip(header, fields))
            key = (row['mutation_type'], row['context'])
            data[key] = int(row['count'])
    return data


def build_category_list():
    cats = []
    for mut in MUT_TYPES:
        for five in BASES:
            for three in BASES:
                cats.append((mut, f'{five}{mut[0]}{three}'))
    return cats


def bootstrap_draws(counts_vec, n_boot, rng):
    """
    Multinomial bootstrap. Returns draws array of shape (n_boot, 96) with
    proportions (count / total) for each category in each draw.
    """
    total = counts_vec.sum()
    if total == 0:
        return np.zeros((n_boot, len(counts_vec)))
    p = counts_vec / total
    draws = rng.multinomial(int(total), p, size=n_boot).astype(float) / total
    return draws


def summarize(draws):
    """Return (mean, 2.5% quantile, 97.5% quantile) along axis 0."""
    return (draws.mean(axis=0),
            np.quantile(draws, 0.025, axis=0),
            np.quantile(draws, 0.975, axis=0))


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


CAT7_SHORT = ['C>A', 'C>G', 'C>T_nonCpG', 'C>T_CpG', 'T>A', 'T>C', 'T>G']


def collapse_counts_7(vec96):
    """Collapse an integer/float length-96 vector (layout: mut x five x three)
    into 7 categories, splitting C>T by 3' G (CpG) vs not."""
    d = np.asarray(vec96, dtype=float).reshape(6, 4, 4)
    ct = MUT_TYPES.index('C>T')
    g = BASES.index('G')
    out = np.zeros(7)
    out[0] = d[MUT_TYPES.index('C>A')].sum()
    out[1] = d[MUT_TYPES.index('C>G')].sum()
    ct_cpg = d[ct, :, g].sum()
    ct_all = d[ct].sum()
    out[2] = ct_all - ct_cpg
    out[3] = ct_cpg
    out[4] = d[MUT_TYPES.index('T>A')].sum()
    out[5] = d[MUT_TYPES.index('T>C')].sum()
    out[6] = d[MUT_TYPES.index('T>G')].sum()
    return out


def fisher_per_category(g_counts, s_counts, labels, mut_ctx=None):
    """Run two-sided Fisher's exact per category on 2x2
    [[g_i, g_total - g_i], [s_i, s_total - s_i]] tables.
    mut_ctx: optional list of (mutation_type, trinuc_context) parallel to labels;
    when None, both fields are left blank in the output rows."""
    g_counts = np.asarray(g_counts, dtype=float)
    s_counts = np.asarray(s_counts, dtype=float)
    g_total = g_counts.sum()
    s_total = s_counts.sum()
    rows = []
    pvals = []
    eps = 0.5
    for i, cat in enumerate(labels):
        a = int(round(g_counts[i])); b = int(round(g_total)) - a
        c = int(round(s_counts[i])); d = int(round(s_total)) - c
        odds, p = stats.fisher_exact([[a, b], [c, d]], alternative='two-sided')
        log2_or = float(np.log2(((a + eps) * (d + eps)) / ((b + eps) * (c + eps))))
        mut, trinuc = ('', '') if mut_ctx is None else mut_ctx[i]
        rows.append({
            'category': cat,
            'mutation_type': mut,
            'trinuc_context': trinuc,
            'germline_count': a,
            'somatic_count': c,
            'germline_prop': a / g_total if g_total else 0.0,
            'somatic_prop': c / s_total if s_total else 0.0,
            'odds_ratio': float(odds) if np.isfinite(odds) else float('nan'),
            'log2_or': log2_or,
            'p_value': float(p),
        })
        pvals.append(p)
    padj = bh_adjust(pvals)
    for i, r in enumerate(rows):
        r['p_adjusted'] = float(padj[i])
        r['significant'] = 'yes' if padj[i] < 0.05 else 'no'
    return rows


def write_tests_tsv(path, group_name, g_total, s_total, n_somatic_files,
                    rows96, rows7):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write(f'# group: {group_name}\n')
        f.write(f'# germline_total={g_total}\n')
        f.write(f'# somatic_total={s_total}\n')
        f.write(f'# n_somatic_files_merged={n_somatic_files}\n')
        f.write(f'# test: two-sided Fishers exact per category; BH-adjusted within each level (96 and 7)\n')
        cols = ['level', 'category', 'mutation_type', 'trinuc_context',
                'germline_count', 'somatic_count',
                'germline_prop', 'somatic_prop',
                'odds_ratio', 'log2_or',
                'p_value', 'p_adjusted', 'significant']
        f.write('\t'.join(cols) + '\n')
        for level, rows in (('96', rows96), ('7', rows7)):
            for r in rows:
                f.write(
                    f"{level}\t{r['category']}\t{r['mutation_type']}\t{r['trinuc_context']}\t"
                    f"{r['germline_count']}\t{r['somatic_count']}\t"
                    f"{r['germline_prop']:.6f}\t{r['somatic_prop']:.6f}\t"
                    f"{r['odds_ratio']:.4g}\t{r['log2_or']:.4f}\t"
                    f"{r['p_value']:.4g}\t{r['p_adjusted']:.4g}\t{r['significant']}\n"
                )


def main():
    parser = argparse.ArgumentParser(description='Germline vs merged-somatic SBS96 comparison with bootstrap CI')
    parser.add_argument('--germline', required=True)
    parser.add_argument('--somatic', nargs='+', required=True)
    parser.add_argument('--group_name', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--output_tsv', required=True)
    parser.add_argument('--n_bootstrap', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)

    germline_counts = read_germline(args.germline)
    merged_somatic = {}
    for path in args.somatic:
        for key, val in read_somatic(path).items():
            merged_somatic[key] = merged_somatic.get(key, 0) + val

    cats = build_category_list()
    g_vec = np.array([germline_counts.get(c, 0) for c in cats], dtype=float)
    s_vec = np.array([merged_somatic.get(c, 0) for c in cats], dtype=float)

    g_total = int(g_vec.sum())
    s_total = int(s_vec.sum())

    # --- Per-category Fisher's exact tests (germline vs merged somatic) ---
    cat96_labels = [f'{c[1][0]}[{c[0]}]{c[1][2]}' for c in cats]
    rows96 = fisher_per_category(g_vec, s_vec, cat96_labels, mut_ctx=cats)
    g_vec_7 = collapse_counts_7(g_vec)
    s_vec_7 = collapse_counts_7(s_vec)
    rows7 = fisher_per_category(g_vec_7, s_vec_7, CAT7_SHORT)

    write_tests_tsv(args.output_tsv, args.group_name, g_total, s_total,
                    len(args.somatic), rows96, rows7)

    sig96_set = {r['category'] for r in rows96 if r['significant'] == 'yes'}
    sig7_set = {r['category'] for r in rows7 if r['significant'] == 'yes'}

    # Bootstrap at 96-category level
    g_draws = bootstrap_draws(g_vec, args.n_bootstrap, rng)  # (B, 96)
    s_draws = bootstrap_draws(s_vec, args.n_bootstrap, rng)

    # Collapse to 7 mutation types, splitting C>T into CpG vs non-CpG.
    # Category layout in the 96-vec: for mut in MUT_TYPES, for five in BASES, for three in BASES.
    # Reshape to (B, 6, 4, 4) with last two axes = (5' base, 3' base).
    # CpG C>T = C>T rows with 3' == G (three index 2).
    def collapse_7(draws):
        d = draws.reshape(-1, 6, 4, 4)
        ct_idx = MUT_TYPES.index('C>T')  # 2
        g_idx = BASES.index('G')         # 2
        ct_cpg = d[:, ct_idx, :, g_idx].sum(axis=-1)              # shape (B,)
        ct_all = d[:, ct_idx].sum(axis=(-2, -1))
        ct_non = ct_all - ct_cpg
        out = np.zeros((d.shape[0], 7))
        out[:, 0] = d[:, MUT_TYPES.index('C>A')].sum(axis=(-2, -1))
        out[:, 1] = d[:, MUT_TYPES.index('C>G')].sum(axis=(-2, -1))
        out[:, 2] = ct_non
        out[:, 3] = ct_cpg
        out[:, 4] = d[:, MUT_TYPES.index('T>A')].sum(axis=(-2, -1))
        out[:, 5] = d[:, MUT_TYPES.index('T>C')].sum(axis=(-2, -1))
        out[:, 6] = d[:, MUT_TYPES.index('T>G')].sum(axis=(-2, -1))
        return out

    g_draws_7 = collapse_7(g_draws)
    s_draws_7 = collapse_7(s_draws)

    COLLAPSED_LABELS = ['C>A', 'C>G', 'C>T\n(non-CpG)', 'C>T\n(CpG)', 'T>A', 'T>C', 'T>G']
    COLLAPSED_COLORS = [MUT_COLORS['C>A'], MUT_COLORS['C>G'],
                        MUT_COLORS['C>T'], MUT_COLORS['C>T'],
                        MUT_COLORS['T>A'], MUT_COLORS['T>C'], MUT_COLORS['T>G']]

    g_mean96, g_lo96, g_hi96 = summarize(g_draws)
    s_mean96, s_lo96, s_hi96 = summarize(s_draws)
    g_mean7, g_lo7, g_hi7 = summarize(g_draws_7)
    s_mean7, s_lo7, s_hi7 = summarize(s_draws_7)

    g_mean96 = g_mean96.reshape(6, 16)
    g_lo96 = g_lo96.reshape(6, 16)
    g_hi96 = g_hi96.reshape(6, 16)
    s_mean96 = s_mean96.reshape(6, 16)
    s_lo96 = s_lo96.reshape(6, 16)
    s_hi96 = s_hi96.reshape(6, 16)

    ymax96 = max(g_hi96.max(), s_hi96.max(), 0.01) * 1.1
    ymax7 = max(g_hi7.max(), s_hi7.max(), 0.01) * 1.1

    # Figure with 2 rows via gridspec: top row = 6 panels, bottom row = 1 panel
    fig = plt.figure(figsize=(20, 7))
    gs = gridspec.GridSpec(2, 6, height_ratios=[2.5, 1.2], hspace=0.5, wspace=0.15)

    top_axes = [fig.add_subplot(gs[0, i]) for i in range(6)]
    bottom_ax = fig.add_subplot(gs[1, :])

    width = 0.4
    x_pos = np.arange(16)

    # --- Top: 96 categories ---
    for i, mut in enumerate(MUT_TYPES):
        ax = top_axes[i]
        color = MUT_COLORS[mut]

        g_err = np.vstack([g_mean96[i] - g_lo96[i], g_hi96[i] - g_mean96[i]])
        ax.bar(x_pos - width / 2, g_mean96[i], width=width,
               color=color, edgecolor='black', linewidth=0.3,
               yerr=g_err, error_kw=dict(ecolor='black', lw=0.6, capsize=1.2))

        s_err = np.vstack([s_mean96[i] - s_lo96[i], s_hi96[i] - s_mean96[i]])
        ax.bar(x_pos + width / 2, s_mean96[i], width=width,
               color=color, edgecolor='black', linewidth=0.3, hatch='///',
               alpha=0.55,
               yerr=s_err, error_kw=dict(ecolor='black', lw=0.6, capsize=1.2))

        title_color = color if mut != 'C>G' else 'black'
        ax.set_title(mut, fontsize=12, color=title_color)
        ax.set_xticks(x_pos)
        trinucs = [f'{five}{mut[0]}{three}' for five in BASES for three in BASES]
        labels_with_stars = [
            f'{t}*' if f'{t[0]}[{mut}]{t[2]}' in sig96_set else t for t in trinucs
        ]
        tick_labels = ax.set_xticklabels(labels_with_stars, rotation=90, fontsize=6, family='monospace')
        for tl, trinuc in zip(tick_labels, trinucs):
            cat_label = f'{trinuc[0]}[{mut}]{trinuc[2]}'
            if cat_label in sig96_set:
                tl.set_fontweight('bold')
        ax.set_ylim(0, ymax96)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if i == 0:
            ax.set_ylabel('Proportion (per context)')
        else:
            ax.tick_params(labelleft=False)

    # --- Bottom: 7 collapsed mutation types (C>T split into CpG vs non-CpG) ---
    x7 = np.arange(7)
    w7 = 0.38
    g_err7 = np.vstack([g_mean7 - g_lo7, g_hi7 - g_mean7])
    bottom_ax.bar(x7 - w7 / 2, g_mean7, width=w7,
                  color=COLLAPSED_COLORS, edgecolor='black', linewidth=0.4,
                  yerr=g_err7, error_kw=dict(ecolor='black', lw=0.8, capsize=3))
    s_err7 = np.vstack([s_mean7 - s_lo7, s_hi7 - s_mean7])
    bottom_ax.bar(x7 + w7 / 2, s_mean7, width=w7,
                  color=COLLAPSED_COLORS, edgecolor='black', linewidth=0.4,
                  hatch='///', alpha=0.55,
                  yerr=s_err7, error_kw=dict(ecolor='black', lw=0.8, capsize=3))

    bottom_ax.set_xticks(x7)
    bottom_labels_with_stars = [
        f'{lbl}*' if short in sig7_set else lbl
        for lbl, short in zip(COLLAPSED_LABELS, CAT7_SHORT)
    ]
    bottom_tick_labels = bottom_ax.set_xticklabels(bottom_labels_with_stars, fontsize=10)
    for tl, short in zip(bottom_tick_labels, CAT7_SHORT):
        if short in sig7_set:
            tl.set_fontweight('bold')
    bottom_ax.set_ylim(0, ymax7)
    bottom_ax.set_ylabel('Proportion (mutation type)')
    bottom_ax.set_title('Collapsed mutation-type comparison (C>T split by CpG context)', fontsize=11)
    bottom_ax.spines['top'].set_visible(False)
    bottom_ax.spines['right'].set_visible(False)

    # Legend
    germline_handle = plt.Rectangle((0, 0), 1, 1, facecolor='grey',
                                     edgecolor='black', linewidth=0.3)
    somatic_handle = plt.Rectangle((0, 0), 1, 1, facecolor='grey',
                                    edgecolor='black', linewidth=0.3,
                                    hatch='///', alpha=0.55)
    fig.legend([germline_handle, somatic_handle],
               [f'Germline {args.group_name} (n={g_total})',
                f'Somatic {args.group_name} (n={s_total})'],
               loc='upper right', fontsize=9, frameon=False)

    fig.suptitle(f'Germline (autosome) vs Somatic (merged) — {args.group_name} — '
                 f'95% bootstrap CI ({args.n_bootstrap} reps)',
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 0.92, 0.95])

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    fig.savefig(args.output, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'Comparison plot -> {args.output}')
    print(f'Tests TSV       -> {args.output_tsv}')
    print(f'Significant categories (BH<0.05): 96-level={len(sig96_set)}, 7-level={len(sig7_set)}')


if __name__ == '__main__':
    main()
