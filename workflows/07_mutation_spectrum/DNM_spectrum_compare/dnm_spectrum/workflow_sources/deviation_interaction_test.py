#!/usr/bin/env python3
"""
deviation_interaction_test.py
Test whether the germline-vs-somatic deviation per mutation category differs
between the social and subsocial species groups.

Per category (7 collapsed categories):
    Breslow-Day test for homogeneity of the germline-vs-somatic odds ratio
    across the two groups (social vs subsocial). BH-adjusted across 7.

Global:
    3-way log-linear interaction test on the 2x2x7 table
    (group x tissue x category), via Poisson GLM with and without the
    group:tissue:category term; LR test on deviance difference.
"""
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

try:
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    import pandas as pd
    HAVE_SM = True
except Exception:
    HAVE_SM = False

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
    """Collapse a {(mut, trinuc)->count} dict to a length-7 numpy array."""
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


def breslow_day(tables):
    """Breslow-Day test for homogeneity of odds ratios across K 2x2 tables.
    Each table is [[a, b], [c, d]]. Uses Tarone correction.
    Returns (chi2, dof, p).
    """
    tables = [np.asarray(t, dtype=float) for t in tables]
    K = len(tables)
    # Mantel-Haenszel common odds ratio
    R = sum(t[0, 0] * t[1, 1] / t.sum() for t in tables)
    S = sum(t[0, 1] * t[1, 0] / t.sum() for t in tables)
    if S == 0 or R == 0:
        return (np.nan, K - 1, np.nan)
    or_mh = R / S

    chi2_sum = 0.0
    var_sum = 0.0
    for t in tables:
        a, b = t[0]
        c, d = t[1]
        n1 = a + b
        n2 = c + d
        m1 = a + c
        m2 = b + d
        N = n1 + n2
        # Solve expected a_hat: or_mh = a_hat*(n2-m1+a_hat) / ((n1-a_hat)*(m1-a_hat))
        # Quadratic: (or_mh-1)*a^2 + (n1+m1 - (or_mh-1)*(n1+m1) ... ) — use standard form
        A = or_mh - 1.0
        B = -(or_mh * (n1 + m1) + (n2 - m1))
        C = or_mh * n1 * m1
        if abs(A) < 1e-12:
            a_hat = n1 * m1 / N
        else:
            disc = B * B - 4 * A * C
            a_hat = (-B - np.sqrt(max(disc, 0))) / (2 * A)
            if a_hat < max(0, m1 - n2) or a_hat > min(n1, m1):
                a_hat = (-B + np.sqrt(max(disc, 0))) / (2 * A)
        # variance of a under common OR
        var = 1.0 / (1.0 / a_hat + 1.0 / (n1 - a_hat) +
                     1.0 / (m1 - a_hat) + 1.0 / (n2 - m1 + a_hat))
        chi2_sum += (a - a_hat) ** 2 / var
        var_sum += var
    # Tarone correction
    diff_sum = sum(t[0, 0] - _a_hat_single(t, or_mh) for t in tables)
    tarone = chi2_sum - diff_sum * diff_sum / var_sum if var_sum > 0 else chi2_sum
    dof = K - 1
    p = 1 - stats.chi2.cdf(tarone, dof)
    return (tarone, dof, p)


def _a_hat_single(t, or_mh):
    a, b = t[0]
    c, d = t[1]
    n1 = a + b
    n2 = c + d
    m1 = a + c
    N = n1 + n2
    A = or_mh - 1.0
    B = -(or_mh * (n1 + m1) + (n2 - m1))
    C = or_mh * n1 * m1
    if abs(A) < 1e-12:
        return n1 * m1 / N
    disc = B * B - 4 * A * C
    root = (-B - np.sqrt(max(disc, 0))) / (2 * A)
    if root < max(0, m1 - n2) or root > min(n1, m1):
        root = (-B + np.sqrt(max(disc, 0))) / (2 * A)
    return root


def odds_ratio(a, b, c, d, eps=0.5):
    """Haldane-Anscombe corrected OR for safety in bootstrap."""
    a2, b2, c2, d2 = a + eps, b + eps, c + eps, d + eps
    return (a2 * d2) / (b2 * c2)


def global_3way_test(g_soc, s_soc, g_sub, s_sub):
    """Poisson log-linear LR test for the 3-way interaction."""
    if not HAVE_SM:
        return _global_3way_ipf(g_soc, s_soc, g_sub, s_sub)
    rows = []
    for i, cat in enumerate(CAT7_SHORT):
        rows.append({'group': 'social', 'tissue': 'germline',
                     'category': cat, 'count': int(g_soc[i])})
        rows.append({'group': 'social', 'tissue': 'somatic',
                     'category': cat, 'count': int(s_soc[i])})
        rows.append({'group': 'subsocial', 'tissue': 'germline',
                     'category': cat, 'count': int(g_sub[i])})
        rows.append({'group': 'subsocial', 'tissue': 'somatic',
                     'category': cat, 'count': int(s_sub[i])})
    df = pd.DataFrame(rows)
    reduced = smf.glm(
        'count ~ group*tissue + group*category + tissue*category',
        data=df, family=sm.families.Poisson()).fit()
    full = smf.glm(
        'count ~ group*tissue*category',
        data=df, family=sm.families.Poisson()).fit()
    lr = 2 * (full.llf - reduced.llf)
    dof = int(reduced.df_resid - full.df_resid)
    p = 1 - stats.chi2.cdf(lr, dof)
    return (lr, dof, p)


def _global_3way_ipf(g_soc, s_soc, g_sub, s_sub, max_iter=500, tol=1e-8):
    """IPF fallback: fit the no-3-way model, return G^2."""
    obs = np.zeros((2, 2, 7))
    obs[0, 0] = g_soc; obs[0, 1] = s_soc
    obs[1, 0] = g_sub; obs[1, 1] = s_sub
    m = np.ones_like(obs, dtype=float)
    GT = obs.sum(axis=2)
    GC = obs.sum(axis=1)
    TC = obs.sum(axis=0)
    for _ in range(max_iter):
        prev = m.copy()
        m *= (GT / m.sum(axis=2))[:, :, None]
        m *= (GC / m.sum(axis=1))[:, None, :]
        m *= (TC / m.sum(axis=0))[None, :, :]
        if np.max(np.abs(m - prev)) < tol:
            break
    with np.errstate(divide='ignore', invalid='ignore'):
        terms = np.where(obs > 0, obs * np.log(obs / m), 0.0)
    g2 = 2 * terms.sum()
    dof = 6
    p = 1 - stats.chi2.cdf(g2, dof)
    return (g2, dof, p)


def multinomial_bootstrap(counts, n_boot, rng):
    total = counts.sum()
    if total == 0:
        return np.zeros((n_boot, len(counts)))
    p = counts / total
    return rng.multinomial(int(total), p, size=n_boot).astype(int)


def compute_per_category(g_soc, s_soc, g_sub, s_sub):
    ng_soc, ns_soc = g_soc.sum(), s_soc.sum()
    ng_sub, ns_sub = g_sub.sum(), s_sub.sum()
    rows, pvals = [], []
    for i in range(7):
        a, b = int(g_soc[i]), int(ns_soc) - int(s_soc[i])  # placeholder
        # Two 2x2 tables stratified by group:
        #   social:    [[g_soc_i, ng_soc - g_soc_i], [s_soc_i, ns_soc - s_soc_i]]
        #   subsocial: [[g_sub_i, ng_sub - g_sub_i], [s_sub_i, ns_sub - s_sub_i]]
        t_soc = [[int(g_soc[i]), int(ng_soc) - int(g_soc[i])],
                 [int(s_soc[i]), int(ns_soc) - int(s_soc[i])]]
        t_sub = [[int(g_sub[i]), int(ng_sub) - int(g_sub[i])],
                 [int(s_sub[i]), int(ns_sub) - int(s_sub[i])]]
        or_soc = odds_ratio(*t_soc[0], *t_soc[1])
        or_sub = odds_ratio(*t_sub[0], *t_sub[1])
        chi2, dof, p = breslow_day([t_soc, t_sub])
        rows.append({
            'category': CAT7_SHORT[i],
            'g_social': int(g_soc[i]), 's_social': int(s_soc[i]),
            'g_subsocial': int(g_sub[i]), 's_subsocial': int(s_sub[i]),
            'OR_social': or_soc,
            'OR_subsocial': or_sub,
            'ratio_of_ORs': or_soc / or_sub,
            'bd_chi2': chi2,
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


def write_tsv(path, global_res, per_cat_rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write('# Global 3-way interaction test (group x tissue x category)\n')
        f.write('test\tstat\tdof\tp_value\n')
        f.write(f'loglinear_LR\t{global_res[0]:.4f}\t{global_res[1]}\t{global_res[2]:.4g}\n')
        f.write('\n# Per-category Breslow-Day test for homogeneity of germline/somatic OR\n')
        f.write('category\tg_social\ts_social\tg_subsocial\ts_subsocial\t'
                'OR_social\tOR_subsocial\tratio_of_ORs\tbd_chi2\tp_value\tp_adjusted\tsignificant\n')
        for r in per_cat_rows:
            f.write(
                f"{r['category']}\t{r['g_social']}\t{r['s_social']}\t"
                f"{r['g_subsocial']}\t{r['s_subsocial']}\t"
                f"{r['OR_social']:.4g}\t{r['OR_subsocial']:.4g}\t{r['ratio_of_ORs']:.4g}\t"
                f"{r['bd_chi2']:.4f}\t{r['p_value']:.4g}\t{r['p_adjusted']:.4g}\t{r['significant']}\n")


def bootstrap_log_ors(g, s, n_boot, rng):
    """Return (B, 7) array of log2 OR(germline/somatic) per category."""
    g_draws = multinomial_bootstrap(g, n_boot, rng)
    s_draws = multinomial_bootstrap(s, n_boot, rng)
    ng = g_draws.sum(axis=1, keepdims=True)
    ns = s_draws.sum(axis=1, keepdims=True)
    eps = 0.5
    num = (g_draws + eps) * (ns - s_draws + eps)
    den = (s_draws + eps) * (ng - g_draws + eps)
    return np.log2(num / den)


def plot_panels(g_soc, s_soc, g_sub, s_sub, per_cat_rows, global_res,
                out_pdf, seed=42, n_boot=1000):
    rng = np.random.default_rng(seed)
    lor_soc = bootstrap_log_ors(g_soc, s_soc, n_boot, rng)
    lor_sub = bootstrap_log_ors(g_sub, s_sub, n_boot, rng)

    soc_mean = lor_soc.mean(axis=0)
    soc_lo = np.quantile(lor_soc, 0.025, axis=0)
    soc_hi = np.quantile(lor_soc, 0.975, axis=0)
    sub_mean = lor_sub.mean(axis=0)
    sub_lo = np.quantile(lor_sub, 0.025, axis=0)
    sub_hi = np.quantile(lor_sub, 0.975, axis=0)

    ratio = lor_soc - lor_sub  # log2 of ratio of ORs
    r_mean = ratio.mean(axis=0)
    r_lo = np.quantile(ratio, 0.025, axis=0)
    r_hi = np.quantile(ratio, 0.975, axis=0)

    fig, axes = plt.subplots(2, 1, figsize=(11, 8))
    x = np.arange(7)
    w = 0.38

    # Top: log2 OR(germline/somatic) for each group
    ax = axes[0]
    soc_err = np.vstack([soc_mean - soc_lo, soc_hi - soc_mean])
    ax.bar(x - w / 2, soc_mean, width=w,
           color=CAT7_COLORS, edgecolor='black', linewidth=0.5,
           yerr=soc_err, error_kw=dict(ecolor='black', lw=0.8, capsize=3),
           label='social')
    sub_err = np.vstack([sub_mean - sub_lo, sub_hi - sub_mean])
    ax.bar(x + w / 2, sub_mean, width=w,
           color=CAT7_COLORS, edgecolor='black', linewidth=0.5,
           hatch='///', alpha=0.55,
           yerr=sub_err, error_kw=dict(ecolor='black', lw=0.8, capsize=3),
           label='subsocial')
    ax.axhline(0, color='black', lw=0.5)

    ymax = max(soc_hi.max(), sub_hi.max())
    ymin = min(soc_lo.min(), sub_lo.min())
    span = max(abs(ymax), abs(ymin))
    span = span * 1.35 if span > 0 else 1
    for i, r in enumerate(per_cat_rows):
        stars = sig_stars(r['p_adjusted'])
        top = max(soc_hi[i], sub_hi[i])
        ax.text(i, top + span * 0.05, stars, ha='center', va='bottom',
                fontsize=10, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(CAT7_LABELS, fontsize=10)
    ax.set_ylabel('log2 OR(germline / somatic)')
    ax.set_ylim(-span, span)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    gchi2, gdof, gp = global_res
    ax.set_title(f'Per-group germline-vs-somatic deviation  —  '
                 f'global 3-way interaction LR \u03c7\u00b2={gchi2:.2f}, '
                 f'df={gdof}, p={gp:.3g}', fontsize=11)

    # Bottom: ratio of ORs (log2)
    ax = axes[1]
    r_err = np.vstack([r_mean - r_lo, r_hi - r_mean])
    ax.bar(x, r_mean, width=0.6, color=CAT7_COLORS,
           edgecolor='black', linewidth=0.5,
           yerr=r_err, error_kw=dict(ecolor='black', lw=0.8, capsize=3))
    ax.axhline(0, color='black', lw=0.5)
    rspan = max(abs(r_lo.min()), abs(r_hi.max())) * 1.3
    rspan = rspan if rspan > 0 else 1
    ax.set_xticks(x)
    ax.set_xticklabels(CAT7_LABELS, fontsize=10)
    ax.set_ylabel('log2 (OR_social / OR_subsocial)')
    ax.set_ylim(-rspan, rspan)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title('Ratio of germline/somatic ORs: social vs subsocial',
                 fontsize=11)

    soc_handle = plt.Rectangle((0, 0), 1, 1, facecolor='grey',
                                edgecolor='black', linewidth=0.5)
    sub_handle = plt.Rectangle((0, 0), 1, 1, facecolor='grey',
                                edgecolor='black', linewidth=0.5,
                                hatch='///', alpha=0.55)
    fig.legend([soc_handle, sub_handle], ['Social', 'Subsocial'],
               loc='upper right', fontsize=10, frameon=False)
    fig.suptitle(
        'Germline-vs-somatic deviation: social vs subsocial\n'
        '(Breslow-Day per category, BH-adjusted: *p<0.05, **p<0.01, ***p<0.001; '
        '95% bootstrap CI)',
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 0.94, 0.93])
    os.makedirs(os.path.dirname(out_pdf), exist_ok=True)
    fig.savefig(out_pdf, dpi=200, bbox_inches='tight')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description='Social vs subsocial germline-vs-somatic deviation interaction test')
    parser.add_argument('--germline_social', required=True)
    parser.add_argument('--germline_subsocial', required=True)
    parser.add_argument('--somatic_social', nargs='+', required=True)
    parser.add_argument('--somatic_subsocial', nargs='+', required=True)
    parser.add_argument('--output_tsv', required=True)
    parser.add_argument('--output_pdf', required=True)
    args = parser.parse_args()

    g_soc = counts_to_7category(read_germline(args.germline_social))
    g_sub = counts_to_7category(read_germline(args.germline_subsocial))
    s_soc = counts_to_7category(merge_counts([read_somatic(p) for p in args.somatic_social]))
    s_sub = counts_to_7category(merge_counts([read_somatic(p) for p in args.somatic_subsocial]))

    per_cat = compute_per_category(g_soc, s_soc, g_sub, s_sub)
    global_res = global_3way_test(g_soc, s_soc, g_sub, s_sub)

    write_tsv(args.output_tsv, global_res, per_cat)
    plot_panels(g_soc, s_soc, g_sub, s_sub, per_cat, global_res, args.output_pdf)

    print(f'TSV -> {args.output_tsv}')
    print(f'PDF -> {args.output_pdf}')
    print(f'Global 3-way interaction: LR={global_res[0]:.3f}, '
          f'df={global_res[1]}, p={global_res[2]:.4g}')


if __name__ == '__main__':
    main()
