#!/usr/bin/env python3
"""
summarize_bootstrap.py - Summarize bootstrap dN/dS results with confidence intervals

Reads the auto_all result (point estimate) and all bootstrap replicate results,
then computes 95% confidence intervals for dN, dS, and dN/dS.

Handles both branch model and pairwise model results, saving to separate subfolders.

Usage:
    python3 summarize_bootstrap.py \
        --bootstrap_dir steps/08_bootstrap/ \
        --n_bootstrap 500 \
        --output_dir steps/09_summary/
"""

import argparse
import os
import sys


def read_branch_tab(filepath):
    """Read a _branch.tab file and return dict: {branch: {dN, dS, dNdS}}."""
    data = {}
    with open(filepath) as f:
        header = next(f)
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue
            rep_id, branch, dN, dS, dNdS = fields
            data[branch] = {
                'dN': float(dN),
                'dS': float(dS),
                'dNdS': float(dNdS) if dNdS != 'NA' else None
            }
    return data


def read_pairwise_tab(filepath):
    """Read a _pairwise.tab file and return dict: {pair: {dN, dS, dNdS}}."""
    data = {}
    with open(filepath) as f:
        header = next(f)
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 6:
                continue
            rep_id, sp1, sp2, dN, dS, dNdS = fields
            pair = f"{sp1}_{sp2}"
            data[pair] = {
                'dN': float(dN),
                'dS': float(dS),
                'dNdS': float(dNdS) if dNdS != 'NA' else None
            }
    return data


def percentile(values, p):
    """Compute percentile using linear interpolation."""
    sorted_vals = sorted(values)
    n = len(sorted_vals)
    k = (n - 1) * p / 100.0
    f = int(k)
    c = f + 1
    if c >= n:
        return sorted_vals[-1]
    return sorted_vals[f] + (k - f) * (sorted_vals[c] - sorted_vals[f])


def summarize_mode(bootstrap_dir, n_bootstrap, output_dir, mode):
    """Summarize results for one mode (branch or pairwise).

    Args:
        bootstrap_dir: path to 08_bootstrap/
        n_bootstrap: number of bootstrap replicates
        output_dir: output subdirectory (e.g. 09_summary/branch/)
        mode: 'branch' or 'pairwise'
    """
    os.makedirs(output_dir, exist_ok=True)

    if mode == 'branch':
        read_func = read_branch_tab
        suffix = '_branch.tab'
        group_label = 'branch'
    else:
        read_func = read_pairwise_tab
        suffix = '_pairwise.tab'
        group_label = 'pair'

    # Read point estimate (auto_all)
    all_tab = os.path.join(bootstrap_dir, 'auto_all', mode, f'auto_all{suffix}')
    if not os.path.isfile(all_tab):
        print(f"WARNING: missing {all_tab}, skipping {mode} summary")
        return
    point_est = read_func(all_tab)
    groups = list(point_est.keys())
    print(f"  {mode} groups: {groups}")

    # Read bootstrap replicates
    bs_data = {g: {'dN': [], 'dS': [], 'dNdS': []} for g in groups}
    n_ok = 0
    n_fail = 0

    for i in range(1, n_bootstrap + 1):
        rep_id = f'auto_bs_{i}'
        tab_path = os.path.join(bootstrap_dir, rep_id, mode, f'{rep_id}{suffix}')
        if not os.path.isfile(tab_path):
            n_fail += 1
            continue
        try:
            rep_data = read_func(tab_path)
            for g in groups:
                if g in rep_data:
                    bs_data[g]['dN'].append(rep_data[g]['dN'])
                    bs_data[g]['dS'].append(rep_data[g]['dS'])
                    if rep_data[g]['dNdS'] is not None:
                        bs_data[g]['dNdS'].append(rep_data[g]['dNdS'])
            n_ok += 1
        except Exception as e:
            print(f"  WARNING: failed to read {tab_path}: {e}")
            n_fail += 1

    print(f"  Bootstrap replicates: {n_ok} ok, {n_fail} failed")

    # Write summary table
    summary_path = os.path.join(output_dir, f'auto_{mode}_summary.tsv')
    with open(summary_path, 'w') as out:
        out.write(f"{group_label}\tmetric\testimate\tbs_mean\tbs_median\tbs_ci_low\tbs_ci_high\tn_bs\n")
        for g in groups:
            for metric in ['dN', 'dS', 'dNdS']:
                est = point_est[g].get(metric)
                vals = bs_data[g][metric]
                if not vals:
                    out.write(f"{g}\t{metric}\t{est}\tNA\tNA\tNA\tNA\t0\n")
                    continue
                bs_mean = sum(vals) / len(vals)
                bs_median = percentile(vals, 50)
                ci_low = percentile(vals, 2.5)
                ci_high = percentile(vals, 97.5)
                est_str = f"{est:.6f}" if est is not None else "NA"
                out.write(f"{g}\t{metric}\t{est_str}\t{bs_mean:.6f}\t{bs_median:.6f}\t"
                          f"{ci_low:.6f}\t{ci_high:.6f}\t{len(vals)}\n")

    print(f"  Summary: {summary_path}")

    # Write all bootstrap values (long format)
    long_path = os.path.join(output_dir, f'auto_{mode}_bootstrap_values.tsv')
    with open(long_path, 'w') as out:
        out.write(f"replicate\t{group_label}\tdN\tdS\tdNdS\n")
        for i in range(1, n_bootstrap + 1):
            rep_id = f'auto_bs_{i}'
            tab_path = os.path.join(bootstrap_dir, rep_id, mode, f'{rep_id}{suffix}')
            if not os.path.isfile(tab_path):
                continue
            try:
                rep_data = read_func(tab_path)
                for g in groups:
                    if g in rep_data:
                        d = rep_data[g]
                        dNdS_str = f"{d['dNdS']:.6f}" if d['dNdS'] is not None else "NA"
                        out.write(f"{rep_id}\t{g}\t{d['dN']:.6f}\t{d['dS']:.6f}\t{dNdS_str}\n")
            except Exception:
                pass

    print(f"  Bootstrap values: {long_path}")


def main():
    parser = argparse.ArgumentParser(description='Summarize bootstrap dN/dS results')
    parser.add_argument('--bootstrap_dir', required=True,
                        help='Directory containing auto_all/ and auto_bs_*/ subdirs')
    parser.add_argument('--n_bootstrap', type=int, required=True,
                        help='Number of bootstrap replicates')
    parser.add_argument('--output_dir', required=True,
                        help='Output directory for summary files')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Summarize branch model results
    print("Summarizing branch model results...")
    summarize_mode(
        args.bootstrap_dir, args.n_bootstrap,
        os.path.join(args.output_dir, 'branch'), 'branch'
    )

    # Summarize pairwise model results
    print("Summarizing pairwise model results...")
    summarize_mode(
        args.bootstrap_dir, args.n_bootstrap,
        os.path.join(args.output_dir, 'pairwise'), 'pairwise'
    )

    print("Done.")


if __name__ == '__main__':
    main()
