#!/usr/bin/env python3
"""
Per-individual autosomal DP report for one species, from VCF DP over ALL sites
(DP=0 included).

Reads the per-chromosome all-sites DP histograms written by vcf_dp_stats.py
(`{chrom}.vcf_dp_hist.tsv`, columns sample/dp/count) for the AUTOSOMAL
chromosomes only, pools the counts per (sample, dp) across those chromosomes, and
writes one row per individual:

    sample   n_sites   mean_DP_all   median_DP_all

mean and median are computed over all autosomal sites including DP=0. The pooled
histogram lets the median be exact (it cannot be reconstructed from per-chromosome
medians).
"""
import argparse
import os
from collections import defaultdict


def median_from_hist(hist, n):
    """Median of a {value: count} histogram covering n observations (n > 0)."""
    if n <= 0:
        return None
    lo_idx = (n - 1) // 2
    hi_idx = n // 2
    lo_val = hi_val = None
    cum = 0
    for v in sorted(hist):
        cum += hist[v]
        if lo_val is None and cum > lo_idx:
            lo_val = v
        if cum > hi_idx:
            hi_val = v
            break
    return (lo_val + hi_val) / 2.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf-dir", required=True, help="dir with {chrom}.vcf_dp_hist.tsv")
    ap.add_argument("--autosomes", required=True,
                    help="comma-separated autosomal chromosome names")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    autosomes = [c for c in args.autosomes.split(",") if c]

    # pooled per-sample histogram over all autosomal sites
    hist = defaultdict(lambda: defaultdict(int))   # sample -> {dp: count}
    samples_order = []
    seen = set()

    for chrom in autosomes:
        path = os.path.join(args.vcf_dir, f"{chrom}.vcf_dp_hist.tsv")
        if not os.path.exists(path):
            raise SystemExit(f"ERROR: missing histogram file {path}")
        with open(path) as fh:
            header = fh.readline()  # sample\tdp\tcount
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                sample, dp, count = line.split("\t")
                if sample not in seen:
                    seen.add(sample)
                    samples_order.append(sample)
                hist[sample][int(dp)] += int(count)

    with open(args.out, "w") as out:
        out.write("sample\tn_sites\tmean_DP_all\tmedian_DP_all\n")
        for sample in samples_order:
            h = hist[sample]
            n = sum(h.values())
            sum_dp = sum(dp * c for dp, c in h.items())
            mean_all = (sum_dp / n) if n > 0 else 0.0
            med_all = median_from_hist(h, n)
            med_s = f"{med_all:.1f}" if med_all is not None else "NA"
            out.write(f"{sample}\t{n}\t{mean_all:.6f}\t{med_s}\n")

    print(f"Wrote {len(samples_order)} individuals to {args.out} "
          f"(pooled over {len(autosomes)} autosomes)")


if __name__ == "__main__":
    main()
