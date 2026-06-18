#!/usr/bin/env python3
"""
Per-sample called-DP statistics for one chromosome from a `bcftools query` stream.

Reads lines on stdin produced by:
    bcftools query -r <chrom> -f '[%DP\\t]\\n' <vcf>
i.e. one line per site, with one tab-separated DP column per sample in VCF sample
order (a trailing tab is fine). Missing values ('.') are treated as 0 (the VCF is
all-sites, so every position is a record).

For each sample, two depth histograms are accumulated and summarised:
  (a) all sites, including DP=0        -> mean_DP_all, median_DP_all   (requirement 3a)
  (b) only sites with DP >= 1          -> mean_DP_ge1, median_DP_ge1   (requirement 3b)
"""
import argparse
import sys
from collections import defaultdict


def median_from_hist(hist, n):
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
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--samples", required=True, help="comma-separated sample names in VCF order")
    ap.add_argument("--out", required=True)
    ap.add_argument("--hist-out", required=False,
                    help="optional TSV (sample<TAB>dp<TAB>count) of the all-sites DP "
                         "histogram; consumed by autosome_dp_report.py to pool autosomes")
    args = ap.parse_args()

    samples = args.samples.split(",")
    nsamp = len(samples)

    hist = [defaultdict(int) for _ in range(nsamp)]   # all sites (incl DP=0)
    sum_all = [0] * nsamp
    n_all = [0] * nsamp

    for line in sys.stdin:
        line = line.rstrip("\n")
        if not line:
            continue
        cols = line.split("\t")
        # tolerate a trailing empty field from the trailing tab in the format string
        if len(cols) > nsamp and cols[-1] == "":
            cols = cols[:nsamp]
        if len(cols) < nsamp:
            sys.exit(f"ERROR: site has {len(cols)} DP columns, expected {nsamp}: {line[:80]!r}")
        for i in range(nsamp):
            tok = cols[i]
            d = 0 if (tok == "." or tok == "") else int(tok)
            hist[i][d] += 1
            sum_all[i] += d
            n_all[i] += 1

    with open(args.out, "w") as out:
        out.write("sample\tchrom\tn_sites_all\tmean_DP_all\tmedian_DP_all\t"
                  "n_sites_DPge1\tmean_DP_ge1\tmedian_DP_ge1\n")
        for i, sample in enumerate(samples):
            na = n_all[i]
            # all-sites stats
            mean_all = (sum_all[i] / na) if na > 0 else 0.0
            med_all = median_from_hist(hist[i], na) if na > 0 else None
            # DP>=1 stats derived from the same histogram
            n_ge1 = na - hist[i].get(0, 0)
            sum_ge1 = sum_all[i]  # DP=0 contributes nothing to the sum
            if n_ge1 > 0:
                mean_ge1 = sum_ge1 / n_ge1
                ge1_hist = {k: v for k, v in hist[i].items() if k >= 1}
                med_ge1 = median_from_hist(ge1_hist, n_ge1)
                mean_ge1_s = f"{mean_ge1:.6f}"
                med_ge1_s = f"{med_ge1:.1f}"
            else:
                mean_ge1_s = "NA"
                med_ge1_s = "NA"
            med_all_s = f"{med_all:.1f}" if med_all is not None else "NA"
            out.write(f"{sample}\t{args.chrom}\t{na}\t{mean_all:.6f}\t{med_all_s}\t"
                      f"{n_ge1}\t{mean_ge1_s}\t{med_ge1_s}\n")

    # Side-output: the raw all-sites DP histogram per sample, so the autosome
    # report can pool counts across chromosomes and compute an exact pooled median.
    if args.hist_out:
        with open(args.hist_out, "w") as hf:
            hf.write("sample\tdp\tcount\n")
            for i, sample in enumerate(samples):
                for dp in sorted(hist[i]):
                    hf.write(f"{sample}\t{dp}\t{hist[i][dp]}\n")


if __name__ == "__main__":
    main()
