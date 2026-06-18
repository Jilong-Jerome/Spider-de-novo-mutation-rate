#!/usr/bin/env python3
"""
Merge per-individual BAM coverage TSVs and per-chromosome VCF DP TSVs into a single
tidy per-species report: one row per (individual x chromosome).

BAM stats are joined to VCF stats on (sample, chrom). The join is an inner join on the
VCF sample set (families 1-4), so BAM-only individuals (e.g. family5) are dropped, as
agreed. Output is sorted by sample then chrom.
"""
import argparse
import glob
import os
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam-dir", required=True, help="dir with *.bam_cov.tsv")
    ap.add_argument("--vcf-dir", required=True, help="dir with *.vcf_dp.tsv")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    bam_files = sorted(glob.glob(os.path.join(args.bam_dir, "*.bam_cov.tsv")))
    vcf_files = sorted(glob.glob(os.path.join(args.vcf_dir, "*.vcf_dp.tsv")))
    if not bam_files:
        raise SystemExit(f"No BAM cov files in {args.bam_dir}")
    if not vcf_files:
        raise SystemExit(f"No VCF dp files in {args.vcf_dir}")

    bam = pd.concat((pd.read_csv(f, sep="\t") for f in bam_files), ignore_index=True)
    vcf = pd.concat((pd.read_csv(f, sep="\t") for f in vcf_files), ignore_index=True)

    merged = pd.merge(bam, vcf, on=["sample", "chrom"], how="inner")
    merged = merged.sort_values(["sample", "chrom"]).reset_index(drop=True)

    cols = ["sample", "chrom", "length",
            "covered_bases", "sum_depth", "mean_all_depth",
            "mean_covered", "median_covered",
            "n_sites_all", "mean_DP_all", "median_DP_all",
            "n_sites_DPge1", "mean_DP_ge1", "median_DP_ge1"]
    merged = merged[cols]
    merged.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {len(merged)} rows ({merged['sample'].nunique()} individuals x "
          f"{merged['chrom'].nunique()} chromosomes) to {args.out}")


if __name__ == "__main__":
    main()
