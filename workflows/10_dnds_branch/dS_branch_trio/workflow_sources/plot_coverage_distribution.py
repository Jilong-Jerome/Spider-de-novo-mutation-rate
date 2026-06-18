#!/usr/bin/env python3
"""
plot_coverage_distribution.py - Per-trio coverage-distribution diagnostic.

Consumes the per-species covered-site depth histograms and coverage-stats TSVs
produced by callable_depth_regions.py (Step 05a) and produces, per trio:

  1. coverage_distribution.{png,pdf} : one panel per species, histogram of
     covered-site depth with reference lines for mean, median, and the current
     callable_min/callable_max bounds.
  2. coverage_distribution_summary.tsv : mean vs median (and percentiles),
     callable bounds, and the fraction of covered bases inside the callable
     range, one row per species.

This is a reference-only diagnostic; it does not feed downstream steps.

Usage:
    python3 plot_coverage_distribution.py \
        --species-distributions MIM=MIM_depth_distribution.tsv,AFR=...,BIC=... \
        --species-stats MIM=MIM_coverage_stats.tsv,AFR=...,BIC=... \
        --output-dir steps/MIM_AFR_BIC/05_callable_depth/coverage_distribution \
        --trio-name MIM_AFR_BIC
"""

import argparse
import os

import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_species_map(value):
    """Parse 'sp1=path1,sp2=path2' into an ordered list of (species, path)."""
    pairs = []
    for item in value.split(","):
        item = item.strip()
        if not item:
            continue
        if "=" not in item:
            raise ValueError(f"expected 'species=path', got '{item}'")
        species, path = item.split("=", 1)
        pairs.append((species.strip(), path.strip()))
    return pairs


def read_stats(path):
    """Read a metric/value coverage-stats TSV into a dict (values kept as str)."""
    stats = {}
    with open(path) as handle:
        next(handle, None)  # header
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            metric, value = line.split("\t")[:2]
            stats[metric] = value
    return stats


def main():
    parser = argparse.ArgumentParser(
        description="Per-trio coverage-distribution diagnostic plots")
    parser.add_argument("--species-distributions", required=True,
                        help="Comma list of species=depth_distribution.tsv")
    parser.add_argument("--species-stats", required=True,
                        help="Comma list of species=coverage_stats.tsv")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--trio-name", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    dist_map = dict(parse_species_map(args.species_distributions))
    stats_map = dict(parse_species_map(args.species_stats))
    species_list = [sp for sp, _ in parse_species_map(args.species_distributions)]

    summary_rows = []
    n = len(species_list)
    fig, axes = plt.subplots(n, 1, figsize=(9, 3.2 * n), squeeze=False)

    for idx, species in enumerate(species_list):
        dist = pd.read_csv(dist_map[species], sep="\t")
        stats = read_stats(stats_map[species])

        mean_depth = float(stats["covered_mean_depth"])
        median_depth = float(stats["covered_median_depth"])
        callable_min = int(stats["callable_min_depth"])
        callable_max = int(stats["callable_max_depth"])
        p05 = float(stats.get("covered_depth_p05", "nan"))
        p25 = float(stats.get("covered_depth_p25", "nan"))
        p75 = float(stats.get("covered_depth_p75", "nan"))
        p95 = float(stats.get("covered_depth_p95", "nan"))

        in_range = dist.loc[
            (dist["depth"] >= callable_min) & (dist["depth"] <= callable_max),
            "frac_of_covered",
        ].sum()

        summary_rows.append({
            "trio": args.trio_name,
            "species": species,
            "covered_mean_depth": round(mean_depth, 4),
            "covered_median_depth": median_depth,
            "p05": p05, "p25": p25, "p75": p75, "p95": p95,
            "callable_min_depth": callable_min,
            "callable_max_depth": callable_max,
            "frac_covered_in_callable_range": round(float(in_range), 6),
        })

        # x-axis cap for readability: max of p95 and callable_max, padded.
        x_cap = max(p95, callable_max) * 1.1 if p95 == p95 else callable_max * 1.5
        ax = axes[idx][0]
        ax.bar(dist["depth"], dist["n_bases"], width=1.0,
               color="#4878a8", linewidth=0)
        ax.set_xlim(0, x_cap)
        ax.set_title(f"{species} - covered-site depth distribution")
        ax.set_xlabel("Depth")
        ax.set_ylabel("Bases")

        ax.axvline(median_depth, color="black", linestyle="-", linewidth=1.5,
                   label=f"median = {median_depth:g}")
        ax.axvline(mean_depth, color="crimson", linestyle="--", linewidth=1.5,
                   label=f"mean = {mean_depth:.1f}")
        ax.axvline(callable_min, color="green", linestyle=":", linewidth=1.5,
                   label=f"callable_min = {callable_min}")
        ax.axvline(callable_max, color="green", linestyle=":", linewidth=1.5,
                   label=f"callable_max = {callable_max}")
        ax.legend(fontsize=8)

    fig.suptitle(f"{args.trio_name} coverage distribution (covered sites)",
                 y=1.0)
    fig.tight_layout()
    png_path = os.path.join(args.output_dir, "coverage_distribution.png")
    pdf_path = os.path.join(args.output_dir, "coverage_distribution.pdf")
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    summary_path = os.path.join(args.output_dir,
                                "coverage_distribution_summary.tsv")
    pd.DataFrame(summary_rows).to_csv(summary_path, sep="\t", index=False)

    print(f"Written: {png_path}")
    print(f"Written: {pdf_path}")
    print(f"Written: {summary_path}")


if __name__ == "__main__":
    main()
