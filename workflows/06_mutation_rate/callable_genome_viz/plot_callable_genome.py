#!/usr/bin/env python3
"""
plot_callable_genome.py - Callable (autosomal) genome size per individual, by species.

Reads the per-individual, per-chromosome callable base counts in
``data/<SP>_callable_minDP26.tsv`` and the reference chromosome lengths in
``data/genome/<SP>_ncbi_chromosome.fa.fai`` and produces, for the 7 spider
species, a single faceted figure (one panel per species) showing each
individual's callable autosomal genome size.

Autosome-only: the callable tsv files contain no sex-chromosome rows, so the
reference total used as the denominator/reference line is autosomal as well
(sex chromosomes ``*_X*`` in the fai are excluded).

Outputs (written next to this script):
  - callable_genome_size.png / .pdf : faceted bar chart, one panel per species.
    Left y-axis = callable size (Mb); right y-axis = % of autosomal genome;
    dashed line = total autosomal genome length.
  - callable_genome_summary.tsv : one row per individual with the numbers
    behind the figure.

Usage:
    python3 plot_callable_genome.py
"""

import os
import re

import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")
GENOME_DIR = os.path.join(DATA_DIR, "genome")

# Species code (callable-file prefix) -> fai filename. The fai codes differ from
# the callable-file codes for BIC/SAR/TEN, so the mapping is explicit.
FAI_FILE = {
    "AFR": "AFR_ncbi_chromosome.fa.fai",
    "BIC": "BI_ncbi_chromosome.fa.fai",
    "DUM": "DUM_ncbi_chromosome.fa.fai",
    "LIN": "LIN_ncbi_chromosome.fa.fai",
    "MIM": "MIM_ncbi_chromosome.fa.fai",
    "SAR": "SARA_ncbi_chromosome.fa.fai",
    "TEN": "TENT_ncbi_chromosome.fa.fai",
}
SPECIES = ["AFR", "BIC", "DUM", "LIN", "MIM", "SAR", "TEN"]

BAR_COLOR = "#4C72B0"
REF_COLOR = "#C44E52"


def family_sample_key(individual):
    """Sort key: (family number, sample number) parsed from the individual name.

    Names look like ``AFR_family1_S3_offspring``; unparsed parts fall back to 0
    so ordering stays stable rather than raising.
    """
    fam = re.search(r"family(\d+)", individual)
    sam = re.search(r"_S(\d+)_", individual)
    return (int(fam.group(1)) if fam else 0, int(sam.group(1)) if sam else 0)


def load_species(sp):
    """Return (per-individual DataFrame, total autosomal length in Mb) for one species."""
    callable_path = os.path.join(DATA_DIR, f"{sp}_callable_minDP26.tsv")
    fai_path = os.path.join(GENOME_DIR, FAI_FILE[sp])

    df = pd.read_csv(
        callable_path,
        sep="\t",
        header=None,
        names=["chrom", "callable_bp", "ind_chrom"],
    )
    # individual = ind_chrom with the trailing "_<chrom>" stripped.
    df["individual"] = df.apply(
        lambda r: r["ind_chrom"][: -(len(r["chrom"]) + 1)]
        if r["ind_chrom"].endswith("_" + r["chrom"])
        else r["ind_chrom"],
        axis=1,
    )

    # Reference total: only the chromosomes that appear in the callable file
    # (autosomes), so sex chromosomes in the fai are dropped automatically.
    fai = pd.read_csv(
        fai_path, sep="\t", header=None, usecols=[0, 1], names=["chrom", "length"]
    )
    autosomes = set(df["chrom"].unique())
    genome_autosome_mb = fai.loc[fai["chrom"].isin(autosomes), "length"].sum() / 1e6

    per_ind = (
        df.groupby("individual")["callable_bp"].sum().div(1e6).rename("callable_mb")
    ).reset_index()
    per_ind["species"] = sp
    per_ind["genome_autosome_mb"] = genome_autosome_mb
    per_ind["callable_pct"] = 100.0 * per_ind["callable_mb"] / genome_autosome_mb
    per_ind = per_ind.sort_values(
        "individual", key=lambda s: s.map(family_sample_key)
    ).reset_index(drop=True)
    return per_ind, genome_autosome_mb


def main():
    frames = []
    totals = {}
    for sp in SPECIES:
        per_ind, total_mb = load_species(sp)
        frames.append(per_ind)
        totals[sp] = total_mb
    summary = pd.concat(frames, ignore_index=True)[
        ["species", "individual", "callable_mb", "genome_autosome_mb", "callable_pct"]
    ]

    summary_path = os.path.join(HERE, "callable_genome_summary.tsv")
    summary.to_csv(summary_path, sep="\t", index=False)

    # --- Figure: one panel per species (4x2 grid; last cell hidden) ---
    ncol = 2
    nrow = 4
    fig, axes = plt.subplots(nrow, ncol, figsize=(15, 18))
    axes = axes.flatten()

    for ax, sp in zip(axes, SPECIES):
        sub = summary[summary["species"] == sp].reset_index(drop=True)
        total_mb = totals[sp]
        x = range(len(sub))
        ax.bar(x, sub["callable_mb"], color=BAR_COLOR, width=0.8)
        ax.axhline(total_mb, color=REF_COLOR, linestyle="--", linewidth=1.2)
        ax.set_xticks(list(x))
        ax.set_xticklabels(sub["individual"], rotation=90, fontsize=6)
        ax.set_ylabel("Callable autosomal size (Mb)")
        ax.set_ylim(0, total_mb * 1.08)
        ax.set_title(
            f"{sp}  (n={len(sub)}, mean callable = {sub['callable_pct'].mean():.1f}%)",
            fontsize=11,
        )
        ax.text(
            0.985,
            0.965,
            f"autosomal genome: {total_mb:,.0f} Mb",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=8,
            color=REF_COLOR,
        )

        # Secondary axis: same data as % of autosomal genome (linear rescale).
        sec = ax.secondary_yaxis(
            "right",
            functions=(lambda v, t=total_mb: v / t * 100.0,
                       lambda p, t=total_mb: p / 100.0 * t),
        )
        sec.set_ylabel("% of autosomal genome")

    # Hide any unused panels (8th cell for 7 species).
    for ax in axes[len(SPECIES):]:
        ax.set_visible(False)

    fig.suptitle(
        "Callable autosomal genome size per individual, by species (minDP26)",
        fontsize=15,
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.99))

    png_path = os.path.join(HERE, "callable_genome_size.png")
    pdf_path = os.path.join(HERE, "callable_genome_size.pdf")
    fig.savefig(png_path, dpi=150)
    fig.savefig(pdf_path)
    plt.close(fig)

    print(f"wrote {summary_path}")
    print(f"wrote {png_path}")
    print(f"wrote {pdf_path}")
    print()
    print("individuals per species:")
    print(summary.groupby("species").size().to_string())


if __name__ == "__main__":
    main()
