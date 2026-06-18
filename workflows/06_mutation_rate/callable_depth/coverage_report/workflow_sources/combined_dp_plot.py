#!/usr/bin/env python3
"""
Combined cross-species barplot of per-individual autosomal DP (all sites).

Reads each species' `{SP}_autosome_dp_per_individual.tsv`
(columns sample/n_sites/mean_DP_all/median_DP_all), concatenates them with a
`species` column, writes a merged table, and draws a panelled barplot:

  - panels stacked one per row, grouped by species pair / outgroup
    (AFR-MIM, DUM-TEN, SAR-BIC, then LIN as the outgroup),
  - within each panel one bar per individual, coloured by species,
  - bar height = mean_DP_all, with median_DP_all overlaid as a dot,
  - a single global grand-mean line repeated in every panel = grand mean of all
    individuals' mean_DP_all across all species.
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Panel layout: species pairs each on their own row, LIN outgroup last.
GROUPS = [
    ("AFR–MIM", ["AFR", "MIM"]),
    ("DUM–TEN", ["DUM", "TEN"]),
    ("SAR–BIC", ["SAR", "BIC"]),
    ("LIN (outgroup)", ["LIN"]),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", nargs="+", required=True,
                    help="SPECIES:tsv_path tokens, one per species")
    ap.add_argument("--merged-out", required=True)
    ap.add_argument("--plot-out", required=True)
    args = ap.parse_args()

    frames = []
    for token in args.inputs:
        sp, path = token.split(":", 1)
        df = pd.read_csv(path, sep="\t")
        df.insert(0, "species", sp)
        frames.append(df)
    merged = pd.concat(frames, ignore_index=True)
    # stable order: by species (input order), then by sample
    sp_order = []
    for token in args.inputs:
        sp = token.split(":", 1)[0]
        if sp not in sp_order:
            sp_order.append(sp)
    merged["species"] = pd.Categorical(merged["species"], categories=sp_order, ordered=True)
    merged = merged.sort_values(["species", "sample"]).reset_index(drop=True)
    merged.to_csv(args.merged_out, sep="\t", index=False)

    # grand mean across all individuals of all species (unweighted over individuals)
    grand_mean = merged["mean_DP_all"].mean()

    # one colour per species (stable across panels)
    cmap = plt.get_cmap("tab10")
    color_of = {sp: cmap(i % 10) for i, sp in enumerate(sp_order)}

    # Resolve panels: keep only species present in the data; append any species
    # not listed in GROUPS as a trailing fallback panel so nothing is dropped.
    present = set(sp_order)
    panels = []
    grouped = set()
    for title, members in GROUPS:
        members = [sp for sp in members if sp in present]
        if members:
            panels.append((title, members))
            grouped.update(members)
    leftover = [sp for sp in sp_order if sp not in grouped]
    if leftover:
        panels.append(("other", leftover))

    # width sized to the largest panel so bar widths are consistent across rows
    max_n = max(
        int((merged["species"].isin(members)).sum())
        for _, members in panels
    )
    fig_w = max(8.0, max_n * 0.22)
    fig, axes = plt.subplots(
        nrows=len(panels), ncols=1, sharey=True,
        figsize=(fig_w, 2.2 * len(panels)),
    )
    if len(panels) == 1:
        axes = [axes]

    for ax, (title, members) in zip(axes, panels):
        sub = merged[merged["species"].isin(members)].copy()
        sub["species"] = pd.Categorical(sub["species"], categories=members, ordered=True)
        sub = sub.sort_values(["species", "sample"]).reset_index(drop=True)
        n = len(sub)
        x = np.arange(n)
        bar_colors = [color_of[sp] for sp in sub["species"]]

        ax.bar(x, sub["mean_DP_all"], color=bar_colors, width=0.85)
        ax.scatter(x, sub["median_DP_all"], color="black", s=10, zorder=3)
        ax.axhline(grand_mean, color="red", linestyle="--", linewidth=1.2)

        # species block tick labels at the centre of each species' bars
        centers, labels, boundaries = [], [], []
        for sp in members:
            idx = np.where(sub["species"].to_numpy() == sp)[0]
            if len(idx) == 0:
                continue
            centers.append(idx.mean())
            labels.append(f"{sp}\n(n={len(idx)})")
            boundaries.append(idx.max() + 0.5)
        for b in boundaries[:-1]:
            ax.axvline(b, color="grey", linewidth=0.4, alpha=0.5)
        ax.set_xticks(centers)
        ax.set_xticklabels(labels)
        ax.set_xlim(-0.5, n - 0.5)
        ax.set_ylabel("Autosomal DP")
        ax.set_title(title, fontsize=10, loc="left")

    # figure-level title + one shared legend
    fig.suptitle("Per-individual mean autosomal sequencing depth (VCF DP, all sites)")
    handles = [plt.Rectangle((0, 0), 1, 1, color=color_of[sp]) for sp in sp_order]
    handles.append(plt.Line2D([0], [0], marker="o", color="black", linestyle="None",
                              markersize=5))
    handles.append(plt.Line2D([0], [0], color="red", linestyle="--"))
    leg_labels = list(sp_order) + ["median DP", f"grand mean = {grand_mean:.2f}"]
    fig.legend(handles, leg_labels, ncol=len(leg_labels), fontsize=8,
               loc="lower center", bbox_to_anchor=(0.5, 0.0))

    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    fig.savefig(args.plot_out, dpi=150)
    print(f"Wrote {args.merged_out} ({len(merged)} individuals) and {args.plot_out}")


if __name__ == "__main__":
    main()
