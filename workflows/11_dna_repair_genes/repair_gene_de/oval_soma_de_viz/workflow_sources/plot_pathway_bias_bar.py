#!/usr/bin/env python3
"""Per-pathway grouped bar plot of orthogroup bias direction.

For each pathway: three bars (positive bias / unbiased / negative bias).
Each directional bar is split by consistency across species:
  - solid (alpha=1.0): n_consistent == n_species  (consistently biased)
  - faded (alpha=0.4):  not consistent across species
The unbiased bar is always rendered solid.
"""
import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--in", dest="input_csv", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--comparison", required=True)
    p.add_argument("--positive-label", required=True)
    p.add_argument("--negative-label", required=True)
    p.add_argument("--log2fc-column", default="log2fc",
                   help="column holding the log2 fold change used for classification")
    p.add_argument("--log2fc-threshold", type=float, default=1.0,
                   help="|log2FC| > threshold counts as biased")
    p.add_argument("--count-mode", default="log2fc_threshold",
                   choices=["log2fc_threshold", "is_sig"],
                   help="how to assign rows to positive/unbiased/negative bins")
    p.add_argument("--is-sig-column", default="is_sig",
                   help="column holding TRUE/FALSE tested significance values")
    p.add_argument("--output-suffix", default="",
                   help="suffix appended before output extensions, e.g. _is_sig")
    p.add_argument("--positive-color", required=True)
    p.add_argument("--negative-color", required=True)
    p.add_argument("--unbiased-color", required=True)
    p.add_argument("--pathway-order", required=True,
                   help="comma-separated pathway order")
    p.add_argument("--consistency-column", default=None,
                   help="if set, rows are marked consistent when this column equals TRUE; "
                        "otherwise consistent = (n_consistent == n_species)")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    df = pd.read_csv(args.input_csv)
    lfc = pd.to_numeric(df[args.log2fc_column], errors="coerce")
    if args.count_mode == "log2fc_threshold":
        thr = args.log2fc_threshold
        df["bias"] = np.where(
            lfc > thr, "positive",
            np.where(lfc < -thr, "negative", "unbiased"),
        )
    else:
        sig = df[args.is_sig_column]
        if sig.dtype == object:
            sig = sig.astype(str).str.upper().eq("TRUE")
        else:
            sig = sig.astype(bool)
        df["bias"] = np.where(
            sig & (lfc > 0), "positive",
            np.where(sig & (lfc < 0), "negative", "unbiased"),
        )

    if args.consistency_column:
        col = df[args.consistency_column]
        if col.dtype == object:
            col = col.astype(str).str.upper().eq("TRUE")
        else:
            col = col.astype(bool)
        df["consistent"] = col & (df["bias"] != "unbiased")
    else:
        df["consistent"] = (df["n_consistent"] == df["n_species"]) & (df["bias"] != "unbiased")
    df.loc[df["bias"] == "unbiased", "consistent"] = True

    counts = (
        df.groupby(["pathway", "bias", "consistent"])
        .size()
        .rename("n")
        .reset_index()
    )
    counts.to_csv(
        os.path.join(args.out_dir, f"pathway_bias_counts{args.output_suffix}.tsv"),
        sep="\t",
        index=False,
    )

    pathways = [p.strip() for p in args.pathway_order.split(",") if p.strip()]
    observed = list(df["pathway"].unique())
    for p in observed:
        if p not in pathways:
            pathways.append(p)

    bias_order = ["positive", "unbiased", "negative"]
    color_map = {
        "positive": args.positive_color,
        "unbiased": args.unbiased_color,
        "negative": args.negative_color,
    }
    label_map = {
        "positive": f"{args.positive_label}-biased",
        "unbiased": "unbiased",
        "negative": f"{args.negative_label}-biased",
    }

    def get(pw, bias, consistent):
        m = (counts["pathway"] == pw) & (counts["bias"] == bias) & (counts["consistent"] == consistent)
        return int(counts.loc[m, "n"].sum())

    n_path = len(pathways)
    n_bias = len(bias_order)
    bar_w = 0.8 / n_bias
    x_base = np.arange(n_path)

    fig, ax = plt.subplots(figsize=(max(7, 1.4 * n_path), 5))

    for i, bias in enumerate(bias_order):
        offsets = x_base + (i - (n_bias - 1) / 2) * bar_w
        consistent_vals = np.array([get(pw, bias, True) for pw in pathways])
        inconsistent_vals = np.array([get(pw, bias, False) for pw in pathways])

        ax.bar(
            offsets, consistent_vals, width=bar_w,
            color=color_map[bias], alpha=1.0,
            edgecolor="black", linewidth=0.5,
        )
        if bias != "unbiased":
            ax.bar(
                offsets, inconsistent_vals, width=bar_w,
                bottom=consistent_vals,
                color=color_map[bias], alpha=0.4,
                edgecolor="black", linewidth=0.5,
            )

    display_pathways = ["FA" if p == "Fanconi_anemia" else p for p in pathways]
    ax.set_xticks(x_base)
    ax.set_xticklabels(display_pathways, rotation=20, ha="right")
    ax.set_xlabel("DNA repair pathways")
    ax.set_ylabel("Number of orthogroups")

    color_handles = [
        mpatches.Patch(facecolor=color_map[b], edgecolor="black", label=label_map[b])
        for b in bias_order
    ]
    ax.legend(handles=color_handles, title="Bias direction",
              loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()

    out_pdf = os.path.join(args.out_dir, f"pathway_bias_bar{args.output_suffix}.pdf")
    out_png = os.path.join(args.out_dir, f"pathway_bias_bar{args.output_suffix}.png")
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"Wrote {out_pdf}")
    print(f"Wrote {out_png}")


if __name__ == "__main__":
    main()
