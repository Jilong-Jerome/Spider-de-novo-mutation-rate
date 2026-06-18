#!/usr/bin/env python3
"""Combine three pathway-bias bar panels into a vertically stacked figure."""
import argparse
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--panels", required=True,
                   help="JSON list of panel specs; each has tsv, label, "
                        "positive_label, negative_label, positive_color, "
                        "negative_color, unbiased_color")
    p.add_argument("--pathway-order", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--width", type=float, default=12.0)
    p.add_argument("--height", type=float, default=8.0)
    p.add_argument("--inconsistent-alpha", type=float, default=0.2)
    p.add_argument("--output-suffix", default="",
                   help="suffix appended before output extensions, e.g. _is_sig")
    return p.parse_args()


def render_panel(ax, counts, pathways, display_pathways, panel, show_xticks,
                 show_ylabel, legend_mode, inconsistent_alpha):
    bias_order = ["positive", "unbiased", "negative"]
    color_map = {
        "positive": panel["positive_color"],
        "unbiased": panel["unbiased_color"],
        "negative": panel["negative_color"],
    }
    label_map = {
        "positive": f"{panel['positive_label']}-biased",
        "unbiased": "unbiased",
        "negative": f"{panel['negative_label']}-biased",
    }

    def get(pw, bias, consistent):
        m = (counts["pathway"] == pw) & (counts["bias"] == bias) & (counts["consistent"] == consistent)
        return int(counts.loc[m, "n"].sum())

    n_path = len(pathways)
    n_bias = len(bias_order)
    bar_w = 0.8 / n_bias
    x_base = np.arange(n_path)

    for i, bias in enumerate(bias_order):
        offsets = x_base + (i - (n_bias - 1) / 2) * bar_w
        consistent_vals = np.array([get(pw, bias, True) for pw in pathways])
        inconsistent_vals = np.array([get(pw, bias, False) for pw in pathways])
        ax.bar(offsets, consistent_vals, width=bar_w,
               color=color_map[bias], alpha=1.0,
               edgecolor="black", linewidth=0.5)
        if bias != "unbiased":
            ax.bar(offsets, inconsistent_vals, width=bar_w,
                   bottom=consistent_vals,
                   color=color_map[bias], alpha=inconsistent_alpha,
                   edgecolor="black", linewidth=0.5)

    ax.margins(y=0.18)
    ax.set_xticks(x_base)
    if show_xticks:
        ax.set_xticklabels(display_pathways, rotation=20, ha="right")
        ax.set_xlabel("DNA repair pathways")
    else:
        ax.set_xticklabels([])
    if show_ylabel:
        ax.set_ylabel("# OGs")
    else:
        ax.set_ylabel("")
    if panel.get("title"):
        ax.text(0.01, 0.97, panel["title"], transform=ax.transAxes,
                ha="left", va="top", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if legend_mode == "top":
        handles = [
            mpatches.Patch(facecolor=color_map[b], edgecolor="black", label=label_map[b])
            for b in bias_order
        ]
        ax.legend(handles=handles, title="Bias direction",
                  loc="upper left", bbox_to_anchor=(1.01, 1.0),
                  frameon=False)


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    plt.rcParams.update({
        "font.size": 13,
        "axes.labelsize": 14,
        "axes.titlesize": 14,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 11,
        "legend.title_fontsize": 12,
    })

    panels = json.loads(args.panels)
    pathways = [p.strip() for p in args.pathway_order.split(",") if p.strip()]
    display_pathways = ["FA" if p == "Fanconi_anemia" else p for p in pathways]

    fig, axes = plt.subplots(
        nrows=len(panels), ncols=1,
        figsize=(args.width, args.height),
        sharex=True,
    )
    if len(panels) == 1:
        axes = [axes]

    n = len(panels)
    middle = n // 2
    bottom_panel_for_shared_legend = None
    for i, (ax, panel) in enumerate(zip(axes, panels)):
        counts = pd.read_csv(panel["tsv"], sep="\t")
        legend_mode = "top" if i == 0 else "none"
        render_panel(ax, counts, pathways, display_pathways, panel,
                     show_xticks=(i == n - 1),
                     show_ylabel=(i == middle),
                     legend_mode=legend_mode,
                     inconsistent_alpha=args.inconsistent_alpha)
        if i == 1 and n >= 3:
            bottom_panel_for_shared_legend = panel

    fig.tight_layout()

    if bottom_panel_for_shared_legend is not None:
        panel = bottom_panel_for_shared_legend
        bias_order = ["positive", "unbiased", "negative"]
        color_map = {
            "positive": panel["positive_color"],
            "unbiased": panel["unbiased_color"],
            "negative": panel["negative_color"],
        }
        label_map = {
            "positive": f"{panel['positive_label']}-biased",
            "unbiased": "unbiased",
            "negative": f"{panel['negative_label']}-biased",
        }
        handles = [
            mpatches.Patch(facecolor=color_map[b], edgecolor="black", label=label_map[b])
            for b in bias_order
        ]
        bbox1 = axes[1].get_position()
        bbox2 = axes[2].get_position()
        legend_y = (bbox1.y0 + bbox2.y1) / 2.0 if bbox1.y0 > bbox2.y1 else (bbox1.y0 + bbox2.y0) / 2.0 + (bbox1.height + bbox2.height) / 4.0
        legend_y = ((bbox1.y0 + bbox1.y1) + (bbox2.y0 + bbox2.y1)) / 4.0
        legend_x = max(bbox1.x1, bbox2.x1) + 0.01
        fig.legend(handles=handles, title="Bias direction",
                   loc="center left", bbox_to_anchor=(legend_x, legend_y),
                   bbox_transform=fig.transFigure, frameon=False)
    out_pdf = os.path.join(args.out_dir, f"combined_pathway_bias_bar{args.output_suffix}.pdf")
    out_png = os.path.join(args.out_dir, f"combined_pathway_bias_bar{args.output_suffix}.png")
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"Wrote {out_pdf}")
    print(f"Wrote {out_png}")


if __name__ == "__main__":
    main()
