#!/usr/bin/env python3
"""
combine_trio_visualizations.py

Step 11 of the dS-branch-trio pipeline. Joins the per-trio Step 10 outputs
into a single 2 x N panel figure:

    Row 1 (heatmaps): pairwise dS, shared colour scale across trios
    Row 2 (trees):    unrooted dS tree, shared branch-length scale

Inputs are passed as parallel `--trio_name / --trio_tree / --pairwise_summary
/ --pairwise_bootstrap` flags (each repeated once per trio, in matching order).
"""
import argparse
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from visualize_pairwise_dS import (
    compute_trio_data,
    draw_heatmap,
    draw_unrooted_tree,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trio_name", action="append", required=True)
    ap.add_argument("--trio_tree", action="append", required=True)
    ap.add_argument("--pairwise_summary", action="append", required=True)
    ap.add_argument("--pairwise_bootstrap", action="append", required=True)
    ap.add_argument("--output_dir", required=True)
    args = ap.parse_args()

    n = len(args.trio_name)
    if not (len(args.trio_tree) == len(args.pairwise_summary)
            == len(args.pairwise_bootstrap) == n):
        sys.exit("ERROR: --trio_name / --trio_tree / --pairwise_summary / "
                 "--pairwise_bootstrap must each be repeated the same number "
                 "of times.")

    os.makedirs(args.output_dir, exist_ok=True)

    trios = []
    for name, tree, summ, boot in zip(args.trio_name, args.trio_tree,
                                       args.pairwise_summary,
                                       args.pairwise_bootstrap):
        try:
            trios.append(compute_trio_data(summ, boot, tree, name))
        except ValueError as e:
            sys.exit(f"ERROR ({name}): {e}")

    # Shared colour scale on the heatmaps; trees use per-trio normalisation
    # so each tree fills its panel — the dS values themselves are still on
    # the annotations, so the comparison stays quantitative.
    vmax = max(np.nanmax(t["matrix"]) for t in trios)
    vmin = 0.0

    # Layout: 2 rows, N+1 columns (last col is the shared colorbar)
    fig = plt.figure(figsize=(5.0 * n + 1.4, 10.0))
    gs = fig.add_gridspec(
        2, n + 1,
        width_ratios=[*([1.0] * n), 0.06],
        height_ratios=[1.0, 1.0],
        hspace=0.05, wspace=0.30,
    )

    cbar_ax = fig.add_subplot(gs[0, n])

    for i, t in enumerate(trios):
        ax_h = fig.add_subplot(gs[0, i])
        draw_heatmap(
            ax_h,
            t["matrix"],
            t["species_order"],
            t["trio_name"],
            vmin=vmin, vmax=vmax,
            cbar=(i == 0),         # one shared cbar drawn from the first heatmap
            cbar_ax=cbar_ax if i == 0 else None,
            cbar_label="pairwise dS",
        )

        ax_t = fig.add_subplot(gs[1, i])
        draw_unrooted_tree(
            ax_t,
            species_order=t["species_order"],
            branch_lengths=[t["branch_point"][sp] for sp in t["species_order"]],
            ci_low=[t["ci"][sp][0] for sp in t["species_order"]],
            ci_high=[t["ci"][sp][1] for sp in t["species_order"]],
            trio_name=t["trio_name"],
            tree_scale=None,        # per-trio normalisation
            show_scale_bar=True,    # each tree carries its own scale bar
            title=False,
        )

    fig.suptitle("Pairwise dS across trios", fontsize=16, y=0.995)
    fig.subplots_adjust(left=0.05, right=0.94, top=0.95, bottom=0.04)

    out_prefix = os.path.join(args.output_dir, "combined_pairwise_dS")
    fig.savefig(f"{out_prefix}.png", dpi=200)
    fig.savefig(f"{out_prefix}.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
