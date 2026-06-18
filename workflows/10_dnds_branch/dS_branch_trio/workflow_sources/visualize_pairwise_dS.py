#!/usr/bin/env python3
"""
visualize_pairwise_dS.py

Step 10 of the dS-branch-trio pipeline. Consumes the pairwise dS summary +
bootstrap-replicate TSVs from Step 09 and produces, per trio:

  1. pairwise_dS_heatmap.{png,pdf} : 3x3 heatmap of pairwise dS.
  2. branch_dS_tree.{png,pdf}      : unrooted tree with dS branch lengths.
  3. branch_dS_lengths.tsv         : recalculated per-species dS branch
                                     lengths with 95% bootstrap CIs.

Branch lengths come from the standard three-point formula for an unrooted
3-leaf tree:
    L(A) = (d_AB + d_AC - d_BC) / 2
    L(B) = (d_AB + d_BC - d_AC) / 2
    L(C) = (d_AC + d_BC - d_AB) / 2

The plotting helpers (`draw_heatmap`, `draw_unrooted_tree`) and the data
loader (`compute_trio_data`) are also imported by the cross-trio combined
visualization step.
"""
import argparse
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


# ---------------------------------------------------------------------------
# Pure helpers
# ---------------------------------------------------------------------------

def parse_trio_tree(newick):
    """Return (sp1, sp2, outgroup) from a Newick string like '((SAR,PAC),BIC)'.

    Falls back to alphabetical order if the Newick cannot be parsed.
    """
    tokens = re.findall(r"[A-Za-z0-9_]+", newick)
    if len(tokens) == 3:
        return tokens[0], tokens[1], tokens[2]
    return tuple(sorted(tokens))


def pair_key(a, b):
    """Pair labels in the input TSVs are alphabetic, joined by underscore."""
    return "_".join(sorted([a, b]))


def build_distance_matrix(dS_by_pair, species_order):
    n = len(species_order)
    mat = np.zeros((n, n))
    for i, a in enumerate(species_order):
        for j, b in enumerate(species_order):
            if i == j:
                continue
            mat[i, j] = dS_by_pair[pair_key(a, b)]
    return mat


def three_point_branches(d_ab, d_ac, d_bc):
    """Branch lengths from leaves A, B, C to the central internal node."""
    return (
        (d_ab + d_ac - d_bc) / 2.0,
        (d_ab + d_bc - d_ac) / 2.0,
        (d_ac + d_bc - d_ab) / 2.0,
    )


def compute_trio_data(pairwise_summary_path, pairwise_bootstrap_path,
                      trio_tree, trio_name):
    """Read summary + bootstrap TSVs, return everything needed to plot a trio.

    Returns a dict with keys:
        trio_name, species_order, dS_point, matrix, branch_point, ci, n_bs
    """
    sp1, sp2, outgroup = parse_trio_tree(trio_tree)
    species_order = [sp1, sp2, outgroup]

    summary = pd.read_csv(pairwise_summary_path, sep="\t")
    dS_summary = summary[summary["metric"] == "dS"].set_index("pair")

    expected_pairs = [pair_key(sp1, sp2), pair_key(sp1, outgroup),
                      pair_key(sp2, outgroup)]
    missing = [p for p in expected_pairs if p not in dS_summary.index]
    if missing:
        raise ValueError(
            f"missing dS rows for pairs {missing} in {pairwise_summary_path}")

    dS_point = {p: float(dS_summary.loc[p, "estimate"]) for p in expected_pairs}
    matrix = build_distance_matrix(dS_point, species_order)

    d_ab = dS_point[pair_key(sp1, sp2)]
    d_ac = dS_point[pair_key(sp1, outgroup)]
    d_bc = dS_point[pair_key(sp2, outgroup)]
    L_sp1, L_sp2, L_out = three_point_branches(d_ab, d_ac, d_bc)
    branch_point = {sp1: L_sp1, sp2: L_sp2, outgroup: L_out}

    bs = pd.read_csv(pairwise_bootstrap_path, sep="\t")
    bs_wide = bs[["replicate", "pair", "dS"]].pivot(
        index="replicate", columns="pair", values="dS")

    n_bs = len(bs_wide)
    ci = {sp: (np.nan, np.nan) for sp in species_order}
    bs_branches = {sp: np.array([]) for sp in species_order}
    if n_bs > 0 and all(p in bs_wide.columns for p in expected_pairs):
        Lb_sp1, Lb_sp2, Lb_out = three_point_branches(
            bs_wide[pair_key(sp1, sp2)].to_numpy(),
            bs_wide[pair_key(sp1, outgroup)].to_numpy(),
            bs_wide[pair_key(sp2, outgroup)].to_numpy(),
        )
        bs_branches = {sp1: Lb_sp1, sp2: Lb_sp2, outgroup: Lb_out}
        for sp, vals in bs_branches.items():
            ci[sp] = (float(np.quantile(vals, 0.025)),
                      float(np.quantile(vals, 0.975)))

    return {
        "trio_name": trio_name,
        "species_order": species_order,
        "dS_point": dS_point,
        "matrix": matrix,
        "branch_point": branch_point,
        "ci": ci,
        "bs_branches": bs_branches,
        "n_bs": n_bs,
    }


# ---------------------------------------------------------------------------
# Social-vs-subsocial lineage shortening
# ---------------------------------------------------------------------------

def _summarise(values):
    """mean / median / 2.5% / 97.5% / n over a 1-D array, ignoring NaN."""
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return {"bs_mean": np.nan, "bs_median": np.nan,
                "bs_ci_low": np.nan, "bs_ci_high": np.nan, "n_bs": 0}
    return {
        "bs_mean": float(np.mean(arr)),
        "bs_median": float(np.median(arr)),
        "bs_ci_low": float(np.quantile(arr, 0.025)),
        "bs_ci_high": float(np.quantile(arr, 0.975)),
        "n_bs": int(arr.size),
    }


def compute_social_shortening(data, social_species):
    """Per-replicate social-vs-subsocial branch contrast.

    The two ingroup species are `data["species_order"][:2]`. Exactly one of
    them must be in `social_species` for the comparison to be defined.
    Returns a list of result rows (one per metric); when the contrast is
    undefined a single not_applicable row is returned instead.
    """
    sp1, sp2, _outgroup = data["species_order"]
    ingroup = [sp1, sp2]
    social = [s for s in ingroup if s in social_species]
    subsocial = [s for s in ingroup if s not in social_species]

    if len(social) != 1 or len(subsocial) != 1:
        reason = ("no_social_ingroup_species"
                  if len(social) == 0 else "two_social_ingroup_species")
        return [{
            "social_species": ",".join(social) if social else "NA",
            "subsocial_species": ",".join(subsocial) if subsocial else "NA",
            "metric": "NA",
            "estimate": np.nan,
            "bs_mean": np.nan, "bs_median": np.nan,
            "bs_ci_low": np.nan, "bs_ci_high": np.nan,
            "n_bs": 0,
            "status": f"not_applicable:{reason}",
        }]

    soc_sp = social[0]
    sub_sp = subsocial[0]
    L_soc_point = data["branch_point"][soc_sp]
    L_sub_point = data["branch_point"][sub_sp]
    L_soc_bs = np.asarray(data["bs_branches"][soc_sp], dtype=float)
    L_sub_bs = np.asarray(data["bs_branches"][sub_sp], dtype=float)

    # abs_diff = L_subsocial - L_social  (positive when social is shorter)
    abs_point = L_sub_point - L_soc_point
    abs_bs = L_sub_bs - L_soc_bs if L_soc_bs.size else np.array([])

    # rel_short = 1 - L_social / L_subsocial
    rel_point = (1.0 - L_soc_point / L_sub_point) if L_sub_point > 0 else np.nan
    if L_sub_bs.size:
        with np.errstate(divide="ignore", invalid="ignore"):
            rel_bs = np.where(L_sub_bs > 0, 1.0 - L_soc_bs / L_sub_bs, np.nan)
    else:
        rel_bs = np.array([])

    rows = []
    for metric, point, bs_vals in (
            ("abs_diff", abs_point, abs_bs),
            ("rel_short", rel_point, rel_bs)):
        stats = _summarise(bs_vals)
        rows.append({
            "social_species": soc_sp,
            "subsocial_species": sub_sp,
            "metric": metric,
            "estimate": float(point) if np.isfinite(point) else np.nan,
            **stats,
            "status": "ok",
        })
    return rows


# ---------------------------------------------------------------------------
# Reusable axes-based renderers
# ---------------------------------------------------------------------------

def draw_heatmap(ax, matrix, species_order, trio_name,
                 vmin=None, vmax=None, cbar=True, cbar_ax=None,
                 cbar_label="pairwise dS"):
    """Render a 3x3 dS heatmap onto a given matplotlib Axes."""
    df = pd.DataFrame(matrix, index=species_order, columns=species_order)
    cbar_kws = {"label": cbar_label} if cbar else None
    sns.heatmap(
        df,
        annot=True,
        fmt=".4f",
        cmap="viridis",
        vmin=vmin, vmax=vmax,
        cbar=cbar,
        cbar_ax=cbar_ax,
        cbar_kws=cbar_kws,
        square=True,
        linewidths=0.5,
        ax=ax,
    )
    ax.set_title(f"Pairwise dS — {trio_name}")


def draw_unrooted_tree(ax, species_order, branch_lengths, ci_low, ci_high,
                       trio_name, tree_scale=None, show_scale_bar=True,
                       title=True, leaf_fontsize=14, annot_fontsize=9):
    """Render an unrooted 3-leaf tree onto a given matplotlib Axes.

    Branch positions are in data coordinates (normalised so the longest branch
    spans 1.0 unit when `tree_scale` is None, or `length / tree_scale`
    otherwise). Leaf labels and branch annotations use display-point offsets
    via `ax.annotate`, so spacing stays consistent regardless of panel size.

    `tree_scale` sets the dS value that maps to normalised length 1.0; pass
    a shared value across panels for visually-comparable branch lengths.
    """
    angles_deg = [90.0, 210.0, 330.0]
    angles = np.deg2rad(angles_deg)

    if tree_scale is None or tree_scale <= 0:
        tree_scale = max(branch_lengths) if max(branch_lengths) > 0 else 1.0
    norm_lengths = [L / tree_scale for L in branch_lengths]

    ax.set_aspect("equal")
    ax.axis("off")

    leaf_pt = 16                  # tip → species label, in display points
    annot_pt = 22                 # species label → branch annotation, points

    for sp, L_data, L_norm, lo, hi, theta in zip(
            species_order, branch_lengths, norm_lengths,
            ci_low, ci_high, angles):
        x_tip = L_norm * np.cos(theta)
        y_tip = L_norm * np.sin(theta)
        ax.plot([0, x_tip], [0, y_tip], color="black", linewidth=2.5)

        # Species label, offset in display points along the branch direction
        ax.annotate(
            sp,
            xy=(x_tip, y_tip), xycoords="data",
            xytext=(leaf_pt * np.cos(theta), leaf_pt * np.sin(theta)),
            textcoords="offset points",
            ha="center", va="center",
            fontsize=leaf_fontsize, fontweight="bold",
        )

        # Branch annotation either ABOVE or BELOW the species label.
        # Pick the side opposite to where the branch enters the leaf:
        #   branch points upward to the leaf  → annotation above species
        #   branch points downward to the leaf → annotation below species
        side = 1.0 if np.sin(theta) >= 0 else -1.0
        if np.isnan(lo) or np.isnan(hi):
            label = f"dS = {L_data:.4f}"
        else:
            label = f"dS = {L_data:.4f}\n[{lo:.4f}, {hi:.4f}]"
        ax.annotate(
            label,
            xy=(x_tip, y_tip), xycoords="data",
            xytext=(leaf_pt * np.cos(theta),
                    leaf_pt * np.sin(theta) + side * annot_pt),
            textcoords="offset points",
            ha="center", va="bottom" if side > 0 else "top",
            fontsize=annot_fontsize, color="0.25",
        )

    ax.plot(0, 0, "o", color="black", markersize=5)

    if show_scale_bar:
        bar_norm = 0.5
        bar_data = bar_norm * tree_scale
        bar_y = -1.25
        ax.plot([-1.0, -1.0 + bar_norm], [bar_y, bar_y],
                color="black", linewidth=2)
        ax.annotate(
            f"dS = {bar_data:.4f}",
            xy=(-1.0 + bar_norm / 2.0, bar_y), xycoords="data",
            xytext=(0, -4), textcoords="offset points",
            ha="center", va="top", fontsize=annot_fontsize,
        )

    if title:
        ax.set_title(f"Unrooted dS tree — {trio_name}", fontsize=13, pad=10)

    # Data range is just the normalised branches plus a small margin; the
    # display-point offsets for labels naturally extend past the axes box.
    ax.set_xlim(-1.25, 1.25)
    ax.set_ylim(-1.40 if show_scale_bar else -1.25, 1.25)


# ---------------------------------------------------------------------------
# File-saving wrappers used by Step 10 (single trio)
# ---------------------------------------------------------------------------

def plot_heatmap(matrix, species_order, trio_name, out_prefix):
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    draw_heatmap(ax, matrix, species_order, trio_name)
    plt.tight_layout()
    fig.savefig(f"{out_prefix}.png", dpi=200)
    fig.savefig(f"{out_prefix}.pdf")
    plt.close(fig)


def plot_unrooted_tree(species_order, branch_lengths, ci_low, ci_high,
                       trio_name, out_prefix):
    fig, ax = plt.subplots(figsize=(8.0, 8.0))
    draw_unrooted_tree(ax, species_order, branch_lengths, ci_low, ci_high,
                       trio_name)
    fig.savefig(f"{out_prefix}.png", dpi=200, bbox_inches="tight")
    fig.savefig(f"{out_prefix}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairwise_summary", required=True)
    ap.add_argument("--pairwise_bootstrap", required=True)
    ap.add_argument("--trio_tree", required=True)
    ap.add_argument("--trio_name", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--social_species", default="SAR,DUM,MIM",
                    help="Comma-separated species names treated as social "
                         "for the lineage-shortening contrast.")
    args = ap.parse_args()
    social_species = {s.strip() for s in args.social_species.split(",")
                      if s.strip()}

    os.makedirs(args.output_dir, exist_ok=True)

    try:
        data = compute_trio_data(
            args.pairwise_summary, args.pairwise_bootstrap,
            args.trio_tree, args.trio_name,
        )
    except ValueError as e:
        sys.exit(f"ERROR: {e}")

    species_order = data["species_order"]
    branch_point = data["branch_point"]
    ci = data["ci"]

    plot_heatmap(
        data["matrix"], species_order, args.trio_name,
        os.path.join(args.output_dir, "pairwise_dS_heatmap"),
    )

    rows = []
    for sp in species_order:
        lo, hi = ci[sp]
        rows.append({
            "trio": args.trio_name,
            "species": sp,
            "branch_dS": branch_point[sp],
            "bs_ci_low": lo,
            "bs_ci_high": hi,
            "n_bs": data["n_bs"],
        })
    out_tsv = os.path.join(args.output_dir, "branch_dS_lengths.tsv")
    pd.DataFrame(rows).to_csv(out_tsv, sep="\t", index=False,
                              float_format="%.6f")

    plot_unrooted_tree(
        species_order=species_order,
        branch_lengths=[branch_point[sp] for sp in species_order],
        ci_low=[ci[sp][0] for sp in species_order],
        ci_high=[ci[sp][1] for sp in species_order],
        trio_name=args.trio_name,
        out_prefix=os.path.join(args.output_dir, "branch_dS_tree"),
    )

    shortening_rows = compute_social_shortening(data, social_species)
    shortening_tsv = os.path.join(args.output_dir,
                                  "social_lineage_shortening.tsv")
    shortening_cols = ["trio", "social_species", "subsocial_species",
                       "metric", "estimate", "bs_mean", "bs_median",
                       "bs_ci_low", "bs_ci_high", "n_bs", "status"]
    for row in shortening_rows:
        row["trio"] = args.trio_name
    pd.DataFrame(shortening_rows, columns=shortening_cols).to_csv(
        shortening_tsv, sep="\t", index=False, float_format="%.6f")


if __name__ == "__main__":
    main()
