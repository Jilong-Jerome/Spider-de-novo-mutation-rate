#!/usr/bin/env python3
"""
plot_callable_segments.py

Combined callable genome segment + de novo mutation visualisation.
One horizontal panel per autosome (shared x-axis, scaled to the longest chromosome).
Within each panel:
  - A grey bar spanning the true chromosome length
  - One thin row per individual showing callable BED segments (steelblue)
  - De novo mutations for each individual overlaid as crimson dots

Callable segments are bridged at visualisation time with --bridge_gap (default
20 kb): nearby segments are merged into one bar.  The alpha transparency of each
bar reflects the fraction of the bridged span that is genuinely callable.

Input BED files: per-offspring TSVs from compute_callable_bed.py
    columns: offspring  chrom  start  end
Input DNM file: species-level classified DNM TSV
    columns include: chrom  pos  child  father  mother  ref  alt  ...
"""
import argparse
import glob
import os
import re
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def load_fai(fai_path):
    """Return {chrom: length} dict."""
    chrom_lengths = {}
    with open(fai_path) as fh:
        for line in fh:
            parts = line.strip().split('\t')
            chrom_lengths[parts[0]] = int(parts[1])
    return chrom_lengths


def chrom_sort_key(chrom):
    """Sort chromosomes numerically by trailing integer."""
    m = re.search(r'(\d+)$', chrom)
    return int(m.group(1)) if m else 0


def bridge_segments(starts, ends, bridge_gap):
    """
    Merge BED segments whose inter-segment gap is <= bridge_gap bp.

    Parameters
    ----------
    starts, ends : array-like of int (0-based half-open, sorted by start)
    bridge_gap   : int, maximum gap to bridge (bp)

    Returns
    -------
    list of (bridged_start, bridged_end, callable_fraction)
        callable_fraction = sum of original callable bp / bridged span
    """
    if len(starts) == 0:
        return []

    result = []
    seg_start   = int(starts[0])
    seg_end     = int(ends[0])
    callable_bp = seg_end - seg_start

    for s, e in zip(starts[1:], ends[1:]):
        s, e = int(s), int(e)
        if s - seg_end <= bridge_gap:
            callable_bp += e - s
            seg_end      = e
        else:
            span = seg_end - seg_start
            result.append((seg_start, seg_end, callable_bp / span))
            seg_start   = s
            seg_end     = e
            callable_bp = e - s

    span = seg_end - seg_start
    result.append((seg_start, seg_end, callable_bp / span))
    return result


def main():
    parser = argparse.ArgumentParser(
        description='Plot callable segments and de novo mutations per individual.'
    )
    parser.add_argument('--bed_dir',        required=True,
                        help='Directory containing per-offspring callable BED TSV files')
    parser.add_argument('--species_prefix', required=True,
                        help='Species prefix used in filenames (e.g. afr)')
    parser.add_argument('--genome_fai',     required=True)
    parser.add_argument('--dnm_file',       required=True,
                        help='DNM TSV (columns: chrom, pos, child, father, mother, '
                             'ref, alt, ...)')
    parser.add_argument('--x_chroms',       nargs='+', default=[])
    parser.add_argument('--bridge_gap',     type=int, default=20000,
                        help='Bridge callable segments within this many bp '
                             '(default: 20000). Alpha encodes callable fraction.')
    parser.add_argument('--output',         required=True,
                        help='Output PDF path')
    args = parser.parse_args()

    x_chroms_set  = set(args.x_chroms)
    chrom_lengths = load_fai(args.genome_fai)

    # ------------------------------------------------------------------
    # Load callable BED files
    # ------------------------------------------------------------------
    pattern   = os.path.join(args.bed_dir, f"{args.species_prefix}_*_callable.bed")
    bed_files = glob.glob(pattern)
    if not bed_files:
        raise FileNotFoundError(f"No callable BED files found matching: {pattern}")

    bed_df = pd.concat(
        [pd.read_csv(f, sep='\t') for f in bed_files],
        ignore_index=True,
    )
    bed_df = bed_df[~bed_df['chrom'].isin(x_chroms_set)].copy()
    bed_df = bed_df.sort_values(['offspring', 'chrom', 'start']).reset_index(drop=True)

    # ------------------------------------------------------------------
    # Load DNMs — per individual, no deduplication (each child's own events)
    # ------------------------------------------------------------------
    dnm_df = pd.read_csv(args.dnm_file, sep='\t')
    dnm_df = dnm_df[~dnm_df['chrom'].isin(x_chroms_set)].copy()

    # ------------------------------------------------------------------
    # Plot dimensions
    # ------------------------------------------------------------------
    autosomes = sorted(
        [c for c in chrom_lengths if c not in x_chroms_set],
        key=chrom_sort_key,
    )
    n_chroms      = len(autosomes)
    max_chrom_len = max(chrom_lengths[c] for c in autosomes)

    offspring_list = sorted(bed_df['offspring'].unique())
    n_offspring    = len(offspring_list)
    offspring_idx  = {name: i for i, name in enumerate(offspring_list)}

    # Layout constants (axis data units)
    CHROM_Y     = 0
    CHROM_H     = 0.6
    IND_H       = 0.75
    IND_SPACING = 0.05
    IND_UNIT    = IND_H + IND_SPACING

    def ind_y_center(i):
        return CHROM_Y - (CHROM_H / 2) - IND_SPACING - IND_H / 2 - i * IND_UNIT

    y_bottom = ind_y_center(n_offspring - 1) - IND_H / 2 - 0.1
    y_top    = CHROM_Y + CHROM_H / 2 + 0.2

    row_inch   = 0.10
    panel_inch = (n_offspring + 2) * row_inch
    fig_height = max(6, n_chroms * panel_inch + 1.0)
    fig_width  = 14

    fig, axes = plt.subplots(
        n_chroms, 1,
        figsize=(fig_width, fig_height),
        sharex=True,
    )
    if n_chroms == 1:
        axes = [axes]

    fig.suptitle(
        f"{args.species_prefix.upper()} — callable genome & de novo mutations "
        f"({n_offspring} individuals, bridge gap {args.bridge_gap // 1000} kb, "
        f"alpha = callable fraction)",
        fontsize=11, fontweight='bold',
    )

    # ------------------------------------------------------------------
    # Pre-set axes limits on all panels, then render the canvas once so
    # that the data→display transform is accurate.  This lets us convert
    # IND_H (data units) to points and compute a DNM marker size that is
    # slightly larger than the individual track height regardless of the
    # figure size or number of individuals.
    # ------------------------------------------------------------------
    for ax in axes:
        ax.set_xlim(0, max_chrom_len)
        ax.set_ylim(y_bottom, y_top)

    fig.canvas.draw()

    # IND_H in display pixels (transform from the first panel)
    ax0       = axes[0]
    px_lo     = ax0.transData.transform((0, 0))[1]
    px_hi     = ax0.transData.transform((0, IND_H))[1]
    ind_h_pts = abs(px_hi - px_lo) * 72 / fig.dpi   # pixels → points

    # Marker diameter = 1.5 × track height (slightly oversize); s = π*(d/2)²
    dnm_marker_s = 3.14159 * (ind_h_pts * 1.5 / 2) ** 2

    # ------------------------------------------------------------------
    # Draw each chromosome panel
    # ------------------------------------------------------------------
    for ax, chrom in zip(axes, autosomes):
        chrom_len  = chrom_lengths[chrom]
        chrom_bed  = bed_df[bed_df['chrom'] == chrom]
        chrom_dnms = dnm_df[dnm_df['chrom'] == chrom]

        # Chromosome bar (grey, true length)
        ax.broken_barh(
            [(0, chrom_len)],
            (CHROM_Y - CHROM_H / 2, CHROM_H),
            facecolors='lightgrey', edgecolors='none',
        )

        for ind_name, ind_i in offspring_idx.items():
            y_cen   = ind_y_center(ind_i)
            y_range = (y_cen - IND_H / 2, IND_H)

            # --- Callable segments (steelblue, alpha = callable fraction) ---
            ind_bed = chrom_bed[chrom_bed['offspring'] == ind_name]
            if not ind_bed.empty:
                bridged = bridge_segments(
                    ind_bed['start'].values,
                    ind_bed['end'].values,
                    args.bridge_gap,
                )
                for bstart, bend, callable_frac in bridged:
                    ax.broken_barh(
                        [(bstart, bend - bstart)],
                        y_range,
                        facecolors='steelblue',
                        edgecolors='none',
                        alpha=max(0.05, callable_frac),
                    )

            # --- De novo mutations (crimson dots, calibrated to track height) ---
            ind_dnms = chrom_dnms[chrom_dnms['child'] == ind_name]['pos'].values
            if len(ind_dnms) > 0:
                ax.scatter(
                    ind_dnms,
                    [y_cen] * len(ind_dnms),
                    s=dnm_marker_s,
                    color='crimson',
                    linewidths=0,
                    zorder=3,
                )

        ax.set_ylabel(chrom, rotation=0, ha='right', va='center', fontsize=8)

        ytick_positions = [ind_y_center(i) for i in range(n_offspring)]
        ax.set_yticks(ytick_positions)
        ax.set_yticklabels(offspring_list, fontsize=4)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.tick_params(axis='x', labelsize=7)

    # Shared x-axis Mb labels on bottom panel only
    axes[-1].xaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, _: f"{x/1e6:.0f} Mb")
    )

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig(args.output, bbox_inches='tight')
    plt.close()
    print(f"Saved: {args.output}")


if __name__ == '__main__':
    main()
