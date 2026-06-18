"""
Paternal bias analysis: unique vs shared de novo mutations (DNMs).

Species: BIC, LIN, TEN  (BIC_family2 and BIC_family3 excluded)
Output:  results/paternal_bias_unique_shared_summary.txt
"""

import re
import os
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
KEEP_SPECIES  = {"BIC", "LIN", "TEN"}
EXCL_FAMILIES = {"BIC_family2", "BIC_family3"}
DATA_FILE     = "data/autosomal_DNM_spread_sheet.tsv"
OUT_FILE      = "results/paternal_bias_unique_shared_summary.txt"

CHILD_RE = re.compile(r"^([A-Z]+)_(family\d+)_S\d+_offspring$")

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def parse_child(child_str):
    """Return (species, family_label) or (None, None) if not parseable."""
    m = CHILD_RE.match(str(child_str).strip())
    if m:
        return m.group(1), f"{m.group(1)}_{m.group(2)}"
    return None, None


def resolve_phase(phases):
    """
    Collapse a list of phase labels for a shared mutation.

    Priority:
      all same           → that phase
      P/U mix (no M)     → P
      M/U mix (no P)     → M
      P/M conflict       → CONFLICT
      all U              → U
      all F              → F
      F mixed with other → treat F as absent, re-evaluate remaining
    """
    # Drop F entries and check if anything is left
    non_f = [p for p in phases if p != "F"]
    if not non_f:
        return "F"

    unique = set(non_f)

    if len(unique) == 1:
        return unique.pop()

    has_p = "P" in unique
    has_m = "M" in unique

    if has_p and has_m:
        return "CONFLICT"
    if has_p:          # P + U
        return "P"
    if has_m:          # M + U
        return "M"
    # Only U left
    return "U"


def phase_breakdown(series):
    """Count P, M, U, F, CONFLICT in a phase Series."""
    counts = series.value_counts()
    return {
        "P":        int(counts.get("P",        0)),
        "M":        int(counts.get("M",        0)),
        "U":        int(counts.get("U",        0)),
        "F":        int(counts.get("F",        0)),
        "CONFLICT": int(counts.get("CONFLICT", 0)),
    }


def paternal_bias(series):
    """Return (n_P, n_M, bias_ratio) where bias = P/M; NaN if M==0."""
    phased = series[series.isin(["P", "M"])]
    n_p = int((phased == "P").sum())
    n_m = int((phased == "M").sum())
    ratio = n_p / n_m if n_m > 0 else float("nan")
    return n_p, n_m, ratio


def run_fisher(uP, uM, sP, sM):
    """2×2 Fisher exact test: [[uP, uM], [sP, sM]]. Returns (OR, p, table)."""
    table = [[uP, uM], [sP, sM]]
    if (uP + uM) == 0 or (sP + sM) == 0:
        return float("nan"), float("nan"), table
    odds_ratio, p_value = fisher_exact(table, alternative="two-sided")
    return float(odds_ratio), float(p_value), table


def fmt_ratio(r):
    return f"{r:.4f}" if not (isinstance(r, float) and np.isnan(r)) else "NaN"


def format_fixed_table(rows, headers):
    """
    Render a list-of-lists as a plain-text fixed-width table.
    rows    : list of lists (all values will be str-cast)
    headers : list of column header strings
    """
    all_rows = [headers] + [[str(v) for v in row] for row in rows]
    col_w = [max(len(r[c]) for r in all_rows) for c in range(len(headers))]
    sep = "  ".join("-" * w for w in col_w)
    lines = []
    lines.append("  ".join(h.ljust(col_w[i]) for i, h in enumerate(headers)))
    lines.append(sep)
    for row in all_rows[1:]:
        lines.append("  ".join(str(v).ljust(col_w[i]) for i, v in enumerate(row)))
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    lines = []  # accumulate report sections

    def section(title):
        lines.append("")
        lines.append("=" * 70)
        lines.append(title)
        lines.append("=" * 70)

    def log(msg=""):
        lines.append(msg)

    # ------------------------------------------------------------------
    # 1. Load data
    # ------------------------------------------------------------------
    df_raw = pd.read_csv(DATA_FILE, sep="\t", dtype=str)
    n_raw = len(df_raw)

    section("FILTERING SUMMARY")
    log(f"Raw rows loaded               : {n_raw}")

    # ------------------------------------------------------------------
    # 2. Parse child column
    # ------------------------------------------------------------------
    parsed = df_raw["child"].apply(parse_child)
    df_raw["species"]       = parsed.apply(lambda x: x[0])
    df_raw["family_label"]  = parsed.apply(lambda x: x[1])

    n_no_parse = df_raw["species"].isna().sum()
    log(f"Rows with unparseable child   : {n_no_parse}  (excluded)")
    df_raw = df_raw[df_raw["species"].notna()].copy()

    # ------------------------------------------------------------------
    # 3. Filter species
    # ------------------------------------------------------------------
    mask_species = df_raw["species"].isin(KEEP_SPECIES)
    n_excl_species = (~mask_species).sum()
    log(f"Rows excluded (wrong species) : {n_excl_species}  "
        f"(kept: {', '.join(sorted(KEEP_SPECIES))})")
    df = df_raw[mask_species].copy()

    # ------------------------------------------------------------------
    # 4. Exclude BIC_family2 / BIC_family3
    # ------------------------------------------------------------------
    mask_family = df["family_label"].isin(EXCL_FAMILIES)
    n_excl_fam  = mask_family.sum()
    excl_list   = sorted(df.loc[mask_family, "family_label"].unique())
    log(f"Rows excluded (excluded fams) : {n_excl_fam}  "
        f"({', '.join(excl_list)})")
    df = df[~mask_family].copy()

    log(f"Rows retained for analysis    : {len(df)}")

    families_kept = sorted(df["family_label"].unique())
    log(f"Families retained             : {', '.join(families_kept)}")

    # ------------------------------------------------------------------
    # 5. Classify unique vs shared
    # ------------------------------------------------------------------
    event_key = ["chrom", "pos", "ref", "alt", "family_label"]

    grp = (
        df.groupby(event_key)
          .agg(
              species      = ("species",      "first"),
              n_offspring  = ("child",        "nunique"),
              phases       = ("Phased",       list),
          )
          .reset_index()
    )

    grp["mutation_type"]   = grp["n_offspring"].apply(
        lambda n: "shared" if n > 1 else "unique"
    )
    grp["resolved_phase"]  = grp["phases"].apply(resolve_phase)

    n_unique_events = (grp["mutation_type"] == "unique").sum()
    n_shared_events = (grp["mutation_type"] == "shared").sum()

    # ------------------------------------------------------------------
    # 6. Mutation counts per species
    # ------------------------------------------------------------------
    section("MUTATION COUNTS (collapsed events)")

    count_rows = []
    species_list = sorted(grp["species"].unique())

    for sp in species_list:
        sub = grp[grp["species"] == sp]
        nu  = (sub["mutation_type"] == "unique").sum()
        ns  = (sub["mutation_type"] == "shared").sum()
        count_rows.append([sp, nu, ns, nu + ns])

    # Overall
    count_rows.append(["OVERALL",
                        n_unique_events, n_shared_events,
                        n_unique_events + n_shared_events])

    log(format_fixed_table(
        count_rows,
        ["Species", "Unique", "Shared", "Total"]
    ))

    # ------------------------------------------------------------------
    # 7. Phasing breakdown
    # ------------------------------------------------------------------
    section("PHASING BREAKDOWN (per species and overall)")

    pb_headers = ["Group", "Class", "P", "M", "U", "F", "CONFLICT", "Total_events"]
    pb_rows = []

    for sp in species_list:
        for cls in ["unique", "shared"]:
            sub   = grp[(grp["species"] == sp) & (grp["mutation_type"] == cls)]
            bd    = phase_breakdown(sub["resolved_phase"])
            total = len(sub)
            pb_rows.append([sp, cls,
                             bd["P"], bd["M"], bd["U"], bd["F"], bd["CONFLICT"],
                             total])

    for cls in ["unique", "shared"]:
        sub   = grp[grp["mutation_type"] == cls]
        bd    = phase_breakdown(sub["resolved_phase"])
        total = len(sub)
        pb_rows.append(["OVERALL", cls,
                         bd["P"], bd["M"], bd["U"], bd["F"], bd["CONFLICT"],
                         total])

    log(format_fixed_table(pb_rows, pb_headers))

    # Check P+M+U+F+CONFLICT == Total_events
    log("")
    log("Integrity check (P+M+U+F+CONFLICT == Total_events):")
    all_ok = True
    for row in pb_rows:
        _sp, _cls, p, m, u, f, c, tot = row
        computed = p + m + u + f + c
        ok = "OK" if computed == tot else f"MISMATCH (got {computed})"
        if computed != tot:
            all_ok = False
        log(f"  {_sp:10s}  {_cls:6s}  {computed} == {tot}  {ok}")
    if all_ok:
        log("  All checks passed.")

    # ------------------------------------------------------------------
    # 8. Paternal bias table
    # ------------------------------------------------------------------
    section("PATERNAL BIAS (P / M)")

    bias_headers = ["Group", "Class", "P", "M", "P/M_ratio"]
    bias_rows    = []

    for sp in species_list:
        for cls in ["unique", "shared"]:
            sub        = grp[(grp["species"] == sp) & (grp["mutation_type"] == cls)]
            n_p, n_m, r = paternal_bias(sub["resolved_phase"])
            bias_rows.append([sp, cls, n_p, n_m, fmt_ratio(r)])

    for cls in ["unique", "shared"]:
        sub        = grp[grp["mutation_type"] == cls]
        n_p, n_m, r = paternal_bias(sub["resolved_phase"])
        bias_rows.append(["OVERALL", cls, n_p, n_m, fmt_ratio(r)])

    log(format_fixed_table(bias_rows, bias_headers))

    # ------------------------------------------------------------------
    # 9. Fisher's exact test (overall)
    # ------------------------------------------------------------------
    section("FISHER'S EXACT TEST (overall, all species combined)")

    ov_unique = grp[grp["mutation_type"] == "unique"]
    ov_shared = grp[grp["mutation_type"] == "shared"]

    uP, uM, _ = paternal_bias(ov_unique["resolved_phase"])
    sP, sM, _ = paternal_bias(ov_shared["resolved_phase"])

    OR, pval, table = run_fisher(uP, uM, sP, sM)

    log("")
    log("2×2 contingency table:")
    log(f"                  Paternal (P)  Maternal (M)")
    log(f"  Unique         {uP:>12d}  {uM:>12d}")
    log(f"  Shared         {sP:>12d}  {sM:>12d}")
    log("")

    if np.isnan(pval):
        log("Fisher test skipped: one or both rows have 0 phased mutations.")
    else:
        log(f"Odds ratio : {OR:.4f}")
        log(f"P-value    : {pval:.4g}")
        sig = "significant (α=0.05)" if pval < 0.05 else "not significant (α=0.05)"
        log(f"Result     : {sig}")

    # ------------------------------------------------------------------
    # 10. CONFLICT list
    # ------------------------------------------------------------------
    conflicts = grp[grp["resolved_phase"] == "CONFLICT"]
    section("PHASING CONFLICTS (P/M disagreement within shared event)")
    if len(conflicts) == 0:
        log("No phasing conflicts detected.")
    else:
        log(f"{len(conflicts)} conflict(s) found (excluded from bias calculation):")
        log("")
        for _, row in conflicts.iterrows():
            log(f"  {row['species']}  {row['family_label']}  "
                f"{row['chrom']}:{row['pos']}  {row['ref']}>{row['alt']}  "
                f"phases={row['phases']}")

    # ------------------------------------------------------------------
    # 11. Notes
    # ------------------------------------------------------------------
    section("NOTES")
    log("1. Species analysed : BIC, LIN, TEN")
    log("2. Excluded families: BIC_family2, BIC_family3")
    log("3. Shared mutations collapsed to one event per (chrom, pos, ref, alt, family).")
    log("4. Phase resolution for shared events (priority order):")
    log("     all same        → keep that phase")
    log("     P + U (no M)    → resolve to P")
    log("     M + U (no P)    → resolve to M")
    log("     P + M           → CONFLICT (excluded from bias/test)")
    log("     all U           → U")
    log("     all F           → F")
    log("     F mixed others  → F labels dropped; remaining labels re-evaluated")
    log("5. Paternal bias = P / M; U, F, CONFLICT excluded.")
    log("6. Fisher's exact test run once (overall), two-sided.")
    log("   Test skipped (NaN) if either class has 0 phased mutations.")
    log(f"7. Analysis run: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # ------------------------------------------------------------------
    # Write output
    # ------------------------------------------------------------------
    os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)
    report = "\n".join(lines).lstrip("\n")
    with open(OUT_FILE, "w") as fh:
        fh.write(report + "\n")

    print(f"Report written to: {OUT_FILE}")


if __name__ == "__main__":
    main()
