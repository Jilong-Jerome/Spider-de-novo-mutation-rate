"""
Fraction of sibling-shared de novo mutations (DNMs) vs. sociality.

Tests whether the fraction of mutation loci that are sibling-shared differs
between social and subsocial species with a pooled 2x2 Fisher's exact test.

Sociality groups (from the study species tree / DNM_alpha config):
  social    : DUM, MIM, SAR
  subsocial : AFR, BIC, LIN, TEN

Locus definition mirrors DNM_alpha/workflow_sources/summarize_shared_mutations.py:
  - drop excluded trios and rows with Phased == 'F'
  - collapse child rows to one locus per (family, chrom, pos, ref, alt)
  - a locus seen in >= 2 offspring of a family is sibling-shared, else unique

Output: results/shared_fraction_sociality_summary.txt
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
SOCIAL        = {"DUM", "MIM", "SAR"}
SUBSOCIAL     = {"AFR", "BIC", "LIN", "TEN"}
KEEP_SPECIES  = SOCIAL | SUBSOCIAL
EXCL_FAMILIES = {"BIC_family2", "BIC_family3", "SAR_family3"}
DATA_FILE     = "data/autosomal_DNM_spread_sheet.tsv"
OUT_FILE      = "results/shared_fraction_sociality_summary.txt"

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


def sociality(species):
    """Return 'social', 'subsocial' or None for a species code."""
    if species in SOCIAL:
        return "social"
    if species in SUBSOCIAL:
        return "subsocial"
    return None


def run_fisher(sub_shared, sub_unique, soc_shared, soc_unique):
    """2x2 Fisher exact test: [[sub_shared, sub_unique], [soc_shared, soc_unique]].
    Returns (OR, p, table)."""
    table = [[sub_shared, sub_unique], [soc_shared, soc_unique]]
    if (sub_shared + sub_unique) == 0 or (soc_shared + soc_unique) == 0:
        return float("nan"), float("nan"), table
    odds_ratio, p_value = fisher_exact(table, alternative="two-sided")
    return float(odds_ratio), float(p_value), table


def pseudocount_or(a, b, c, d, k=1.0):
    """Pseudocount-corrected odds ratio: add k to each cell of [[a, b], [c, d]].
    Gives a finite OR when a cell is 0. Returns (a+k)(d+k)/((b+k)(c+k))."""
    return ((a + k) * (d + k)) / ((b + k) * (c + k))


def fmt_frac(shared, total):
    return f"{(shared / total * 100):.2f}%" if total else "NaN"


def format_fixed_table(rows, headers):
    """Render a list-of-lists as a plain-text fixed-width table."""
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
    df_raw["species"]      = parsed.apply(lambda x: x[0])
    df_raw["family_label"] = parsed.apply(lambda x: x[1])

    n_no_parse = df_raw["species"].isna().sum()
    log(f"Rows with unparseable child   : {n_no_parse}  (excluded)")
    df_raw = df_raw[df_raw["species"].notna()].copy()

    # ------------------------------------------------------------------
    # 3. Keep social + subsocial species
    # ------------------------------------------------------------------
    mask_species = df_raw["species"].isin(KEEP_SPECIES)
    n_excl_species = (~mask_species).sum()
    log(f"Rows excluded (wrong species) : {n_excl_species}  "
        f"(kept: {', '.join(sorted(KEEP_SPECIES))})")
    df = df_raw[mask_species].copy()

    # ------------------------------------------------------------------
    # 4. Exclude flagged families
    # ------------------------------------------------------------------
    mask_family = df["family_label"].isin(EXCL_FAMILIES)
    n_excl_fam  = mask_family.sum()
    excl_list   = sorted(df.loc[mask_family, "family_label"].unique())
    log(f"Rows excluded (excluded fams) : {n_excl_fam}  "
        f"({', '.join(excl_list) if excl_list else 'none present'})")
    df = df[~mask_family].copy()

    # ------------------------------------------------------------------
    # 5. Drop unphaseable (F) rows  (matches DNM_alpha sharing summary)
    # ------------------------------------------------------------------
    mask_f = df["Phased"] == "F"
    log(f"Rows excluded (Phased == F)   : {int(mask_f.sum())}")
    df = df[~mask_f].copy()

    log(f"Rows retained for analysis    : {len(df)}")

    # ------------------------------------------------------------------
    # 6. Classify unique vs sibling-shared (collapse to loci)
    # ------------------------------------------------------------------
    event_key = ["chrom", "pos", "ref", "alt", "family_label"]

    grp = (
        df.groupby(event_key)
          .agg(
              species     = ("species", "first"),
              n_offspring = ("child",   "nunique"),
          )
          .reset_index()
    )
    grp["sociality"]     = grp["species"].apply(sociality)
    grp["mutation_type"] = grp["n_offspring"].apply(
        lambda n: "shared" if n > 1 else "unique"
    )

    # ------------------------------------------------------------------
    # 7. Per-species counts (transparency / cross-check vs DNM_alpha)
    # ------------------------------------------------------------------
    section("PER-SPECIES LOCUS COUNTS (collapsed events, F dropped)")

    sp_rows = []
    for grp_name in ["subsocial", "social"]:
        for sp in sorted(s for s in grp["species"].unique()
                         if sociality(s) == grp_name):
            sub = grp[grp["species"] == sp]
            nu  = int((sub["mutation_type"] == "unique").sum())
            ns  = int((sub["mutation_type"] == "shared").sum())
            sp_rows.append([grp_name, sp, ns, nu, ns + nu, fmt_frac(ns, ns + nu)])

    log(format_fixed_table(
        sp_rows,
        ["Sociality", "Species", "Shared", "Unique", "Total", "Shared_frac"]
    ))

    # ------------------------------------------------------------------
    # 8. Pooled counts per sociality group
    # ------------------------------------------------------------------
    section("POOLED COUNTS BY SOCIALITY")

    pooled = {}
    pool_rows = []
    for grp_name in ["subsocial", "social"]:
        sub = grp[grp["sociality"] == grp_name]
        ns  = int((sub["mutation_type"] == "shared").sum())
        nu  = int((sub["mutation_type"] == "unique").sum())
        pooled[grp_name] = (ns, nu)
        pool_rows.append([grp_name, ns, nu, ns + nu, fmt_frac(ns, ns + nu)])

    log(format_fixed_table(
        pool_rows,
        ["Sociality", "Shared", "Unique", "Total", "Shared_frac"]
    ))

    # ------------------------------------------------------------------
    # 9. Fisher's exact test (pooled, sociality x sharing class)
    # ------------------------------------------------------------------
    section("FISHER'S EXACT TEST (sociality x sharing class)")

    sub_shared, sub_unique = pooled["subsocial"]
    soc_shared, soc_unique = pooled["social"]
    OR, pval, table = run_fisher(sub_shared, sub_unique, soc_shared, soc_unique)

    log("")
    log("2x2 contingency table:")
    log(f"                  Shared        Unique")
    log(f"  Subsocial    {sub_shared:>10d}  {sub_unique:>12d}")
    log(f"  Social       {soc_shared:>10d}  {soc_unique:>12d}")
    log("")

    if np.isnan(pval):
        log("Fisher test skipped: one group has 0 loci.")
    else:
        or_str = "inf" if np.isinf(OR) else f"{OR:.4f}"
        log(f"Odds ratio (raw)            : {or_str}"
            + ("  (zero social shared cell)" if np.isinf(OR) else ""))
        or_pc = pseudocount_or(sub_shared, sub_unique, soc_shared, soc_unique, k=1.0)
        log(f"Odds ratio (pseudocount +1) : {or_pc:.4f}")
        log(f"P-value                     : {pval:.4g}")
        sig = "significant (a=0.05)" if pval < 0.05 else "not significant (a=0.05)"
        log(f"Result                      : {sig}")

    # ------------------------------------------------------------------
    # 10. Notes
    # ------------------------------------------------------------------
    section("NOTES")
    log("1. Sociality groups : social = DUM, MIM, SAR; subsocial = AFR, BIC, LIN, TEN.")
    log("2. Excluded families: BIC_family2, BIC_family3, SAR_family3")
    log("   (SAR_family3 is not present in the data; listed for completeness).")
    log("3. Rows with Phased == 'F' (unphaseable) dropped before locus counting,")
    log("   matching DNM_alpha/workflow_sources/summarize_shared_mutations.py.")
    log("4. Loci collapsed to one event per (chrom, pos, ref, alt, family).")
    log("5. Sibling-shared = locus seen in >= 2 offspring of a family; else unique.")
    log("6. Fisher's exact test pooled across species within each sociality group,")
    log("   two-sided. The raw OR is infinite when the social shared cell is 0;")
    log("   the pseudocount OR adds 1 to each cell for a finite estimate.")
    log("   The pseudocount affects only the OR, not the Fisher p-value.")
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
