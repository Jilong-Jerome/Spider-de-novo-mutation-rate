"""
Annotate filtered DNM table with SBS96 pyrimidine-normalised trinucleotide category.

Species: BIC, LIN, TEN  (BIC_family2 and BIC_family3 excluded)
Output:  results/parsed_mutations_SBS96.tsv
"""

import re
import os

import pandas as pd

# ---------------------------------------------------------------------------
# Constants  (mirrored from paternal_bias_analysis.py)
# ---------------------------------------------------------------------------
KEEP_SPECIES  = {"BIC", "LIN", "TEN"}
EXCL_FAMILIES = {"BIC_family2", "BIC_family3"}
DATA_FILE     = "data/autosomal_DNM_spread_sheet.tsv"
OUT_FILE      = "results/parsed_mutations_SBS96.tsv"

CHILD_RE = re.compile(r"^([A-Z]+)_(family\d+)_S\d+_offspring$")

COMPLEMENT = str.maketrans("ACGT", "TGCA")

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def parse_child(child_str):
    """Return (species, family_label) or (None, None) if not parseable."""
    m = CHILD_RE.match(str(child_str).strip())
    if m:
        return m.group(1), f"{m.group(1)}_{m.group(2)}"
    return None, None


def sbs96(n_before, ref, alt, n_after):
    """
    Return SBS96 pyrimidine-normalised trinucleotide category string.
    If ref is C or T, use as-is.  If ref is A or G, flip to minus strand.
    Returns 'NA' if any input is missing or contains unexpected characters.
    """
    try:
        ref   = str(ref).strip().upper()
        alt   = str(alt).strip().upper()
        nb    = str(n_before).strip().upper()
        na    = str(n_after).strip().upper()
        valid = set("ACGT")
        if not all(c in valid for c in (ref, alt, nb, na)):
            return "NA"
    except Exception:
        return "NA"

    if ref in ("C", "T"):
        return f"{nb}[{ref}>{alt}]{na}"
    else:  # A or G → flip to pyrimidine strand
        return (
            f"{na.translate(COMPLEMENT)}"
            f"[{ref.translate(COMPLEMENT)}>{alt.translate(COMPLEMENT)}]"
            f"{nb.translate(COMPLEMENT)}"
        )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    # ------------------------------------------------------------------
    # 1. Load data
    # ------------------------------------------------------------------
    df = pd.read_csv(DATA_FILE, sep="\t", dtype=str)
    print(f"Rows loaded              : {len(df)}")

    # ------------------------------------------------------------------
    # 2. Parse child column
    # ------------------------------------------------------------------
    parsed = df["child"].apply(parse_child)
    df["species"]      = parsed.apply(lambda x: x[0])
    df["family_label"] = parsed.apply(lambda x: x[1])

    n_no_parse = df["species"].isna().sum()
    print(f"Rows with unparseable child: {n_no_parse}  (excluded)")
    df = df[df["species"].notna()].copy()

    # ------------------------------------------------------------------
    # 3. Exclude BIC_family2 / BIC_family3 (all species retained)
    # ------------------------------------------------------------------
    mask_family = df["family_label"].isin(EXCL_FAMILIES)
    print(f"Rows excluded (excluded fams): {mask_family.sum()}")
    df = df[~mask_family].copy()

    print(f"Rows retained for annotation : {len(df)}")

    # ------------------------------------------------------------------
    # 5. Classify unique vs shared (individual-level, not collapsed)
    # ------------------------------------------------------------------
    event_key = ["chrom", "pos", "ref", "alt", "family_label"]

    n_offspring = (
        df.groupby(event_key)["child"]
        .transform("nunique")
    )
    df["mutation_type"] = n_offspring.apply(lambda n: "shared" if n > 1 else "unique")

    # ------------------------------------------------------------------
    # 6. Compute SBS96 category
    # ------------------------------------------------------------------
    df["SBS96"] = df.apply(
        lambda row: sbs96(row["n_before"], row["nuc_ref"], row["alt"], row["nuc_after"]),
        axis=1,
    )

    # ------------------------------------------------------------------
    # 7. Write output
    # ------------------------------------------------------------------
    os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)
    df.to_csv(OUT_FILE, sep="\t", index=False)

    # ------------------------------------------------------------------
    # 8. Summary
    # ------------------------------------------------------------------
    sbs_counts = df["SBS96"].value_counts()
    sub_types  = sorted(
        {cat.split("[")[1].split("]")[0] for cat in sbs_counts.index if cat != "NA"}
    )

    print(f"\nRows written             : {len(df)}")
    print(f"Unique SBS96 categories  : {df['SBS96'].nunique()}")
    print(f"Substitution types       : {', '.join(sub_types)}")
    print(f"\nSBS96 category counts (top 10):")
    print(sbs_counts.head(10).to_string())
    print(f"\nOutput written to: {OUT_FILE}")


if __name__ == "__main__":
    main()
