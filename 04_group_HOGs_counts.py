import pandas as pd
import re
from collections import defaultdict
from sys import argv

OMARK_CLASS = [
    "Consistent_Full", "Consistent_Partial", "Consistent_Fragment",
    "Inconsistent_Full", "Inconsistent_Partial", "Inconsistent_Fragment",
    "Contamination_Full", "Contamination_Partial", "Contamination_Fragment",
    "Unknown"
]

TE_CLASS = {"PteMte", "P0Mte"}  # unique values from the list

def parse_cell(cell):
    """Parse a cell like 'PcpM0:6;P0M0:13' → dict of {class: count}."""
    if pd.isna(cell) or str(cell).strip() == "":
        return {}
    result = {}
    for item in str(cell).strip().split(";"):
        item = item.strip()
        if not item:
            continue
        if ":" in item:
            key, _, val = item.partition(":")
            try:
                result[key.strip()] = int(val.strip())
            except ValueError:
                result[key.strip()] = 1
        else:
            result[item] = 1
    return result

def classify_species(row):
    """
    Given a row (one species in one HOG), collect all DETENGA classes
    across all OMARK columns and determine if it is:
      - OnlyCoding: all classes are NOT in TE_CLASS
      - TEOnly:     all classes ARE in TE_CLASS
      - Mixed:      mix of both
    Returns one of 'coding', 'te', 'mixed', or None if no data.
    """
    all_classes = set()
    for col in OMARK_CLASS:
        cell_dict = parse_cell(row.get(col, ""))
        all_classes.update(cell_dict.keys())

    if not all_classes:
        return None  # no data for this species in this HOG

    has_te = bool(all_classes & TE_CLASS)
    has_coding = bool(all_classes - TE_CLASS)

    if has_te and has_coding:
        return "mixed"
    elif has_te:
        return "te"
    else:
        return "coding"


def main():
    print("Loading CSV...")
    df = pd.read_csv(argv[1], dtype=str)

    # Fill NaN with empty string for OMARK columns
    for col in OMARK_CLASS:
        if col in df.columns:
            df[col] = df[col].fillna("")

    print(f"Rows: {len(df)}, HOGs: {df['HOGID'].nunique()}")

    results = []

    for hog_id, group in df.groupby("HOGID", sort=False):
        n_species = len(group)

        counts = defaultdict(int)
        for _, row in group.iterrows():
            cat = classify_species(row)
            if cat is not None:
                counts[cat] += 1

        total_classified = counts["coding"] + counts["te"] + counts["mixed"]

        if total_classified == 0:
            only_coding_ratio = 0.0
            te_only_ratio = 0.0
            mixed_ratio = 0.0
        else:
            only_coding_ratio = counts["coding"] / n_species
            te_only_ratio = counts["te"] / n_species
            mixed_ratio = counts["mixed"] / n_species

        results.append({
            "HOGID": hog_id,
            "NumSpecies": n_species,
            "OnlyCodingRatio": round(only_coding_ratio, 6),
            "TEOnlyRatio": round(te_only_ratio, 6),
            "MixedRatio": round(mixed_ratio, 6),
        })

    out_df = pd.DataFrame(results)
    #Removing HOGs with only PchMte
    out_df = out_df.loc[out_df["MinCounts"] < 1]
    print(out_df.loc[out_df["MinCounts"] < 1])
    out_path = "HOG_ratios.csv"
    out_df.to_csv(out_path, index=False)
    print(f"Saved {len(out_df)} rows to {out_path}")
    print(out_df.head(10).to_string())

if __name__ == "__main__":
    main()
