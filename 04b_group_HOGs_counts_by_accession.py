import pandas as pd
import sys

OMARK_CLASS = [
    "Consistent_Full", "Consistent_Partial", "Consistent_Fragment",
    "Inconsistent_Full", "Inconsistent_Partial", "Inconsistent_Fragment",
    "Contamination_Full", "Contamination_Partial", "Contamination_Fragment",
    "Unknown"
]

TE_CLASS = {"PteMte", "P0Mte"}
IGNORE = {"PchMte", "PchM0", "P0M0"}


def parse_cell(cell):
    """Parse 'PcpM0:6;P0M0:13' -> {'PcpM0': 6, 'P0M0': 13}"""
    if pd.isna(cell) or str(cell).strip() == "":
        return {}
    result = {}
    for item in str(cell).strip().split(";"):
        item = item.strip()
        if not item:
            continue
        if ":" in item:
            key, _, val = item.partition(":")
            if key in IGNORE:
                continue
            try:
                result[key.strip()] = int(val.strip())
            except ValueError:
                result[key.strip()] = 1
        else:
            result[item] = 1
    return result


def total_counts(row):
    """Sum all numeric values across all OMARK columns for a single row."""
    total = 0
    for col in OMARK_CLASS:
        total += sum(parse_cell(row.get(col, "")).values())
    return total


def classify_row(row):
    """
    Classify a single row (one HOG x Accession) based on OMARK columns.
    Returns 'coding', 'te', 'mixed', or None (no data).
    Since each HOG+Accession combination has exactly one row, no merging needed.
    """
    all_classes = set()
    for col in OMARK_CLASS:
        all_classes.update(parse_cell(row.get(col, "")).keys())
    if not all_classes:
        return None
    has_te = bool(all_classes & TE_CLASS)
    has_coding = bool(all_classes - TE_CLASS)
    if has_te and has_coding:
        return "mixed"
    if has_te:
        return "te"
    return "coding"


def main():
    print("Loading CSV...")
    df = pd.read_csv(
        sys.argv[1],
        dtype=str
    )
    for col in OMARK_CLASS:
        if col in df.columns:
            df[col] = df[col].fillna("")

    print("Rows:", len(df), "| HOGs:", df["HOGID"].nunique())

    # Classify every row and compute total counts up front
    df["_cat"] = df.apply(classify_row, axis=1)
    df["_counts"] = df.apply(total_counts, axis=1)
    print(df)
    print(len(df))
    df = df.loc[df["_counts"] > 0]
    print(len(df))
    results = []
    for hog_id, hog_group in df.groupby("HOGID", sort=False):
        n_acc = len(hog_group)

        # One row per accession guaranteed — count directly
        counts = hog_group["_cat"].value_counts()

        results.append({
            "HOGID": hog_id,
            "NumAccessions": n_acc,
            "OnlyCodingRatio": round(counts.get("coding", 0) / n_acc, 6),
            "TEOnlyRatio":     round(counts.get("te",     0) / n_acc, 6),
            "MixedRatio":      round(counts.get("mixed",  0) / n_acc, 6),
            "MinCounts":       int(hog_group["_counts"].min()),
            "MaxCounts":       int(hog_group["_counts"].max()),
        })

    out = pd.DataFrame(results)
    #Removing HOGs with only PchMte
    #out = out.loc[out["MaxCounts"] >= 1]
    #print(out.loc[out["MaxCounts"] >= 1])
    out_path = "HOG_ratios.csv"
    out.to_csv(out_path, index=False)
    print("Saved", len(out), "rows ->", out_path)
    print(out.head(10).to_string())


main()