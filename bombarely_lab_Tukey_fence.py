import csv
import pandas as pd
import sys

OMARK_DUPLICATION_CUTOFF = 50


def tukey_filter(group, k=1.5):
    value_col = "CDS_Models (N)"
    q1 = group[value_col].quantile(0.25)
    q3 = group[value_col].quantile(0.75)
    iqr = q3 - q1
    
    lower = q1 - k * iqr
    upper = q3 + k * iqr
    
    return group[(group[value_col] >= upper)]


def main():
    dataset = {"CDS_Models (N)": [], "Species": [],
               "Family": [], "Annotation_Methodology":[],
               "OMArk 2.0.3 Completeness Results": [],
               "DETENGA_FP%": []}
    fields_required = ["Annotation_Methodology", "Species", "Family", "CDS_Models (N)", "OMArk 2.0.3 Completeness Results", "DETENGA_FP%"]
    input = sys.argv[1]
    with open(input) as input_fhand:
        for row in csv.DictReader(input_fhand, delimiter="\t"):
                valid = True
                for field in fields_required:
                     if not row[field]:
                          valid = False
                if valid:
                    omark_stats = row["OMArk 2.0.3 Completeness Results"]
                    omark_duplication_levels = omark_stats.split(",")[1]
                    omark_duplication_levels = float(omark_duplication_levels.split("%")[0].replace("D:", ""))
                    if omark_duplication_levels >= OMARK_DUPLICATION_CUTOFF:
                        valid = False
                if valid:
                     for field in fields_required:
                          if field == "CDS_Models (N)":
                                dataset[field].append(int(row[field]))
                          else:
                            dataset[field].append(row[field])

    df = pd.DataFrame(dataset)
    df.to_csv("Bombarely_Benchmark.tsv", sep="\t")
    group_col = "Family"

    df_filtered = df.groupby(group_col, group_keys=False).apply(tukey_filter)

    print("Original:", df.shape)
    print("Filtrado:", df_filtered.shape)
                     
    df_filtered.to_csv("Bombarely_Benchmark_CDS_models_top_quantile.tsv", sep="\t")





if __name__ == "__main__":
    main()