import pandas as pd
import numpy as np

df = pd.read_csv("HOGS_plants_datasets.csv")

ALL_PCOLS = ["P0Mte","PcpM0","PcpMte","PchM0","PchMte","PteM0","PteMte","P0M0"]

df["is_consistent"] = df["OMArk_class"].astype(str).str.contains("Consistent", na=False)
df["is_inconsistent"] = df["OMArk_class"].astype(str).str.contains("Inconsistent", na=False)

df["TotalP"] = df[ALL_PCOLS].sum(axis=1)

df["Consistent_total"] = np.where(df["is_consistent"], df["TotalP"], 0)
df["Inconsistent_total"] = np.where(df["is_inconsistent"], df["TotalP"], 0)

df["TE_score"] = df["PteMte"] + df["P0Mte"]
df["Gene_score"] = df["PcpM0"]
df["Pseudo_score"] = df["P0M0"]

hog_species = df.groupby(["HOG_ID", "Species"])[
    ["TE_score","Gene_score","Pseudo_score","Consistent_total","Inconsistent_total"]
].sum().reset_index()

hog_species["AnyInconsistent"] = (hog_species["Inconsistent_total"] > 0).astype(int)

hog_species.to_csv("hog_species_for_GLMM.tsv", sep="\t", index=False)
print("Saved hog_species_for_GLMM.tsv")