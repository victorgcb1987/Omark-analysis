from csv import DictReader
from sys import argv


TEs_classes = ["PteMte", "P0Mte"]
OMARK_CLASSES = ["Consistent_Full", "Consistent_Partial", 
               "Consistent_Fragment", "Inconsistent_Full", 
               "Inconsistent_Partial", "Inconsistent_Fragment", 
               "Contamination_Full", "Contamination_Partial", 
               "Contamination_Fragment", "Unknown"]

def TEs_found_in_HOG(row):
    for omark_class in OMARK_CLASSES:
        for te_class in TEs_classes:
            if te_class in row[omark_class]:
                return True
    return False


def main():
    input = argv[1]
    with open(input) as input_fhand:
        current_accession = ""
        hogs_found = []
        hogs_with_tes_found = []
        num_accessions = [0]
        num_HOGS = [0]
        num_HOGs_with_TEs = [0]
        new_hogs = 0
        new_hogs_with_Tes = 0
        first = True
        for row in DictReader(input_fhand, delimiter=","):
            if "." in row["HOGID"]:
                row["HOGID"] = row["HOGID"].split(".")[0]
            if first:
                current_accession = row["Accession"]
                first = False
            if row["Accession"] != current_accession:
                num_accessions.append(num_accessions[-1] + 1)
                num_HOGS.append(num_HOGS[-1]+new_hogs)
                num_HOGs_with_TEs.append(num_HOGs_with_TEs[-1] + new_hogs_with_Tes)
                new_hogs = 0
                new_hogs_with_Tes = 0
                current_accession = row["Accession"]
            if row["HOGID"] not in hogs_found:
                hogs_found.append(row["HOGID"])
                new_hogs += 1
            if TEs_found_in_HOG(row):
                if row["HOGID"] not in hogs_with_tes_found:
                    hogs_with_tes_found.append(row["HOGID"])
                    new_hogs_with_Tes += 1
        print("Num_Accessions", num_accessions)
        print("Num_HOGS", num_HOGS)
        print("Num_HOGS_with_TEs", num_HOGs_with_TEs)
        with open("saturation_analysis.tsv", "w") as out_fhand:
            out_fhand.write("NumAccessions,NumHOGS,NumTEsHOGs\n")
            for pos in range(0, len(num_accessions)):
                out_fhand.write(f'{num_accessions[pos]},{num_HOGS[pos]},{num_HOGs_with_TEs[pos]}\n')



if __name__ == "__main__":
    main()