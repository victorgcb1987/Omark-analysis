import argparse
from glob import glob
from pathlib import Path
from csv import DictReader


OMARK_CLASS = ["Consistent_Full", "Consistent_Partial",
               "Consistent_Fragment", "Inconsistent_Full",
               "Inconsistent_Partial", "Inconsistent_Fragment",
               "Contamination_Full", "Contamination_Partial",
               "Contamination_Fragment", "Unknown"]

DETENGA_CLASS = ["PcpM0", "PteMte", "P0Mte", "PchMte", "PchM0", "PteM0", "PcpMte", "P0M0"]

def parse_arguments():
    parser = argparse.ArgumentParser(prog="",
                                    description="Summarize AGAT, DeTEnGA and OMA HOGS from GAQET results",
                                    )
    parser.add_argument("--input_dir", "-i", help="GAQET input directory")
    parser.add_argument("--metadata_file", "-m", help="Metadata file")
    parser.add_argument("--out_filename", "-o", help="Summary output file")

    return parser.parse_args()


def get_arg_values():
    parser = parse_arguments()
    return {"input_dir": Path(parser.input_dir),
            "out_prefix": Path(parser.out_filename),
            "metadata": Path(parser.metadata_file)}



def get_detenga_results(input_dir):
    detenga_results = {}
    detenga_results_dir = input_dir / "DETENGA_run"
    for file in list(detenga_results_dir.glob("*TE_summary.csv")):
        with open(file) as fhand:
            for line in fhand:
                line = line.rstrip()
                if line.startswith("Transcript_ID"):
                    continue
                if line:
                    line = line.split(";")
                    seqID = line[0]
                    detenga = line[-1]
                    detenga_results[seqID] = detenga
    return detenga_results


def get_omamer_results(input_dir):
    omamer_results = {}
    omamer_results_folder = input_dir / "OMARK_run"
    for file in list(omamer_results_folder.glob("*.omamer")):
        with open(file) as input_fhand:
            for line in input_fhand:
                if line.startswith("!") or line.startswith("qseqid"):
                    continue
                if line:
                    line = line.rstrip().split("\t")
                    seqID = line[0]
                    HOGID = line[1]
                    if "N/A" in HOGID:
                        HOGID = "Unknown"
                    if HOGID not in omamer_results:
                        omamer_results[HOGID] = [seqID]
                    else:
                        omamer_results[HOGID].append(seqID)
    return omamer_results 
                    

def get_omark_results(input_dir):
    omark_results = {}
    omark_consistency_folder = input_dir / "OMARK_run" / "omark" 
    for file in list(omark_consistency_folder.glob("*.ump")):
        with open(file) as fhand:
            for line in fhand:
                if line.startswith(">"):
                    category = line.rstrip().replace(">", "")
                else:
                    seqID = line.rstrip()
                    omark_results[seqID] = category
    return omark_results


def get_summary_by_class(hog_classification, 
                         omark_classification, 
                        detenga_classification):
    summary = {}
    for hog, seqIDs in hog_classification.items():
        summary[hog] = {omark: {detenga: 0 for detenga in DETENGA_CLASS} for omark in OMARK_CLASS}
        for seqID in seqIDs:
            detenga = detenga_classification.get(seqID, "P0M0")
            omark = omark_classification[seqID]
            summary[hog][omark][detenga] += 1
    return summary

def detenga_line(detenga):
    line = []
    for category in DETENGA_CLASS:
        line.append(f"{category}:{detenga[category]}")
        return ";".join(line)

def main():
    '''Needs a dir generated with script 00_run_gaqet_commands'''
    #header1
    #hogID,species,[OMARKS_status]
    arguments = get_arg_values()
    with open(f'{arguments["out_prefix"]}_classification_summary.csv', "w") as summary_out_fhand:
        summary_out_fhand.write(f'HOGID,Species,taxID,{",".join(OMARK_CLASS)}\n')
        for row in DictReader(open(arguments["metadata"]), delimiter=","):
            species = row["species"]
            accession = row["accession"]
            taxid = row["taxid"]
            input_dir = arguments["input_dir"] / "ncbi_dataset" / "data" / accession
            print(input_dir)
            #It assumes that only a subdir by accession!!!!!!
            gaqet_dir = [dir for dir in input_dir.glob("*") if not dir.is_file()][0]
            hog_classification = get_omamer_results(gaqet_dir)
            omark_classification = get_omark_results(gaqet_dir)
            detenga_classification = get_detenga_results(gaqet_dir)
            print(detenga_classification)
            summary_by_class = get_summary_by_class(hog_classification, 
                                                    omark_classification, 
                                                    detenga_classification)
            print(summary_by_class)
            for hog, omark_classification in summary_by_class.items():
                results = []
                for category, detenga in omark_classification.items():
                    results.append(detenga_line(detenga))
                line = f'{hog},{species},{taxid}{",".join(results)}\n'
                summary_out_fhand.write(line)
                





if __name__ == "__main__":
    main()