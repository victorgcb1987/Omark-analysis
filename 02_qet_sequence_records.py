from pathlib import Path

import argparse
from csv import DictReader


def parse_arguments():
    parser = argparse.ArgumentParser(prog="",
                                    description="Summarize AGAT, DeTEnGA and OMA HOGS from GAQET",
                                    )
    parser.add_argument("--input_dir", "-i", help="input dir")
    parser.add_argument("--metadata", "-m", help="input metadata")
    parser.add_argument("--out_filename", "-o", help="out filename")
    parser.add_argument("--descriptions", "-d", help="HOGs description file")
    return parser.parse_args()


def get_arg_values():
    parser = parse_arguments()
    return {"input_dir": Path(parser.input_dir),
            "metadata": Path(parser.metadata),
            "out_filename": Path(parser.out_filename),
            "descriptions": Path(parser.descriptions)}


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
                    pfams = line[3]
                    description = line[4]
                    tesorter_class = line[2]
                    detenga_results[seqID] = {"status": detenga,
                                              "pfams": pfams,
                                              "description": description,
                                              "tesorter": tesorter_class}
    return detenga_results

def get_description(HOG, descriptions):
    check = HOG
    if check == "Unknown":
        return "Unknown"
    if "." in check:
        check = check.split(".")[0]
    return descriptions[check]

def get_omamer_results(input_dir,descriptions):
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
                    omamer_results[seqID] = {"HOG": HOGID,
                                             "description": get_description(HOGID, descriptions)}
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


def main():
    args = get_arg_values()
    out_filename = args["out_filename"]
    descriptions = {line.split(",")[0]: line.split(",")[-1].rstrip() for line in open(args["descriptions"])}
    with open (out_filename, "w") as summary_out_fhand:
        summary_out_fhand.write("SeqID,Accession,Species,TaxID,HOGID,HOG_Description,OMARK_status,DeTEnGA_status,PFAMs,PFAMs_Description,TesorterClass\n")
        for row in DictReader(open(args["metadata"]), delimiter=","):
            species = row["species"]
            accession = row["accession"]
            taxid = row["taxid"]
            input_dir = args["input_dir"] / "ncbi_dataset" / "data" / accession
            gaqet_dir = [dir for dir in input_dir.glob("*") if not dir.is_file()][0]
            hog_classification = get_omamer_results(gaqet_dir, descriptions)
            omark_classification = get_omark_results(gaqet_dir)
            detenga_classification = get_detenga_results(gaqet_dir)
            
            for seqID, hog in hog_classification.items():
                seq_OMArk_class = omark_classification[seqID]
                if seqID in detenga_classification:
                    detenga_status = detenga_classification[seqID]["status"]
                    pfams = detenga_classification[seqID]["pfams"]
                    description = detenga_classification[seqID]["description"]
                    tesorter = detenga_classification[seqID]["tesorter"]
                else:
                    detenga_status = "P0M0"
                    pfams = ""
                    description = ""
                    tesorter = ""

                line = f'{seqID},{accession},{species},{taxid},{hog["HOG"]},'
                line += f'{hog["description"]},{seq_OMArk_class},'
                line += f'{detenga_status},{pfams},{description},{tesorter}\n'
                summary_out_fhand.write(line)
                summary_out_fhand.flush()


if __name__ == "__main__":
    main()