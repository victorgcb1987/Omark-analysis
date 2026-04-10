import argparse
from glob import glob
from pathlib import Path


def parse_arguments():
    parser = argparse.ArgumentParser(prog="",
                                    description="Summarize AGAT, DeTEnGA and OMA HOGS from GAQET results",
                                    )
    parser.add_argument("--input_dir", "-i", help="GAQET input directory")
    parser.add_argument("--out_filename", "-o", help="Summary output file")

    return parser.parse_args()
        

def get_arg_values():
    parser = parse_arguments()
    return {"input_dir": Path(parser.input_dir),
            "out_filename": Path(parser.out_filename)}


def get_detenga_classification(summary, input):
    detenga_results_dir = input / "DETENGA_run"
    for file in list(detenga_results_dir.glob("*TE_summary.csv")):
        with open(file) as fhand:
            for line in fhand:
                line = line.rstrip()
                if line.startswith("Transcript_ID"):
                    continue
                if line:
                    line = line.split(";")
                    id = line[0]
                    if id not in summary:
                        continue
                        # summary[id] = {"HOG_ID_assignation": "",
                        #                "HOG_level": "",
                        #                "DeTEnGA_status": "",
                        #                "PFAM_IDs" :"",
                        #                "PFAMs_descriptions": "",
                        #                "mrna_TE_classification": "",
                        #                "OMArk_status": ""}
                    summary[id]["DeTEnGA_status"] = line[-1]
                    summary[id]["PFAM_IDs"] = line[3]
                    summary[id]["PFAMs_descriptions"] = line[4]
                    summary[id]["mrna_TE_classification"] = line[2]


def parse_omamer_results(summary, input):
    omamer_results_folder = input / "OMARK_run"
    for file in list(omamer_results_folder.glob("*.omamer")):
        with open(file) as input_fhand:
            for line in input_fhand:
                if line.startswith("!") or line.startswith("qseqid"):
                    continue
                if line:
                    line = line.rstrip().split("\t")
                    print(line)
                    seqID = line[0]
                    HOGID = line[1]
                    HOGLevel = line[2]
                    summary[seqID]["HOG_ID_assignation"] = HOGID
                    summary[seqID]["HOG_level"] = HOGLevel
                    

def summary_init(input):
    summary = {}
    seqs_folder = input / "input_sequences"
    for file in list(seqs_folder.glob("*.proteins_longest_isoform.fasta")):
        with open(file) as fhand:
            for line in fhand:
                if line.startswith(">"):
                    id = line.split()[0].replace(">", "")
                    summary[id] = {"HOG_ID_assignation": "",
                                   "HOG_level": "",
                                   "OMArk_status": "",
                                   "DeTEnGA_status": "",
                                   "PFAM_IDs" :"",
                                   "PFAMs_descriptions": "",
                                   "mrna_TE_classification": ""}
    return summary


def parse_omark_consistency_results(summary, input):
    omark_consistency_folder = input / "OMARK_run" / "omark" 
    for file in list(omark_consistency_folder.glob("*.ump")):
        with open(file) as fhand:
            for line in fhand:
                if line.startswith(">"):
                    category = line.rstrip().replace(">", "")
                else:
                    id = line.rstrip()
                    summary[id]["OMArk_status"] = category




def main():
    arguments = get_arg_values()
    summary = summary_init(arguments["input_dir"])
    parse_omamer_results(summary, arguments["input_dir"])
    get_detenga_classification(summary, arguments["input_dir"])
    parse_omark_consistency_results(summary, arguments["input_dir"])
    with open(arguments["out_filename"], "w") as out_fhand:
        out_fhand.write("SeqID\tHOG_ID\tHOG_level\tOMArk_status\tDeTEnGA_status\tPFAM_IDs\tPFAM_descriptions\tmRNA_classification\n")
        for seq, values in summary.items():
            out_fhand.write(f'{seq}\t{values["HOG_ID_assignation"]}\t{values["HOG_level"]}\t{values["OMArk_status"]}\t{values["DeTEnGA_status"]}\t{values["PFAM_IDs"]}\t{values["PFAMs_descriptions"]}\t{values["mrna_TE_classification"]}\n')








if __name__ == "__main__":
    main()