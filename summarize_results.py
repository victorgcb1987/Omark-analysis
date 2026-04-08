import argparse
from glob import glob
from pathlib import Path


def get_arguments():
    parser = argparse.ArgumentParser(prog="",
                                    description="Summarize AGAT, DeTEnGA and OMA HOGS from GAQET results",
                                    )
    parser.add_argument("--input_dir", "-i", help="GAQET input directory")
    parser.add_argumnt("--out_filename", "-o", help="Summary output file")

    return {"input_dir": Path(parser.input_dir),
            "out_filename": Path(parser.out_filename)}


def get_detenga_classification():
    pass

def parse_omamer_results(summary, input):
    omamer_results_folder = input / "OMARK_run"
    for file in list(omamer_results_folder.glob("*.omamer")):
        with open(file) as input_fhand:
            for line in input_fhand:
                if line.startswith("!") or line.startswith("qseqid"):
                    continue
                if line:
                    line.rstrip().split("\t")
                    seqID = line[0]
                    HOGID = line[1]
                    HOGLevel = line[2]
                    summary[seqID] = {"HOGID": HOGID,
                                      "HOGLevel": HOGLevel}

def summary_init(input):
    summary = {}
    seqs_folder = input / "input_sequences"
    for file in list(seqs_folder.glob("*.proteins.fasta")):
        for line in file:
            if line.startwith(">"):
                id = line.split()[0].replace(">", "")
                summary[id] = {"HOG_ID_assignation": "",
                               "HOG_level": ""}
    return summary



def main():
    arguments = get_arguments()
    summary = summary_init(arguments["input_dir"])
    parse_omamer_results(summary, arguments["input_dir"])
    print(summary)








if __name__ == "__main__":
    main()