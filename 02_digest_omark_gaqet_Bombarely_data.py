from pathlib import Path

import argparse
import json
import subprocess
import yaml


OMARK_CLASS = ["Consistent_Full", "Consistent_Partial",
               "Consistent_Fragment", "Inconsistent_Full",
               "Inconsistent_Partial", "Inconsistent_Fragment",
               "Contamination_Full", "Contamination_Partial",
               "Contamination_Fragment", "Unknown"]

IGNORE = ["Lupinus_albus", "Citrullus_lanatus", "Chelidonium_majus", 
          "Aegilops_comosa", "Mangifera_indica", "Gossypium_arboreum",
          "Humulus_lupulus", "Theobroma_cacao", "Parasponia_andersonii"]

DETENGA_CLASS = ["PcpM0", "PteMte", "P0Mte", "PchMte", "PchM0", "PteM0", "PcpMte", "P0M0"]


def parse_arguments():
    parser = argparse.ArgumentParser(prog="",
                                    description="Summarize AGAT, DeTEnGA and OMA HOGS from GAQET from Bombarely et al",
                                    )
    parser.add_argument("--YAML", "-i", help="input YAML")
    parser.add_argument("--out_filename", "-o", help="out filename")

    return parser.parse_args()


def get_arg_values():
    parser = parse_arguments()
    return {"YAML": Path(parser.YAML),
            "out_filename": Path(parser.out_filename)}


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
        cat_count = detenga[category]
        if cat_count > 0:
            line.append(f"{category}:{detenga[category]}")
    return ";".join(line)


def get_metadata(values):
    accession = values["report"].split("/")[8]
    accession = accession.split("_")
    print(accession)
    try:
        accession = f'{accession[1]}_{accession[2]}.{accession[3]}'
    except:
        accession = "N/A"
    print(accession)
    gaqet_dir = Path(values["gaqet_results"]).parents[0]
    return accession, gaqet_dir


def get_taxid(species):    
        cmd = "datasets summary taxonomy taxon \"{}\"".format(species)
        metadata = subprocess.run(cmd, shell=True, capture_output=True)
        
        if metadata.returncode == 0:
            metadata = json.loads(metadata.stdout)
            tax_metadata = metadata["reports"][0]["taxonomy"]["classification"]
            taxid = tax_metadata["species"]["id"]
            taxon_id = taxid
        else:
            taxon_id = "N/A"
        
        return taxon_id


def main():
    args = get_arg_values()
    input_yaml = args["YAML"]
    out_filename = args["out_filename"]
    with open (out_filename, "w") as summary_out_fhand:
        summary_out_fhand.write(f'HOGID,Accession,Species,TaxID,{",".join(OMARK_CLASS)}\n')
        records = yaml.safe_load(open(input_yaml))
        for species, features in records.items():
            if species in IGNORE:
                continue
            print(species)
            for feature, values in features.items():
                if "NCBI" in feature:
                    accession, gaqet_dir = get_metadata(values)
                    hog_classification = get_omamer_results(gaqet_dir)
                    omark_classification = get_omark_results(gaqet_dir)
                    detenga_classification = get_detenga_results(gaqet_dir)
                    summary_by_class = get_summary_by_class(hog_classification, 
                                                    omark_classification, 
                                                    detenga_classification)
                    taxid = get_taxid(features["species"])
                    for hog, omark_classification in summary_by_class.items():
                        results = []
                        for category, detenga in omark_classification.items():
                            results.append(detenga_line(detenga))
                        line = f'{hog},{accession},{species},{taxid},{",".join(results)}\n'
                        summary_out_fhand.write(line)
                        summary_out_fhand.flush()


if __name__ == "__main__":
    main()