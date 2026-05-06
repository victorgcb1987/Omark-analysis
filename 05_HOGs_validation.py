import argparse
from omadb import Client
from pathlib import Path
import csv
import statistics




def parse_arguments():
    parser = argparse.ArgumentParser(prog="",
                                    description="Summarize AGAT, DeTEnGA and OMA HOGS from GAQET",
                                    )
    parser.add_argument("--input_file", "-i", help="input dir")
    parser.add_argument("--out_prefix", "-o", help="out filename")
    return parser.parse_args()


def get_arg_values():
    parser = parse_arguments()
    return {"input_file": Path(parser.input_file),
            "prefix": Path(parser.out_prefix)}


def get_hog_values(hog):
    connection = Client()
    hog_data = connection.hogs[hog]
    return hog_data['completeness_score'], hog_data["description"]


def main():
    fail_log = []
    analized_hogs = {}
    args = get_arg_values()
    log_fhand = open(f'{args["prefix"]}_errors.txt', "w")
    with open(args["input_file"]) as input_fhand:
        for row in csv.DictReader(input_fhand, delimiter=","):
            if row["HOGID"] == "Unknown":
                continue
            if float(row["TEOnlyRatio"]) >= 1:
                hog = row["HOGID"]
                if "." in hog:
                    hog = hog.split(".")[0]
                proteins_in_hog = {}
                if hog not in analized_hogs:
                    try:
                        print(hog)
                        connection = Client()
                        hogs = connection.hogs
                        proteins = hogs.members(hog)
                        try:
                            for protein in proteins:
                                species = protein["species"]["species"]
                                protein_id = protein["omaid"]   
                                protein_data = connection.entries[protein_id]
                                domains = protein_data.domains["regions"]
                                domains = [domain["name"] for domain in domains]
                                isoforms = protein_data.isoforms
                                for isoform in isoforms:
                                    if isoform["is_main_isoform"]:
                                        number_of_exons = isoform["nr_exons"]
                                        feats = {"protein_id": protein_id,
                                                 "domains": set(domains),
                                                 "num_exons": int(number_of_exons)}
                                        if species not in proteins_in_hog:
                                            proteins_in_hog[species] = [feats]
                                        else:
                                            proteins_in_hog[species].append(feats)  
                            analized_hogs[hog] = proteins_in_hog     
                        except:
                                msg = f'Protein failed for species {species}: {protein_id} {hog}\n'
                                log_fhand.write(msg)
                                log_fhand.flush()
                    except:
                        msg = f'HOG failed: {hog}\n'
                        log_fhand.write(msg)
                        log_fhand.flush()
                        continue
               

    with open(f'{args["prefix"]}_results.tsv', "w") as out_fhand:
        header = "HOGID,Description,CompletnessScore,NumberOfSpecies,Min,MedianBySpecies,Max,Total,MinExons,MedianExons,MaxExons,DomainsFound\n"
        out_fhand.write(header)
        for hog, values in analized_hogs.items():
            completness, description = get_hog_values(hog)
            num_species = len(values.keys())
            domains_found = {}
            num_proteins = []
            num_exons = []
            for species, proteins in values.items():
                num_proteins.append(len(proteins))
                for protein in proteins:
                    for domain in protein["domains"]:
                        num_exons.append(protein["num_exons"])
                        if domain not in domains_found:
                            domains_found[domain] = 1
                        else:
                            domains_found[domain] += 1
            min_prots = min(num_proteins)
            max_prots = max(num_proteins)
            median_prots = statistics.median(num_proteins)
            min_exons = min(num_exons)
            max_exons = max(num_exons)
            total = sum(num_proteins)
            median_exons = statistics.median(num_exons)
            line = f'{hog},{description},{completness},{num_species},'
            line += f'{min_prots},{median_prots},{max_prots},{total},'
            line += f'{min_exons},{median_exons},{max_exons}'
            line += f'{";".join([f"{name}:{value}" for name, value in domains_found.items()])}'
            out_fhand.write(line+"\n")
            

    
if __name__ == "__main__":
    main()