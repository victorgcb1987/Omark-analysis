import gzip
import json
import sys
import time

from statistics import median
from pathlib import Path
from omadb import Client


def get_names_from_CATH_name_file(names_fhand):
    cath_names = {}
    for line in names_fhand:
        line = line.rstrip()
        if line.startswith("#"):
            continue
        elif line:
            cath_id = line.split()[0]
            cath_name = line.split(":")[-1]
            cath_names[cath_id] = cath_name
    return cath_names


def get_hog_descriptions_from_file(hog_descriptions_fhand):
    return {line.split(",")[0]: line.split(",")[-1].rstrip() for line in hog_descriptions_fhand if line}


def get_completness_scores_from_file(completness_fhand):
    return {line.split(",")[0]: line.split(",")[-1].rstrip() for line in completness_fhand if line}


def get_results(entries, cath_descriptions, HOG_descriptions, 
                completness_scores, merge_subhogs=True):
    results = {}
    for entry in entries:
        hog = entry["hogid"]
        if "." in hog and merge_subhogs:
            hog = hog.split(".")[0]
        if hog not in results:
            results[hog] = {"completness_score": completness_scores[hog],
                            "num_exons": [], "domains": {"n/a": 0},
                            "description": HOG_descriptions[hog].replace(",", " ")
                            }
        domains = set([domain[0] for domain in entry["domains"]])
        if not domains:
            results[hog]["domains"]["n/a"] += 1
        else:
            for domain in domains:
                domain_name = cath_descriptions.get(domain, "-").replace(",", " ")
                if domain_name not in results[hog]["domains"]:
                    results[hog]["domains"][domain_name] = 1
                else:
                    results[hog]["domains"][domain_name] += 1
        results[hog]["num_exons"].append(entry["nr_exons"])
    return results
    




def main():
    main_dir = Path(sys.path[0])
    docs_dir = main_dir / "docs"
    cath_names_input = docs_dir / "cath-names_12_06_2026.txt.gz"
    HOG_members_input = docs_dir / "TE_OMA_database_12_05_2026.json.gz"    
    HOG_descriptions_input = docs_dir / "HOGs_descriptions.csv"
    HOGs_TEs_completness_scores = docs_dir / "hog_TEs_completness_score_12_05_2026.csv"
    
    with gzip.open(HOG_members_input, "rt") as hog_members_fhand:
        entries = json.loads(hog_members_fhand.read())
    
    with gzip.open(cath_names_input, "rt") as names_fhand:
        cath_descriptions = get_names_from_CATH_name_file(names_fhand)
        print(cath_descriptions)
    
    with open(HOG_descriptions_input) as hog_descriptions_fhand:
        HOG_descriptions = get_hog_descriptions_from_file(hog_descriptions_fhand)
    
    with open(HOGs_TEs_completness_scores) as completness_fhand:
        completness_scores = get_completness_scores_from_file(completness_fhand)
    
    results = get_results(entries, cath_descriptions, HOG_descriptions, completness_scores)
    with open(docs_dir / "TE_validation.csv", "w") as out_fhand:
        lines = ["HOGID,HOG_Description,CompletnessScore,NumMembers,MaxExons,MinExons,MedianExons,DomainsFound"]
        for hog, values in results.items():
            exons = values["num_exons"]
            max_exons = max(exons)
            min_exons = min(exons)
            median_exons = median(exons)
            domains = ";".join([f"{name}:{occurrences}" for name, occurrences in values["domains"].items()])
            line = f'{hog},{values["description"]},{values["completness_score"]},'
            line += f'{len(exons)},{max_exons}({exons.count(max_exons)}),{min_exons}({exons.count(min_exons)}),{median_exons},{domains},'
            lines.append(line)
        out_fhand.write("\n".join(lines))
    for entry in entries:
        if entry["hogid"] == "HOG:E0801531":
            if entry["nr_exons"] == 0:
                print(entry["omaid"], entry["nr_exons"])



if __name__ == "__main__":
    main()