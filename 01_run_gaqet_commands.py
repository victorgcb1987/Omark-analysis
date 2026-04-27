
from sys import argv as args
from pathlib import Path
from csv import DictReader

import subprocess


def main():
    '''Èxample of metadata
    species,accession,taxid
    Mus musculus,GCA_000002165.1,10090
    Rattus norvegicus,GCA_000002265.1,10116
    Cricetulus griseus,GCA_000223135.1,10029
    Neotoma lepida,GCA_001675575.1,56216'''
    main_directory = Path(args[1]) / "ncbi_dataset" / "data"
    if not main_directory.exists():
        raise RuntimeError(f"{main_directory} does not exists")
    metadata = args[2]
    log = args[3]
    gaqet_config = args[4]
    with open(log, "w") as log_fhand:
        with open(metadata) as metadata_fhand:
            for row in DictReader(metadata_fhand, delimiter=","):
                filepath = main_directory / row["accession"]
                if not filepath.exists():
                    print(row)
                    raise RuntimeError(f"{filepath} does not exists")
                else:
                    assembly_file = list(filepath.glob("*fna"))[0]
                    annotation_file = list(filepath.glob("*.gff"))[0]
                    spp = f'{"_".join(row["species"].split())}_{row["accession"]}'
                    gaqet_out = filepath / spp
                    cmd = f'GAQET -i {gaqet_config} -o {gaqet_out} -s {spp} '
                    cmd += f'-g {assembly_file} -a {annotation_file} '
                    cmd += f'-t {row["taxid"]}'
                    print(cmd)
                    results = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.DEVNULL)
                    if results.returncode == 0:
                        msg = "DONE,"
                    else:
                        msg = "FAILED,"
                    msg += f'{row["accession"]},{row["species"]},{cmd}\n'
                    log_fhand.write(msg)
                    
                
if __name__ == "__main__":
    main()