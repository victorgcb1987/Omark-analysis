from pathlib import Path
import argparse
import pandas as pd
import matplotlib.pyplot as plt




def parse_arguments():
    parser = argparse.ArgumentParser(prog="",
                                    description="Plot DeTEnGA and OMArk results",
                                    )
    parser.add_argument("--input_dir", "-i", 
                        help="Input summary file")
    parser.add_argument("--keep_hierarchy", 
                        "-k", help="Keep HOG subdivisons hierarchy. Default: False",
                        action="store_true")
    parser.add_argument("--out_dir", "-o", help="Directory for plots")
    parser.add_argument("--species", "-s", help="Species Name")

    return parser.parse_args()
        

def get_arg_values():
    parser = parse_arguments()
    return {"input_dir": Path(parser.input_dir),
            "out_dir": Path(parser.out_dir),
            "keep_hierarchy": parser.keep_hierarchy,
            "species": parser.species}


def main():
    arguments = get_arg_values()
    input_fpath = arguments["input_dir"]
    out_dir = arguments["out_dir"]
    species = " ".join(arguments["species"].split("_"))
    if not out_dir.exists():
        out_dir.mkdir(out_dir)
    df = pd.read_csv(input_fpath, sep="\t")
    df["OMArk_status"] = df["OMArk_status"].fillna("Unknown")
    df["DeTEnGA_status"] = df["DeTEnGA_status"].fillna("P0M0")
    if not arguments["keep_hierarchy"]:
        print("Grouping all HOG subidivions into the ancestral one")
    df["HOG_ID"] = df["HOG_ID"].str.split(".").str[0]

    # Tabla cruzada: OMArk en eje X, DeTEnGA como capas apiladas
    ct = pd.crosstab(df["OMArk_status"], df["DeTEnGA_status"])

    # Normalizar a porcentaje
    ct_pct = ct.div(ct.sum(axis=1), axis=0) * 100

    # Etiquetas con el número de genes por categoría
    counts = ct.sum(axis=1)
    ct_pct.index = [f"{label}\n(n={counts[label]:,})" for label in ct_pct.index]

    # Plot
    ax = ct_pct.plot(
    kind="bar",
    stacked=True,
    figsize=(10, 6),
    colormap="tab10",
    edgecolor="white"
    )
    ax.set_xlabel("OMArk status")
    ax.set_ylabel("Genes per OMArk category(%)")
    ax.set_title(f"{species}, DeTEnGA status distribution by OMArk status")
    ax.legend(title="DeTEnGA status", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(out_dir / f"{arguments["species"]}_DeTEnGA_distribution_by_OMArk_status.svg", dpi=300)
    plt.show()

    
    




if __name__ == "__main__":
    main()