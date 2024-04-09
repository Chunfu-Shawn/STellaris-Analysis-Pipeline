import argparse
import os
import pandas as pd
import scanpy as sc
import numpy as np
import scipy.sparse as sp

st_dir="/home/user/data/qijt/stdataset/h5ad_files_only_count"

def main(argsv):
    parser = argparse.ArgumentParser(description="Prepare ST input for CytoSPACE.")
    parser.add_argument("--dataset", help="ST dataset")
    parser.add_argument("--section", help="ST section")
    parser.add_argument("--outDir",type=str,help="Path to output directory",default="./")

    args = parser.parse_args()

    st_h5ad_file = ".".join([args.section, "input.h5ad"])
    st_h5ad_file = os.path.join(st_dir, args.dataset, args.section, st_h5ad_file)

    if not os.path.exists(st_h5ad_file):
        raise ValueError("The h5ad file for the requested ST section doesn't exist!")

    # Read data
    adata = sc.read_h5ad(st_h5ad_file)

    # Save ST count
    if sp.isspmatrix(adata.X):
        count_cytospace = pd.DataFrame(adata.X.toarray().T, index=list(adata.var_names), columns=list(adata.obs_names))
    else:
        count_cytospace = pd.DataFrame(adata.X.T, index=list(adata.var_names), columns=list(adata.obs_names))

    count_cytospace.index.name = "gene"
    count_cytospace.to_csv(os.path.join(args.outDir, 'st_counts.txt'), sep="\t")

    # Save ST coordinate
    coord_cytospace = pd.DataFrame(adata.obsm['spatial'], index=list(adata.obs_names), columns=['row', 'col'])
    coord_cytospace.index.name = "SpotID"
    coord_cytospace.to_csv(os.path.join(args.outDir, 'st_coordinates.txt'), sep="\t")


if __name__ == "__main__":
    import sys

    main(sys.argv)