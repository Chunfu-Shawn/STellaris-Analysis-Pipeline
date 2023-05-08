import os
import sys
import scanpy as sc
import time
import utils
import argparse
import Tangram as tg


def main(argsv):
    parser = argparse.ArgumentParser(description="Perform spatial mapping using Tangram method")
    parser.add_argument("--sc_h5ad", help="Path to sc h5ad file")
    parser.add_argument("--st_h5ad", help="Path to st h5ad file selected by the user")
    parser.add_argument("--key_celltype", help="Names of cell type")
    parser.add_argument("--num_epochs", help="Number of epochs for training")
    parser.add_argument("--outDir", help="Path to output directory")
    args = parser.parse_args()

    # Read scRNA-seq and ST data
    utils.info("Read scRNA-seq and ST data...")
    adata_sc = sc.read_h5ad(args.sc_h5ad)
    adata_sp = utils.read_sp(input_h5ad=args.st_h5ad, use_raw=True)

    # 添加GPU环境变量
    os.environ['CUDA_VISIBLE_DEVICES'] = '0,1'
    os.environ['CUDA_LAUNCH_BLOCKING'] = '1'

    # Initiate Tangram instance
    utils.info("Initiate tangram instance...")
    tg_ins = tg.Tangram(ad_sc=adata_sc, ad_sp=adata_sp, output_dir=args.outDir, cell_type_key=args.key_celltype)

    # Perform tangram mapping
    utils.info("Perform tangram mapping...")
    start_time = time.time()
    tg_ins.run_tangram(mode="cells", num_epochs=int(args.num_epochs))
    end_time = time.time()
    elapsed_time = end_time - start_time
    print("Training time: " + str(elapsed_time / 60) + " minutes")

    # Assign single cell positions
    tg_ins.project_cell_positions()
    # save registered sc h5ad
    tg_ins.save()


if __name__ == "__main__":
    main(sys.argv)
