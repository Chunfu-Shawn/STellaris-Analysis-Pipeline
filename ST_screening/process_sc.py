import os
import re
import sys
import argparse
import logging
import json
import pandas as pd
import numpy as np
import scanpy as sc
import warnings
warnings.filterwarnings("ignore")

logger = logging.getLogger("process_sc")
logger.setLevel(logging.INFO)
handler=logging.StreamHandler(sys.stdout)
handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s:%(step)s - %(levelname)s - %(message)s"))
logger.addHandler(handler)


def match_ens_id(string):
    pattern = r'^ENS.*000.*$'
    return bool(re.match(pattern, string))


def identify_ens_id(gene_names):
    for g in gene_names:
        if match_ens_id(g):
            return True
    return False


def read_count(count_path):

    if count_path.endswith('.txt.gz') or count_path.endswith('.tsv.gz') or count_path.endswith('.txt') or count_path.endswith('.tsv'):

        try:
            return sc.read_text(count_path)
        except:
            raise ValueError("Can't read the raw count matrix in .txt.gz, .tsv.gz or .txt format! \
                              Check if the length of the header (genes) matches the number of columns in the matrix. \
                              Please upload raw count matrix that meets the requirements.")

    elif count_path.endswith('.csv.gz') or count_path.endswith('.csv'):

        try:
            return sc.read_csv(count_path)
        except:
            raise ValueError("Can't read the raw count matrix in .csv.gz or .csv format! \
                              Check if the length of the header (genes) matches the number of columns in the matrix. \
                              Please upload raw count matrix that meets the requirements.")

    elif count_path.endswith('.txt.zip') or count_path.endswith('.tsv.zip'):
        import zipfile
        import anndata
        data = []
        row_names = []

        try:
            with zipfile.ZipFile(count_path, 'r') as z:
                with z.open(z.namelist()[0], 'r') as f:
                    # gene names (var)
                    col_names = f.readline().decode().split('\t')[1:]
                    for line in f:
                        line_lst = line.rstrip().decode().split('\t')
                        row_names.append(line_lst[0])
                        data.append(line_lst[1:])
        except:
            raise ValueError("Can't read the raw count matrix in .txt.zip or .tsv.zip format! \
                              Please upload raw count matrix that meets the requirements.")

        try:
            data = np.array(data, dtype=int)
        except ValueError:
            try:
                data = np.array(data, dtype=float)
            except:
                raise ValueError("Can't read the raw count matrix in .txt.zip or .tsv.zip format! Please upload raw count matrix that meets the requirements.")

        if len(col_names) != np.shape(data)[1]:
            raise ValueError("The length of the header (genes) doesn't match the number of columns in the matrix!")

        try:
            return anndata.AnnData(data, obs=dict(obs_names=row_names), var=dict(var_names=col_names))
        except:
            raise ValueError('Unsupported raw count matrix file format! Please upload tab-delimited (txt/tsv) or comma-delimited (csv) file compressed with gzip (.gz) or zip (.zip)!')


    elif count_path.endswith('.csv.zip'):
        import zipfile
        import anndata
        data = []
        row_names = []

        try:
            with zipfile.ZipFile(count_path, 'r') as z:
                with z.open(z.namelist()[0], 'r') as f:
                    # gene names (var)
                    col_names = f.readline().decode().split(',')[1:]
                    for line in f:
                        line_lst = line.rstrip().decode().split(',')
                        row_names.append(line_lst[0])
                        data.append(line_lst[1:])
        except:
            raise ValueError("Can't read the raw count matrix in .csv.zip format! Please upload raw count matrix that meets the requirements.")

        try:
            data = np.array(data, dtype=int)
        except ValueError:
            try:
                data = np.array(data, dtype=float)
            except:
                raise ValueError("Can't read the raw count matrix in .csv.zip format! Please upload raw count matrix that meets the requirements.")

        if len(col_names) != np.shape(data)[1]:
            raise ValueError("The length of the header (genes) doesn't match the number of columns in the matrix!")

        try:
            return anndata.AnnData(data, obs=dict(obs_names=row_names), var=dict(var_names=col_names))
        except:
            raise ValueError('Unsupported raw count matrix file format! Please upload tab-delimited (txt/tsv) or comma-delimited (csv) file compressed with gzip (.gz) or zip (.zip)!')

    else:
        raise ValueError('Unsupported raw count matrix file format! Please upload tab-delimited (txt/tsv) or comma-delimited (csv) file compressed with gzip (.gz) or zip (.zip)!')


def read_label(label_path):
    
    if label_path.endswith('.txt.gz') or label_path.endswith('.tsv.gz') or label_path.endswith('.txt.zip') or label_path.endswith('.tsv.zip') or label_path.endswith('.txt') or label_path.endswith('.tsv'):
        try:
            return pd.read_csv(label_path,index_col=0,header=0,sep='\t')
        except:
            raise ValueError("Can't read the label file! Please upload label file that meets the requirements.")

    elif label_path.endswith('.csv.gz') or label_path.endswith('.csv.zip') or label_path.endswith('.csv'):
        try:
            return pd.read_csv(label_path,index_col=0,header=0)
        except:
            raise ValueError("Can't read the label file! Please upload label file that meets the requirements.")

    else:
        raise ValueError('Unsupported label file format! Please upload tab-delimited (txt/tsv) or comma-delimited (csv) file compressed with gzip (.gz) or zip (.zip)!')


def read_sc(count,label,key_celltype='cell_type'):

    adata = read_count(count)

    if identify_ens_id(adata.var_names):
        raise ValueError('The gene names you provided include Ensembl ID. Please provide a matrix that only contains gene symbols!')

    label = read_label(label)

    if not all(adata.obs.index == label.index):
        raise ValueError('The cell ids of count matrix and label file are not consistent!')

    if not key_celltype in label.columns:
        raise  ValueError(key_celltype + ' is not found in column names of label file, pleas specify the column name of cell type.')

    adata.obs = label

    #special_characters = "!@#$%^&*()[]{};:,./<>?\|`~-=_+"
    if pd.api.types.is_numeric_dtype(adata.obs[key_celltype]):
        adata.obs[key_celltype] = "celltype_" + adata.obs[key_celltype].astype(str)
    num_cell_types = len(adata.obs[key_celltype].unique())
    if num_cell_types > 1000:
        raise ValueError("The number of cell type in your scRNA-seq data is {} (>1000). It seems that you have mistakenly used continuous values as cell types!".format(num_cell_types))
    else:
        special_characters = "!@#$%^&*;:,/<>?\|`~="
        adata.obs[key_celltype] = [ z.translate ({ord(c): "_" for c in special_characters}) for z in adata.obs[key_celltype] ]

    return adata


def preprocess(adata,key_celltype='cell_type',min_genes=200,min_mt_pct=20):
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    num_raw = adata.obs[key_celltype].value_counts().sort_values(ascending=False)

    sc.pp.filter_cells(adata, min_genes=min_genes)
    adata = adata[adata.obs.pct_counts_mt < min_mt_pct,:].copy()

    if len(adata.obs[key_celltype].unique()) <= 1: # Number of cell types is 1, raise error.
        raise ValueError("The number of cell types is 1 or 0 such that no cell-cell interactions can be assessed! \
                          Please provide scRNA-seq data with more cell types.")
    elif adata.n_obs < 100:
        raise ValueError("The number of cells after filtration is less than 100, please provide more cells!")

    sc.pp.normalize_total(adata, target_sum=10000, inplace=True)
    sc.pp.log1p(adata)

    num_kept = adata.obs[key_celltype].value_counts().reindex(num_raw.index)
    num_filter = num_raw.subtract(num_kept,fill_value=0).reindex(num_raw.index)

    filter_summary = {'Cell_type': list(num_filter.index), 'Low_quality': list(num_filter.values), 'Passed':list(num_kept.values)}

    meta_info = {'Cell_number':len(adata.obs_names), 'Gene_number': len(adata.var_names)}

    return meta_info, filter_summary, adata


def calculate_markers(adata,key_celltype):
    key_added = 'rank_genes_groups' + '_' + key_celltype
    value_counts = adata.obs[key_celltype].value_counts()
    filtered_value_counts = value_counts[value_counts >= 3]
    if len(filtered_value_counts) >= 2:
        print("Computing markers for {}".format(key_celltype))
        sc.tl.rank_genes_groups(
            adata,
            groupby=key_celltype,
            key_added=key_added,
            method="t-test",
            groups=filtered_value_counts.index.to_list(),
            use_raw=False
        )
    else:
        print("Skip computing marker genes due to few cells!")

    return adata


def extract_markers(adata,de_key='rank_genes_groups_cell_type',logfoldchange=1,pval_adj=0.01):
    marker = {}
    marker['all_genes'] = list(adata.var.index)
    groups = adata.uns[de_key]['names'].dtype.names
    for group in groups:
        de_df = sc.get.rank_genes_groups_df(adata, key=de_key,group=group)
        marker[group] = list(de_df[(de_df['logfoldchanges'] >= logfoldchange) & (de_df['pvals_adj'] >= pval_adj)]['names'])
    return marker


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)


def main(argsv):
    parser = argparse.ArgumentParser(description="Calculate marker genes for scRNA-seq grouped by annotated cell type.")
    parser.add_argument("--count", help="Count matrix of scRNA-seq")
    parser.add_argument("--label", help="Metadata including cell type of scRNA-seq")
    parser.add_argument("--key_celltype", help="Column name of cell type",default="cell_type")
    parser.add_argument("--min_genes",type=int, help="Minimum number of genes expressed required for a cell to pass filtering",default=200)
    parser.add_argument("--max_mt_pct", type=int, help="Maximum percent of mitochondrial UMIs for a cell to pass filtering",default=20)
    parser.add_argument("--outDir",type=str,help="Path to output directory",default="./")

    args = parser.parse_args()

    if (
            os.path.exists(os.path.join(args.outDir,'out','json','meta_info.preprocessing.json')) and
            os.path.exists(os.path.join(args.outDir,'out','json','filter_summary.preprocessing.json')) and
            os.path.exists(os.path.join(args.outDir,'sc.h5ad')) and
            os.path.exists(os.path.join(args.outDir,'out','json','sc_markers.json'))
    ):
        logger.info("All Finished!", extra={'step': 'Main'})
        sys.exit(0)

    #### Read data
    logger.info("Reading data...", extra={'step': 'Main'})
    adata = read_sc(args.count,args.label,key_celltype=args.key_celltype)

    #### Preprocessing
    logger.info("Preprocessing...", extra={'step': 'Main'})
    meta_info, filter_summary, adata = preprocess(adata)

    #### Calculate marker genes
    logger.info("Calculate marker genes...", extra={'step': 'Main'})
    adata = calculate_markers(adata,key_celltype=args.key_celltype)

    #### Save results
    # meta info (cell number & gene number)
    with open(os.path.join(args.outDir,'out','json','meta_info.preprocessing.json'),'w') as f:
        json.dump(meta_info,f,cls=NpEncoder)
    # filter summary
    with open(os.path.join(args.outDir,'out','json','filter_summary.preprocessing.json'),'w') as f:
        json.dump(filter_summary,f,cls=NpEncoder)
    # h5ad
    adata.write_h5ad(os.path.join(args.outDir,'sc.h5ad'))
    # marker gene
    de_key = 'rank_genes_groups' + '_' + args.key_celltype
    markers = extract_markers(adata,de_key=de_key)
    with open(os.path.join(args.outDir,'out','json','sc_markers.json'),'w') as f:
        json.dump(markers,f)

    logger.info("All Finished!", extra={'step': 'Main'})

if __name__ == "__main__":
    import sys

    main(sys.argv)