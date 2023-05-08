import scanpy as sc
import snapatac2 as snap
import numpy as np
import pandas as pd
import os
import argparse
import re
import logging
import sys

logger = logging.getLogger("construct_multiomics")
logger.setLevel(logging.INFO)
handler=logging.StreamHandler(sys.stdout)
handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s:%(step)s - %(levelname)s - %(message)s"))
logger.addHandler(handler)

genome_dict = {
    'GRCh37': snap.genome.GRCh37,
    'GRCh38': snap.genome.GRCh38,
    'GRCm38': snap.genome.GRCm38,
    'GRCm39': snap.genome.GRCm39,
    'hg19': snap.genome.hg19,
    'hg38': snap.genome.hg38,
    'mm10': snap.genome.mm10,
    'mm39': snap.genome.mm39,

}

def read_fragment(fragment_file, genome, outDir):
    genome = genome_dict[genome]
    if os.path.exists(os.path.join(outDir, 'temp.h5ad')):
        os.remove(os.path.join(outDir, 'temp.h5ad'))

    if fragment_file.endswith('.bed.gz') or fragment_file.endswith('.bed'):
        adata = snap.pp.import_data(
            fragment_file,
            genome=genome,
            file=os.path.join(outDir, 'temp.h5ad'),
            sorted_by_barcode=False,
        )
    elif fragment_file.endswith('.bed.zip'):
        import zipfile
        fragment_dirname=os.path.dirname(fragment_file)
        fragment_file_unzip=os.path.splitext(fragment_file)[0]
        with zipfile.ZipFile(fragment_file, "r") as zip_ref:
            zip_ref.extractall(fragment_dirname)
        adata = snap.pp.import_data(
            fragment_file_unzip,
            genome=genome,
            file=os.path.join(outDir,'temp.h5ad'),
            sorted_by_barcode=False,
        )
    else:
        raise ValueError('Unsupported fragment file format! Please upload tab-delimited bed file compressed with gzip (.gz) or zip (.zip)!')

    return adata

def read_peak(peak_file):
    if peak_file.endswith('.bed.gz') or peak_file.endswith('.bed.zip') or peak_file.endswith('.bed'):
        peak = pd.read_csv(peak_file, index_col=None, header=None, sep='\t')
    else:
        raise ValueError('Unsupported peak file format! Please upload tab-delimited bed file compressed with gzip (.gz) or zip (.zip)!')

    if peak.shape[1] > 3:
        peak = peak.iloc[:,0:4]
        peak.columns = ['chr','start','end','id']
    elif peak.shape[1] == 3:
        peak.columns = ['chr','start','end']
        peak['id'] = [i + ':' + str(j) + '-' + str(z) for i, j, z in zip(peak['chr'], peak['start'], peak['end'])]
    else:
        raise ValueError('Bed file should contains no more than 3 columns! Please check your input format!')

    return peak

def snapatac_analysis(fragment_file, genome, outDir):
    adata = read_fragment(fragment_file,genome,outDir)

    snap.pp.add_tile_matrix(adata)
    snap.pp.select_features(adata)

    ## Doublet removal
    snap.pp.scrublet(adata)
    snap.pp.call_doublets(adata)
    snap.pl.scrublet(adata, interactive=False)
    adata.subset(~adata.obs["is_doublet"])

    return adata

def call_peak(adata, peak_file_path):
    peak = read_peak(peak_file_path)
    peak_list = [i + ':' + str(j) + '-' + str(z) for i, j, z in zip(peak['chr'], peak['start'], peak['end'])]
    peak_mat = snap.pp.make_peak_matrix(adata, use_rep=peak_list)
    peak_mat.obs_names = adata.obs_names
    peak_mat.var['peak'] = peak_mat.var.index
    peak_mat.var.index = peak['id']

    return peak_mat

def make_names(char_list):
    # Prepend "X" to string
    char_list = [ "X" + str(i) if bool(re.match(r'^\d+',i)) else i for i in char_list ]
    # Convert space to dot
    char_list = [ re.sub(r' ','.',i) for i in char_list ]
    # Convert special characters to dot
    special_characters='!"#$%&\'()*+,-./:;<=>?@[\\]^`{|}~'
    char_list = [ z.translate ({ord(c): "." for c in special_characters}) for z in char_list ]

    return char_list

def add_spatial(adata_omics,adata_sc_registered):
    # make name valid in a R manner
    adata_omics.obs_names = make_names(adata_omics.obs_names)
    # extract spatial coordinates (sc)
    columns = ['id', 'id_new', 'cell_type', 'x_noise', 'y_noise']
    spatial_info = adata_sc_registered.obs.loc[:, columns]
    # subset cells that exist in single-cell epigenomics data
    spatial_info = spatial_info.loc[spatial_info['id'].isin(adata_omics.obs_names), :]
    # Build spatial patterning of single-cell epigenomics data
    adata_omics1 = adata_omics[spatial_info['id'].tolist(), :].copy()
    adata_omics1.obs_names = spatial_info['id_new']
    adata_omics1.obs = pd.concat([adata_omics1.obs, spatial_info], axis=1)
    adata_omics1.obsm['spatial'] = np.array(adata_omics1.obs.loc[:, ['x_noise', 'y_noise']])

    return adata_omics1

def main(argsv):

    parser = argparse.ArgumentParser(description="Spatial construction of multiomics data.")
    parser.add_argument("--fragment_file", type=str, help="Input fragment file path")
    parser.add_argument("--genome", type=str, choices=["GRCh37", "GRCh38", "GRCm38", "GRCm39", "hg19", "hg38", "mm10", "mm39"], default='hg38', help="Genome version")
    parser.add_argument("--peak_file", type=str, help='The path to the selected genome region file')
    parser.add_argument("--top_n", type=int, default=10000, help="Top n genomic regions to retain for visualization if peak file is not provided")
    parser.add_argument("--sc_registered", type=str, help='The path to Celltrek result h5ad file')
    parser.add_argument("--outDir", type=str, default='./', help='The path of output files')
    args = parser.parse_args()

    # Basic processing
    logger.info("Perform basic SnapATAC2 analysis...",extra={'step':'Main'})
    adata = snapatac_analysis(fragment_file=args.fragment_file,
                              genome=args.genome,
                              outDir=args.outDir)

    # Aggregate peaks if provided
    logger.info("Aggregate peaks if provided...",extra={'step':'Main'})
    if args.peak_file:
        adata1 = call_peak(adata=adata, peak_file_path=args.peak_file)
    else:
        logger.info("Skip aggregating peaks, only keep first {} genomic regions...".format(args.top_n), extra={'step': 'Main'})
        adata1 = adata[:,0:args.top_n].copy()

    # Spatial construction by integrating registered scRNA-seq coordinates
    logger.info("Spatial construction of single-cell multiomics data...",extra={'step':'Main'})
    adata_sc = sc.read_h5ad(args.sc_registered)
    adata2 = add_spatial(adata_omics=adata1,
                         adata_sc_registered=adata_sc)

    # Normalization
    logger.info("Normalize...",extra={'step':'Main'})
    sc.pp.normalize_total(adata2)
    sc.pp.log1p(adata2)

    # Save
    logger.info("Save results...",extra={'step':'Main'})
    adata2.write_h5ad(os.path.join(args.outDir,'sc_multiomics.h5ad'))
    os.remove(os.path.join(args.outDir,'temp.h5ad'))

if __name__ == '__main__':
    import sys
    main(sys.argv)