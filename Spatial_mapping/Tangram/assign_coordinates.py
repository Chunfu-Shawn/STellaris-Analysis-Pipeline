import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors
import os
import json
import argparse


def make_elements_unique(lst):
    count_dict = {}
    unique_lst = []
    for elem in lst:
        if elem not in count_dict:
            count_dict[elem] = 0
            unique_lst.append(elem)
        else:
            count_dict[elem] += 1
            unique_lst.append(f"{elem}.{count_dict[elem]}")
    return unique_lst


def find_k_nearest_neighbors(data, k):
    """
    Finds the k-nearest neighbors for a 2D array.

    Parameters:
        data (numpy array): The input 2D array.
        k (int): The number of neighbors to return for each point.

    Returns:
        numpy array: A 2D array containing the distances of the k-nearest neighbors for each point.
    """
    # create the NearestNeighbors object
    nbrs = NearestNeighbors(n_neighbors=k+1, algorithm='auto').fit(data)

    # find the k-nearest neighbors for each data point
    distances, indices = nbrs.kneighbors(data)

    # return the indices of the k-nearest neighbors for each data point
    return distances[:, 1:]  # exclude the first column which corresponds to the point itself


def spot_dis_intp(data, k):
    distances = find_k_nearest_neighbors(data, k)
    return np.median(distances)


def add_coord_noise(data, k):
    # add noise to x and y coordinates
    theta = np.random.uniform(0, 2*np.pi, size=data.shape[0])
    alpha = np.sqrt(np.random.uniform(0, 1, size=data.shape[0]))
    dist_intp = spot_dis_intp(data, k)
    noise = np.empty(data.shape)
    noise[:,0] = data[:,0] + (dist_intp/2) * np.sin(theta) * alpha
    noise[:,1] = data[:,1] + (dist_intp/2) * np.cos(theta) * alpha
    return noise


def main(argsv):
    parser = argparse.ArgumentParser(description="Assign coordinate generated from CytoSPACE.")
    parser.add_argument("--sc_reduction", help="Path to sc_reduction.h5ad file")
    parser.add_argument("--assigned_locations", help="Path to assigned_locations.csv file")
    parser.add_argument("--outDir",type=str,help="Path to output directory",default="./")

    args = parser.parse_args()

    # Read data
    adata = sc.read_h5ad(args.sc_reduction)
    assigned_locations = pd.read_csv(args.assigned_locations)

    # Generate mapping_summary

    map_factor = adata.obs.loc[:,'cell_type'].value_counts().index
    map_summary_mapped = adata.obs.loc[adata.obs['id'].isin(assigned_locations['OriginalCID']),'cell_type'].value_counts()
    map_summary_mapped = map_summary_mapped.reindex(map_factor,fill_value=0)
    map_summary_failed = adata.obs.loc[~adata.obs['id'].isin(assigned_locations['OriginalCID']),'cell_type'].value_counts()
    map_summary_failed = map_summary_failed.reindex(map_factor,fill_value=0)
    map_summary = {'Cell_type': list(map_factor),
                   'Failed': list(map_summary_failed),
                   'Mapped':list(map_summary_mapped)}

    # Assign coordinates

    adata.obs['id_interm'] = adata.obs_names

    assigned_locations1 = assigned_locations[['OriginalCID','SpotID','row','col']].rename(columns={'OriginalCID':'id','SpotID':'id_st','row':'x','col':'y'})

    coordinates = pd.merge(adata.obs, assigned_locations1, on='id', how='inner')
    coordinates.index = make_elements_unique(coordinates['id_interm'])

    adata1 = adata[coordinates['id_interm']].copy()
    adata1.obs = coordinates
    adata1.obs['id_interm'] = adata1.obs_names
    adata1.obs.rename(columns={'id_interm':'id_new'},inplace=True)
    adata1.var = pd.DataFrame({'name':adata1.var_names},index=adata1.var_names)

    adata1.obsm = None
    adata1.obsm['spatial'] = add_coord_noise(np.array(adata1.obs[['x','y']]), 10)
    adata1.obs = adata1.obs.assign(x_noise=adata1.obsm['spatial'][:,0], y_noise=adata1.obsm['spatial'][:,1])

    # Save
    with open(os.path.join(args.outDir, 'out', 'json', 'mapping_summary.json'),'w') as f:
        json.dump(map_summary, f)
    adata1.write_h5ad(os.path.join(args.outDir, 'sc_registered.h5ad'))
    adata1.obs.to_csv(os.path.join(args.outDir, 'out', 'table', 'sc_coordinate.csv'),index=False)
    adata1.obs[['id_new', 'cell_type']].to_csv(os.path.join(args.outDir, 'meta.txt'), index=False, sep="\t")

if __name__ == "__main__":
    import sys

    main(sys.argv)