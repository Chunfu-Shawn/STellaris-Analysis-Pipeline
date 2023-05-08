import pandas as pd
import numpy as np
import scanpy as sc
from pandas import DataFrame
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import KMeans
import argparse
import os
from itertools import combinations_with_replacement
import scipy.stats as ss
import logging
import sys
import seaborn as sns
from plotnine import ggplot, aes, geom_boxplot
import json
import re


logger = logging.getLogger("calculate_euclidean_distance")
logger.setLevel(logging.INFO)
handler=logging.StreamHandler(sys.stdout)
handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s:%(step)s - %(levelname)s - %(message)s"))
logger.addHandler(handler)


def get_dis(list1,list2):
    l1 = len(list1)
    l2 = len(list2)
    list_a = list1 + list2
    result = squareform(pdist(list_a))
    result = result[0:l1, l1:(l1+l2)]
    result = pd.DataFrame(result)
    return result


def cal_dis(input_df):
    input_df = input_df.replace(0,np.nan)
    row_min = input_df.min(axis=1)
    row_avg = np.mean(row_min)
    col_min = input_df.min(axis=0)
    col_avg = np.mean(col_min)
    res = (row_avg + col_avg)/2
    if pd.isnull(res):  # In this case, all values in input_df was turned into np.nan. That means distance is 0.
        res = 0
    return res


def calculate_euclidean(coord,annotation,annotation_set):
    combs = [list(i) for i in list(combinations_with_replacement(annotation_set, 2))]
    for i, p in enumerate(combs):
        coord_1 = coord[annotation == p[0]].tolist()
        coord_2 = coord[annotation == p[1]].tolist()
        dist_per_cell: DataFrame = get_dis(coord_1, coord_2)
        combs[i].append(cal_dis(dist_per_cell))
    dist_pair = pd.DataFrame(combs, columns=['cell_type_1', 'cell_type_2', 'euclidean_distance'])
    dist_pair.sort_values(['cell_type_1', 'cell_type_2'], inplace=True)
    dist_mat = dist_pair.pivot(index="cell_type_1", columns="cell_type_2", values="euclidean_distance")
    for i in range(dist_mat.shape[0]):
        for j in range(dist_mat.shape[1]):
            if pd.isnull(dist_mat.iloc[i, j]):
                dist_mat.iloc[i, j] = dist_mat.iloc[j, i]

    return dist_pair, dist_mat


def kmeans_clustering(dist_pair):
    # Extract euclidean distance
    dist = np.array(dist_pair['euclidean_distance']).reshape(-1, 1)
    # Perform k means clustering
    km = KMeans(n_clusters=3, max_iter=1000).fit(dist)
    # Build label correspondence
    label_mapper = {1: 'near', 2: 'medium', 3: 'far'}
    dict_label = {i: label_mapper[j] for i, j in zip(range(len(km.cluster_centers_)), ss.rankdata(km.cluster_centers_))}
    # Assign label (near, medium or far)
    dist_pair['group'] = [dict_label[i] for i in km.labels_]

    return dist_pair


def save_pair_json(dist_pair,output_path):
    dist_pair_json = [{'cell_type_pair':'-'.join([c1,c2]),'euclidean_distance':dist, 'group':grp} \
                      for c1,c2,dist,grp in zip(dist_pair['cell_type_1'],
                                                dist_pair['cell_type_2'],
                                                dist_pair['euclidean_distance'],
                                                dist_pair['group'])]
    with open(output_path, "w") as f:
        json.dump(dist_pair_json, f)


def save_mat_json(dist_mat,dist_mat_heatmap,output_path):
    labels = [i for i in dist_mat_heatmap.ax_heatmap.yaxis.get_majorticklabels()]
    labels = [re.search(".*\'(.*)\'", str(l)).group(1) for l in labels]
    dist_mat_log2 = np.log2(dist_mat+1)
    dist_mat_log2 = dist_mat_log2[labels].reindex(labels)
    log_distance = [[i,j,dist_mat_log2.iloc[i,j]] for j in range(dist_mat_log2.shape[1]) for i in range(dist_mat_log2.shape[0])]
    dist_mat_json = {
        'cell_types':list(dist_mat_log2.index),
        'log_distance':log_distance
    }
    with open(output_path, "w") as f:
        json.dump(dist_mat_json, f)


def main(argsv):
    parser = argparse.ArgumentParser(description="Calculate the euclidean distance between each pair of cell types.")
    parser.add_argument("--sc_coordinate", type=str, help="scRNA-seq metadata with predicted spatial coordinates")
    parser.add_argument("--outDir", type=str, default="./", help="Path to output directory")
    args = parser.parse_args()

    os.makedirs(os.path.join(args.outDir,'out','table'),exist_ok=True)
    os.makedirs(os.path.join(args.outDir,'out','pdf'),exist_ok=True)
    os.makedirs(os.path.join(args.outDir,'out','json'),exist_ok=True)

    # Read data
    logger.info("Reading data...", extra={'step': 'Main'})
    sc_coordinate = pd.read_csv(args.sc_coordinate,header=0,index_col=False)
    coord = np.array(sc_coordinate[['x_noise','y_noise']])
    annotation = sc_coordinate['cell_type'].astype(str)
    annotation_set = list(set(annotation))

    # Calculate euclidean distance
    logger.info("Calculating euclidean distance...", extra={'step': 'Main'})
    dist_pair, dist_mat = calculate_euclidean(coord=coord,
                                              annotation=annotation,
                                              annotation_set=annotation_set)

    # K-means clustering
    logger.info("Performing K-means clustering...", extra={'step': 'Main'})
    dist_pair = kmeans_clustering(dist_pair)

    # Plot
    dist_mat_heatmap = sns.clustermap(np.log2(dist_mat+1), row_cluster=True, col_cluster=True, cmap='YlOrRd_r')
    dist_pair['group'] = pd.Categorical(dist_pair['group'], categories=['near', 'medium', 'far'])
    dist_pair_boxplot = (ggplot(dist_pair,aes(x='group',y='euclidean_distance')) + geom_boxplot()).draw()
    dist_mat_heatmap.savefig(os.path.join(args.outDir, 'out', 'pdf', 'cell_types_eucDist.heatmap.pdf'))
    dist_pair_boxplot.savefig(os.path.join(args.outDir, 'out', 'pdf', 'cell_types_eucDist.boxplot.pdf'))

    # Save results
    logger.info("Saving...", extra={'step': 'Main'})
    # save files
    dist_pair.to_csv(os.path.join(args.outDir,'out','table','cell_types_eucDist_group.csv'),index=False,header=True)
    dist_mat.to_csv(os.path.join(args.outDir,'out','table','cell_types_eucDist.csv'),index=True,header=True)
    # save json
    save_pair_json(dist_pair=dist_pair,
                   output_path=os.path.join(args.outDir,'out','json','cell_types_eucDist_group.json'))
    save_mat_json(dist_mat=dist_mat,
                  dist_mat_heatmap=dist_mat_heatmap,
                  output_path=os.path.join(args.outDir, 'out', 'json', 'cell_types_eucDist.json'))


if __name__ == '__main__':
    import sys

    main(sys.argv)