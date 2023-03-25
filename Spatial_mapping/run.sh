################################################
#File Name: run.sh
#Author: Up Lee    
#Mail: uplee@pku.edu.cn
#Created Time: Tue 07 Mar 2023 03:35:07 PM CST
################################################

#!/bin/sh 

#### 2023-03-07 ####

############################
#### Prepare scripts
############################

# Move latest script here and modify them.

cp /home/user/data3/uplee/projects/spatialTransWeb/bin/calculate_euclidean_distance.py ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/celltrek_utils.R ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/compress_results.sh ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/construct_multiomics.py ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/nicheAnchor-cellInteraction.sh ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/orthologs_mouse2human.json ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/plot_interaction.R ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/prepare_data.py ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/prepare_data.sc_registered.py ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/prepare_utils.py ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/run_celltrek.R ./

# Modify nicheAnchor-cellInteraction.sh
scriptDir=/home/user/data2/rbase/STellaris/scripts/Spatial_mapping

# Modify cellphonedb script
vi /home/user/BGM/uplee/anaconda3/envs/cellphonedb/lib/python3.9/site-packages/cellphonedb/utils/utils.py

def _read_h5ad(path: str, species: str = 'Homo sapiens') -> pd.DataFrame:
    adata = read_h5ad(path)
    df = adata.to_df().T

    if species == 'Mus musculus':
        app_logger.info('We detected your query species is mouse, performing orthologous gene name conversion...')
        import json
        with open('/home/user/data2/rbase/STellaris/scripts/Spatial_mapping/orthologs_mouse2human.json','r') as f:
            ortho_dict = json.load(f)
        genes = [ortho_dict.get(gene, gene) for gene in df.index.values]
        df.index = genes

    return df

# Rename script
mv nicheAnchor-cellInteraction.sh  spatial_mapping.sh

#### EO 2023-03-07 ####

#### 2023-03-17 ####

## Change calculate_euclidean_distance.py 

# Solve issue about Euclidean distance equal to 0

cp /home/user/data3/uplee/projects/spatialTransWeb/manuscript/manuscript_submission/debug/scripts/calculate_euclidean_distance.py ./

## Change compress_results.sh

# 解决压缩存在多级目录的问题

cp /home/user/data3/uplee/projects/spatialTransWeb/manuscript/manuscript_submission/debug/scripts/compress_results.sh ./

#### EO 2023-03-17 ####


#### 2023-03-22 ####

# Change calculate_euclidean_distance.py (debug7)

cp /home/user/data3/uplee/projects/spatialTransWeb/manuscript/manuscript_submission/debug/scripts/calculate_euclidean_distance.py ./

# Change celltrek_utils.R & celltrek_utils.R & spatial_mapping.sh (debug9)

cp /home/user/BGM/uplee/R_projects/debug_STellaris/scripts/celltrek_utils.R ./
cp /home/user/BGM/uplee/R_projects/debug_STellaris/scripts/run_celltrek.R ./
cp /home/user/BGM/uplee/R_projects/debug_STellaris/scripts/spatial_mapping.sh ./

#### EO 2023-03-22 ####

#### 2023-03-23 ####

# Change calculate_euclidean_distance.py (debug11)

cp /home/user/data3/uplee/projects/spatialTransWeb/manuscript/manuscript_submission/debug/scripts/calculate_euclidean_distance.py ./

#### 2023-03-23 ####
