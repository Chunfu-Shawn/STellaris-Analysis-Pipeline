################################################
#File Name: run.sh
#Author: Up Lee    
#Mail: uplee@pku.edu.cn
#Created Time: Tue 07 Mar 2023 03:45:31 PM CST
################################################

#!/bin/sh 

#### 2023-03-07 ####

################################
#### Prepare scripts
################################

# Move xiaocf's script to backup/
mkdir -p bakcup/
mv  ST_screening.sh backup # xiaocf's script

# Move latest script here

cp /home/user/data3/uplee/projects/spatialTransWeb/bin/ST_screening.sh ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/process_sc.py ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/MIA_align.py ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/MIA_utils.py ./

# Modify ST_screening.sh
scriptDir=/home/user/data2/rbase/STellaris/scripts/ST_screening

#### EO 2023-03-07 ####

#### 2023-03-17 ####

## Modify MIA_align.py

# Modified contents:

    parser.add_argument("--cluster_key", help="The key of ST clusters for identifying marker genes",default="leiden1",choices=["leiden0.5","leiden1","leiden1.5","leiden2","leiden0.5_STAGATE","leiden1_STAGATE","leiden1.5_STAGATE","leiden2_STAGATE"])

    #### Perform MIA
    logger.info("Performing MIA...",extra={'step':'Main'})
    MIA_dict = MIA_multiprocess(st_marker_dict=st_marker_dict,
                     sc_marker=sc_marker,
                     cluster_key=args.cluster_key,
                     n_threads=args.n_threads)

#### EO 2023-03-17 ####

#### 2023-03-21 ####

# Modify process_sc.py (debug3 & debug4)

cp /home/user/data3/uplee/projects/spatialTransWeb/manuscript/manuscript_submission/debug/scripts/process_sc.py ./

#### EO 2023-03-21 ####

#### 2023-03-22 ####

# Modify MIA_align.py (debug5)

cp /home/user/data3/uplee/projects/spatialTransWeb/manuscript/manuscript_submission/debug/scripts/MIA_align.py ./

# Modify process_sc.py (debug6)

cp /home/user/data3/uplee/projects/spatialTransWeb/manuscript/manuscript_submission/debug/scripts/process_sc.py ./

# Modify process_sc.py (debug10)

cp /home/user/data3/uplee/projects/spatialTransWeb/manuscript/manuscript_submission/debug/scripts/process_sc.py ./

#### EO 2023-03-22 ####

#### 2023-03-23 ####

# Modify process_sc.py (debug13)

cp /home/user/data3/uplee/projects/spatialTransWeb/manuscript/manuscript_submission/debug/scripts/process_sc.py ./

# Modify MIA_align.py & MIA_utils.py (debug14)

cp /home/user/data3/uplee/projects/spatialTransWeb/manuscript/manuscript_submission/debug/scripts/MIA_align.py ./
cp /home/user/data3/uplee/projects/spatialTransWeb/manuscript/manuscript_submission/debug/scripts/MIA_utils.py ./

#### EO 2023-03-23 ####


