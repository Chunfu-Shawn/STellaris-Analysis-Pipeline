#!/bin/sh 

usage(){
  echo "Usage: bash $(basename $0) --sc_h5ad sc.h5ad --key_celltype cell_type --dataset ST_dataset_id --section ST_section_id --outDir output_directory [-h]"
  echo "Author: Kevin Lee"
  echo "Description: This script assigns single cells to spatial locations based on ST data, and then portrays cell-cell interactions in the spatial context."
  echo "Date: 2022-11-28"
  echo "------------------------------------------------"
  echo "OPTIONS"
  echo -e "     --sc_h5ad \t\tPath to h5ad of scRNA-seq (sc.h5ad) [Required]"
  echo -e "     --sc_counts \t\tPath to counts of scRNA-seq for CytoSPACE (sc_counts.txt) [Required]"
  echo -e "     --sc_labels \t\tPath to labels of scRNA-seq for CytoSPACE (sc_labels.txt) [Required]"
  echo -e "     --dataset \t\tST dataset ID chosen by user for spatial mapping [Required]"
  echo -e "     --section \t\tST section ID chosen by user for spatial mapping [Required]"
  echo -e "     --key_celltype \t\tColumn name of cell type [default: cell_type]"
  echo -e "     --n_threads \t\tNumber of threads available to perform random forest prediction [default: 30]"
  echo -e "     --nosss \t\tNumber of selected subspots from ST data to limit the number of mapped cells for CytoSPACE [default: 5000]"
  echo -e "     --species \t\tQuery species, could be 'Mus musculus' or 'Homo sapiens' [Required]"  
  echo -e "     --outDir \t\tPath to output directory [Required]"
  echo -e "     -h|--help \t\tprint this help page"
  exit 1
}

#### Default argument
key_celltype=cell_type
n_threads=30
nosss=5000

while [[ $# -gt 0 ]]; do
    case $1 in
        --sc_h5ad)            sc_h5ad=$2;shift;;
        --sc_counts)          sc_counts=$2;shift;;
        --sc_labels)          sc_labels=$2;shift;;
        --dataset)            dataset=$2;shift;;
        --section)            section=$2;shift;;
        --key_celltype)       key_celltype=$2;shift;;
        --n_threads)          n_threads=$2;shift;;
        --nosss)              nosss=$2;shift;;
        --species)            species=$2;shift;;        
        --outDir)             outDir=$2;shift;;
        -h)                   usage;exit 1;;
        --)                   shift; break;;
        *)                    usage; echo -e "\n[ERR] $(date -u) Unkonwn option: $1"; exit 1;;
   esac
    shift
done

#### Check necessary arguments
if [ -z $sc_h5ad ];then
   echo "The sc.h5ad must be provided!" && usage
fi

if [ -z $sc_counts ];then
   echo "The counts of scRNA-seq for CytoSPACE must be provided!" && usage
fi

if [ -z $sc_labels ];then
   echo "The labels of scRNA-seq for CytoSPACE must be provided!" && usage
fi

if [ -z $dataset ];then
   echo "No ST datasets were provided!" && usage
fi

if [ -z $section ];then
   echo "No sections of ST datasets were retrieved!" && usage
fi

if [ -z "$species" ];then
   echo "Query species was not provided!" && usage
fi

if [ -z $outDir ];then
   echo "Please specify output directory attributed to the current job!" && usage
fi

#### Configuration
scriptDir=/home/user/data2/rbase/STellaris-Test/scripts/Spatial_mapping/CytoSPACE
source /home/user/BGM/uplee/anaconda3/bin/activate spatialWeb
mkdir -p $outDir/log $outDir/out/table $outDir/out/json $outDir/out/pdf $outDir/cytospace_results

echo -e "*** Execution"

#### Prepare ST data for CytoSPACE
echo -e "`date -u`\tPrepare ST data for CytoSPACE..."

python $scriptDir/prepare_st.py \
 --dataset $dataset \
 --section $section \
 --outDir $outDir 1>$outDir/log/prepare_st.log

if [ ! $? -eq 0  ];then
    echo -e "`date -u`\tPrepare ST data for CytoSPACE failed for some reason, please check log files!" >&2
    exit 1
else
    echo -e "`date -u`\tSuccessfully preparing ST data for CytoSPACE!"
fi

#### Coembedding filtering
echo -e "`date -u`\tPerforming coembedding filtering..."

Rscript $scriptDir/coembedding_filtering.R \
 -sc_h5ad=$sc_h5ad \
 -key_celltype=$key_celltype \
 -dataset=$dataset \
 -section=$section \
 -n_threads=$n_threads \
 -outDir=$outDir 1>$outDir/log/coembedding_filtering.log

if [ ! $? -eq 0  ];then
    echo -e "`date -u`\tPerforming coembedding filtering failed for some reason, please check log files!" >&2
    exit 1
else
    echo -e "`date -u`\tSuccessfully performing coembedding filtering!"
fi

#### Run CytoSPACE
echo -e "`date -u`\tRun CytoSPACE..."

source /home/user/BGM/uplee/anaconda3/bin/activate cytospace

cytospace \
   --scRNA-path $outDir/sc_counts.txt \
   --cell-type-path $outDir/sc_labels.txt \
   --st-path $outDir/st_counts.txt \
   --coordinates-path $outDir/st_coordinates.txt \
   -sss -nosss $nosss \
   --number-of-processors $n_threads \
   --output-folder $outDir/cytospace_results 1>$outDir/log/cytospace.log

if [ ! $? -eq 0  ];then
    echo -e "`date -u`\tRun CytoSPACE failed for some reason, please check log files!" >&2
    exit 1
else
    echo -e "`date -u`\tSuccessfully running CytoSPACE!"
fi

#### Assign CytoSPACE mapping coordinates
echo -e "`date -u`\tAssign CytoSPACE mapping coordinates..."

source /home/user/BGM/uplee/anaconda3/bin/activate spatialWeb

python $scriptDir/assign_coordinates.py \
 --sc_reduction $outDir/sc_reduction.h5ad \
 --assigned_locations $outDir/cytospace_results/assigned_locations.csv \
 --outDir $outDir 1>$outDir/log/assign_coordinates.log

if [ ! $? -eq 0  ];then
    echo -e "`date -u`\tAssign CytoSPACE mapping coordinates failed for some reason, please check log files!" >&2
    exit 1
else
    echo -e "`date -u`\tSuccessfully assigning CytoSPACE mapping coordinates!"
fi

#### Convert mapping results to jsonl for visualization

echo -e "`date -u`\tPrepare visualization..."
 
### Original scRNA-seq
echo -e "`date -u +'%H:%M:%S'` Original cells..."

python $scriptDir/prepare_data.py \
 --dataset $outDir/sc_reduction.h5ad \
 --name sc_reduction \
 --group cell_type \
 --outDir $outDir 1>$outDir/log/prepare_data.sc_reduction.log

if [ ! $? -eq 0  ];then
    echo -e "`date -u`\tVisualization preparation of original scRNA-seq failed for some reason, please check log files!" >&2
    exit 1
else
    echo -e "`date -u`\tSuccessfully performed visualization preparation of original scRNA-seq!"
fi

### Registered cells
echo -e "`date -u +'%H:%M:%S'` Registered cells..."

python $scriptDir/prepare_data.sc_registered.py \
 --dataset $dataset \
 --section $section \
 --sc_registered $outDir/sc_registered.h5ad \
 --name sc_registered \
 --outDir $outDir 1>$outDir/log/prepare_data.sc_registered.log

if [ ! $? -eq 0  ];then
    echo -e "`date -u`\tVisualization preparation of registered cells failed for some reason, please check log files!" >&2
    exit 1
else
    echo -e "`date -u`\tSuccessfully performed visualization preparation of registered cells!"
fi
 
#### Calculate Euclidean distance
echo -e "`date -u`\tCalculate Euclidean distance..."

sc_coordinate=$outDir/out/table/sc_coordinate.csv

python $scriptDir/calculate_euclidean_distance.py \
 --sc_coordinate $sc_coordinate \
 --outDir $outDir  1>$outDir/log/calculate_euclidean_distance.log

if [ ! $? -eq 0  ];then
    echo -e "`date -u`\tCalculating Euclidean distance failed for some reason, please check log files!" >&2
    exit 1
else
    echo -e "`date -u`\tSuccessfully performed calculating Euclidean distance!"
fi

## Run cellphonedb
echo -e "`date -u`\tDetect ligand-receptor interactions for each group of cell type pairs..."

source /home/user/BGM/uplee/anaconda3/bin/activate cellphonedb

meta=$outDir/meta.txt
counts=$outDir/sc_registered.h5ad
output_path_cpdb=$outDir/out/table

cellphonedb method statistical_analysis $meta $counts \
        --counts-data hgnc_symbol \
        --species "$species" \
        --threads $n_threads \
        --output-path $output_path_cpdb 1>$outDir/log/cellphonedb.log

if [ ! $? -eq 0  ];then
    echo -e "`date -u`\tDetecting cell-cell ligand-receptor interactions failed for some reason, please check log files!" >&2
    exit 1
else
    echo -e "`date -u`\tSuccessfully performed detecting cell-cell ligand-receptor interactions!"
fi

## Plot cell-cell interactions
echo -e "`date -u`\tPlot cell-cell interactions..."

Rscript $scriptDir/plot_interaction.R -w $outDir 1>$outDir/log/plot_interaction.log

if [ ! $? -eq 0  ];then
    echo -e "`date -u`\tPlotting cell-cell interactions failed for some reason, please check log files!" >&2
    exit 1
else
    echo -e "`date -u`\tSuccessfully performed plotting cell-cell interactions!"
fi

echo -e "`date -u`\tAll finished!"
