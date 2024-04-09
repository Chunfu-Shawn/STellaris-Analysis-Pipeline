#!/bin/sh 

usage(){
  echo "Usage: bash $(basename $0) --sc_h5ad sc.h5ad --key_celltype cell_type --dataset ST_dataset_id --section ST_section_id --outDir output_directory [-h]"
  echo "Author: Chunfu Shawn"
  echo "Description: This script assigns single cells to spatial locations based on ST data, and then portrays cell-cell interactions in the spatial context."
  echo "Date: 2023-4-22"
  echo "------------------------------------------------"
  echo "OPTIONS"
  echo -e "     --sc_h5ad \t\tPath to h5ad of scRNA-seq (sc.h5ad) [Required]"
  echo -e "     --key_celltype \t\tColumn name of cell type [default: cell_type]"
  echo -e "     --dataset \t\tST dataset ID chosen by user for spatial mapping [Required]"
  echo -e "     --section \t\tST section ID chosen by user for spatial mapping [Required]"
  echo -e "     --num_epochs \t\tNumber of epochs for Tangram training [default: 1000]"
  echo -e "     --fragment_file \t\tThe fragment file of other modality, e.g., single-cell epigenomics data. This is a bed-lik file with first 5 columns to be used. [Optional]"
  echo -e "     --peak_file \t\tThe peak regions to aggregate epigenomic signals, e.g., active promoter region. THis is a bed-like file with first 4 columns to be used. If this file is provided, the fragment file MUST be also provided. [Optional]"
  echo -e "     --genome \t\tGenome version of the additional modality. Should match specified species. [default: hg38]"
  echo -e "     --top_n \t\tTop n genomic regions to retain for visualization if peak file is not provided. [default: 10000]"
  echo -e "     --species \t\tQuery species, could be 'Mus musculus' or 'Homo sapiens'"  
  echo -e "     --n_threads \t\tNumber of threads available to perform random forest prediction [default: 30]"
  echo -e "     --outDir \t\tPath to output directory [Required]"
  echo -e "     -h|--help \t\tprint this help page"
  exit 1
}

## Default argument
key_celltype=cell_type
n_threads=30
genome=hg38
top_n=10000
num_epochs=1000

while [[ $# -gt 0 ]]; do
    case $1 in
        --sc_h5ad)            sc_h5ad=$2;shift;;
        --key_celltype)       key_celltype=$2;shift;;
        --dataset)            dataset=$2;shift;;
        --section)            section=$2;shift;;
        --num_epochs)         num_epochs=$2;shift;;
        --fragment_file)      fragment_file=$2;shift;;
        --peak_file)          peak_file=$2;shift;;
        --genome)             genome=$2;shift;;
        --top_n)              top_n=$2;shift;;
        --species)            species=$2;shift;;        
        --n_threads)          n_threads=$2;shift;;
        --outDir)             outDir=$2;shift;;
        -h)                   usage;exit 1;;
        --)                   shift; break;;
        *)                    usage; echo -e "\n[ERR] $(date -u) Unkonwn option: $1"; exit 1;;
   esac
    shift
done

## Check necessary arguments
if [ -z $sc_h5ad ];then
   echo "The sc.h5ad must be provided!" && usage
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

if [[ ( ! -z $peak_file  ) && ( -z $fragment_file ) ]];then
   echo "the fragment file MUST be also provided if peak file provided!" && usage
fi

## Print argument
#echo -e "*** Arguments"
#echo -e "sc_h5ad\t$sc_h5ad"
#echo -e "key_celltype\t$key_celltype"
#echo -e "dataset\t$dataset"
#echo -e "section\t$section"
#echo -e "knn_num\t$knn_num"
#echo -e "n_spots\t$n_spots"
#echo -e "n_cells\t$n_cells"
#echo -e "n_redundancy\t$n_redundancy"
#echo -e "species\t$species"
#echo -e "n_threads\t$n_threads"
#echo -e "outDir\t$outDir"

## Configuration
scriptDir=/home/user/data2/rbase/STellaris/scripts/Spatial_mapping/Tangram
st_h5ad=/home/user/data/qijt/stdataset/h5ad_files_only_count/${dataset}/${section}/${section}.input.h5ad
source /home/user/BGM/rbase/env/anaconda3/bin/activate st_ann_anls
mkdir -p $outDir/log $outDir/out/table $outDir/out/json $outDir/out/pdf 

echo -e "*** Execution"

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

## Run Tangram
echo -e "`date -u`\tRun spatial mapping using Tangram..."

python $scriptDir/run_tangram.py \
 --sc_h5ad=$sc_h5ad \
 --st_h5ad=$st_h5ad \
 --key_celltype=$key_celltype \
 --num_epochs=$num_epochs\
 --outDir=$outDir 1>$outDir/log/run_tangram.log 2>&1

if [ ! $? -eq 0  ];then
    echo -e "`date -u`\tRun Tangram failed for some reason, please check log files!" >&2
    exit 1
else
    echo -e "`date -u`\tSuccessfully performed spatial mapping using Tangram!"
fi

## Convert celltrek results to jsonl for visualization
echo -e "`date -u`\tPrepare visualization..."
 
# Original scRNA-seq
echo -e "`date -u +'%H:%M:%S'` Original cells..."

python $scriptDir/prepare_data.py \
 --dataset $outDir/sc_reduction.h5ad \
 --name sc_reduction \
 --group $key_celltype \
 --outDir $outDir 1>$outDir/log/prepare_data.sc_reduction.log

if [ ! $? -eq 0  ];then
    echo -e "`date -u`\tVisualization preparation of original scRNA-seq failed for some reason, please check log files!" >&2
    exit 1
else
    echo -e "`date -u`\tSuccessfully performed visualization preparation of original scRNA-seq!"
fi

# Registered cells
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

## Process single-cell multiomics data if provided

if [ ! -z $fragment_file ];then

 echo -e "`date -u`\tProcess single-cell multiomics data..."
 
 if [ ! -z $peak_file ];then

  python $scriptDir/construct_multiomics.py \
   --fragment_file $fragment_file \
   --genome $genome \
   --peak_file $peak_file \
   --top_n $top_n \
   --sc_registered $outDir/sc_registered.h5ad \
   --outDir $outDir 1>$outDir/log/construct_multiomics.log
 
 else

  python $scriptDir/construct_multiomics.py \
   --fragment_file $fragment_file \
   --genome $genome \
   --top_n $top_n \
   --sc_registered $outDir/sc_registered.h5ad \
   --outDir $outDir 1>$outDir/log/construct_multiomics.log

 fi

if [ ! $? -eq 0  ];then
    echo -e "`date -u`\tProcessing single-cell multiomics data failed for some reason, please check log files!" >&2
    exit 1
else
    echo -e "`date -u`\tSuccessfully processed single-cell multiomics data!"
fi

 echo -e "`date -u`\tPrepare visualization for single-cell multiomics data..."
 python $scriptDir/prepare_data.sc_registered.py \
  --dataset $dataset \
  --section $section \
  --sc_registered $outDir/sc_multiomics.h5ad \
  --name sc_multiomics \
  --outDir $outDir 1>$outDir/log/prepare_data.sc_multiomics.log

fi

if [ ! $? -eq 0  ];then
    echo -e "`date -u`\tVisualization preparation of single-cell multiomics data failed for some reason, please check log files!" >&2
    exit 1
else
    echo -e "`date -u`\tSuccessfully performed visualization preparation of single-cell multiomics data!"
fi

## Calculate Euclidean distance
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
