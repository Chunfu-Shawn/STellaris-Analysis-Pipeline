#!/bin/sh 

usage(){
  echo "Usage: bash $(basename $0) --count countMatrix --label annotation --datasets ST_datasets --sections ST_sections --outDir output_directory [-h]"
  echo "Author: Kevin Lee"
  echo "Description: This script detects marker genes of scRNA-seq and perform MIA to find the most suitable ST section."
  echo "Date: 2022-11-21"
  echo "------------------------------------------------"
  echo "OPTIONS"
  echo -e "     --count \t\tCount matrix of scRNA-seq [Required]"
  echo -e "     --label \t\tMetadata including cell type of scRNA-seq [Required]"
  echo -e "     --key_celltype \t\tColumn name of cell type [default: cell_type]"
  echo -e "     --min_genes \t\tMinimum number of genes expressed required for a cell to pass filtering [default: 200]"
  echo -e "     --max_mt_pct \t\tMaximum percent of mitochondrial UMIs for a cell to pass filtering [default: 20]"
  echo -e "     --datasets \t\tA comma delimited list of ST datasets for MIA (allowed to be duplicated if one dataset containing mulitple sections) [Required]"
  echo -e "     --sections \t\tA comma delimited list of sections corresponding to ST datasets [Required]"
  echo -e "     --cluster_key \t\tThe key of ST clusters for identifying marker genes [default: leiden1] [Choice: leiden0.5, leiden1, leiden1.5, leiden2]"
  echo -e "     --st_dir \t\tThe directory to access ST marker genes [default: /home/user/data3/uplee/projects/spatialTransWeb/spatial/marker]"  
  echo -e "     --n_threads \t\tNumber of threads available to perform MIA [default: 30]"
  echo -e "     --outDir \t\tPath to output directory [Required]"
  echo -e "     -h|--help \t\tprint this help page"
  exit 1
}

## Default argument
key_celltype=cell_type
min_genes=200
max_mt_pct=20
cluster_key=leiden1
st_dir=/home/user/data3/uplee/projects/spatialTransWeb/spatial/marker
n_threads=30

while [[ $# -gt 0 ]]; do
    case $1 in
        --count)            count=$2;shift;;
        --label)            label=$2;shift;;
        --key_celltype)     key_celltype=$2;shift;;
        --min_genes)        min_genes=$2;shift;;
        --max_mt_pct)       max_mt_pct=$2;shift;;
        --datasets)          datasets=$2;shift;;
        --sections)          sections=$2;shift;;
        --cluster_key)      cluster_key=$2;shift;;
        --st_dir)           st_dir=$2;shift;;
        --n_threads)        n_threads=$2;shift;;
        --outDir)           outDir=$2;shift;;
        -h)                 usage;exit 1;;
        --)                 shift; break;;
        *)                  usage; echo -e "\n[ERR] $(date -u) Unkonwn option: $1"; exit 1;;
   esac
    shift
done

## Check mandatory arguments

if [ -z $count ];then
   echo "The count matrix of scRNA-seq must be provided!" && usage
fi

if [ -z $label ];then
   echo "The annotation of scRNA-seq must be provided!" && usage
fi

if [ -z $datasets ];then
   echo "No ST datasets were retrieved!" && usage
fi

if [ -z $sections ];then
   echo "No sections of ST datasets were retrieved!" && usage
fi

if [ -z $outDir ];then
   echo "Please specify output directory attributed to the current job!" && usage
fi

## Print argument
#echo -e "*** Arguments"
#echo -e "count\t$count"
#echo -e "label\t$label"
#echo -e "key_celltype\t$key_celltype"
#echo -e "min_genes\t$min_genes"
#echo -e "max_mt_pct\t$max_mt_pct"
#echo -e "datasets\t$datasets"
#echo -e "sections\t$sections"
#echo -e "cluster_key\t$cluster_key"
#echo -e "st_dir\t$st_dir"
#echo -e "n_threads\t$n_threads"
#echo -e "outDir\t$outDir"

## Configuration

scriptDir=/home/user/data2/rbase/STellaris/scripts/ST_screening
source /home/user/BGM/uplee/anaconda3/bin/activate spatialWeb
datasets=$( echo $datasets | sed 's@,@ @g' )
sections=$( echo $sections | sed 's@,@ @g' )

mkdir -p $outDir/log $outDir/out/table $outDir/out/pdf $outDir/out/json

echo -e "*** Execution"

## Check feasibility of ST sections
echo -e "`date -u`\tChecking feasibility of ST sections..."

arr_datasets=($datasets)
arr_sections=($sections)
if (("${#arr_datasets[@]}" != "${#arr_sections[@]}"));then
   echo -e "`date -u`\tNumber of datasets and sections are differ! Exit..." >&2
   exit 1
fi

datasets=""
sections=""

for idx in "${!arr_datasets[@]}";do
   marker_file=$st_dir/${arr_datasets[$idx]}/${arr_sections[$idx]}/${arr_sections[$idx]}.maker_gene.pkl
   if [ -f $marker_file ];then
      datasets=$( echo $datasets ${arr_datasets[$idx]} )
      sections=$( echo $sections ${arr_sections[$idx]} )
   fi
done

if [ -z "$datasets" ];then
   echo -e "`date -u`\tNo ST sections match your scRNA-seq data! Exit..." >&2
   exit 1
fi

## Process scRNA-seq
echo -e "`date -u`\tProcessing scRNA-seq..."

python $scriptDir/process_sc.py \
 --count $count \
 --label $label \
 --key_celltype $key_celltype \
 --min_genes $min_genes \
 --max_mt_pct $max_mt_pct \
 --outDir $outDir/ 1>$outDir/log/process_sc.log

if [ ! $? -eq 0  ];then
    echo -e "`date -u`\tscRNA-seq preprocessing failed for some reason, please check log files!" >&2
    exit 1
else
    echo -e "`date -u`\tSuccessfully performed scRNA-seq preprocessing!"
fi

## Perform MIA
echo -e "`date -u`\tPerforming MIA..."

python $scriptDir/MIA_align.py \
 --datasets $datasets \
 --sections $sections \
 --cluster_key $cluster_key \
 --st_dir $st_dir \
 --sc_path $outDir/out/json/sc_markers.json \
 --n_threads $n_threads \
 --outDir $outDir/ 1>$outDir/log/MIA_align.log

 if [ ! $? -eq 0  ];then
    echo -e "`date -u`\tMIA failed for some reason, please check log files!" >&2
    exit 1
else
    echo -e "`date -u`\tSuccessfully performed MIA!"
fi

echo -e "`date -u`\tAll finished!"
