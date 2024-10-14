#### Usage
usage = function(){
  cat(paste0("Usage: ", scriptName) )
  cat(" -sc_h5ad=sc_h5ad -key_celltype='cell_type' -dataset=dataset -section=section -n_threads=30 -outDir=outDir
Author: Uplee
Date: 2022-11-23
Option:
  -sc_h5ad       FILE       sc.h5ad generated from process_sc.py
  -key_celltype  STR        Column of cell type, default: 'cell_type'
  -dataset       STR        ST dataset id selected by the user
  -section       STR        ST section id selected by the user
  -knn_num       INT        Number of nearest neighboring cells to determine coembedding filtering cutoff, default: 50
  -n_threads     INT        Number of threads available for celltrek, default: 30
  -outDir        DIR        Output directory
")
  q(save = 'no')
}

#### Environment configuration
args <- commandArgs()
scriptPath = strsplit(args[4], '=', fixed = T)[[1]][2]
scriptName = basename(scriptPath)
scriptDir = dirname(scriptPath)
args = args[-(1:5)]
source("/home/user/data2/uplee/bin/R/common.R")
source("/home/user/data2/uplee/bin/R/ggplot.R")
source("/home/user/data3/uplee/projects/spatialTransWeb/manuscript/manuscript_revision/scripts/integration/Spatial_mapping_cytospace/celltrek_utils.R")

#### Default necessary arguments
key_celltype='cell_type'
knn_num=0
#### Argument passing

if(length(args) >= 1){
  for(i in 1:length(args)){
    arg = args[i]
    if(arg == '-h' || arg == '-help') usage()
    
    tmp = parseArg(arg, 'sc_h5ad', 'sc_h5ad'); if(!is.null(tmp)) sc_h5ad  = tmp
    tmp = parseArg(arg, 'key_celltype', 'key_celltype'); if(!is.null(tmp)) key_celltype  = tmp
    tmp = parseArg(arg, 'dataset', 'dataset'); if(!is.null(tmp)) dataset  = tmp
    tmp = parseArg(arg, 'section', 'section'); if(!is.null(tmp)) section  = tmp
    tmp = parseArgAsNum(arg, 'knn_num', 'knn_num'); if(!is.null(tmp)) knn_num  = tmp
    tmp = parseArgAsNum(arg, 'n_threads', 'n_threads'); if(!is.null(tmp)) n_threads  = tmp
    tmp = parseArg(arg, 'outDir', 'outDir'); if(!is.null(tmp)) outDir  = tmp
    
  }
}

#### Require library
suppressPackageStartupMessages(library(spatstat.core))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Matrix,quietly = T))
suppressPackageStartupMessages(library(magrittr,quietly = T))
suppressPackageStartupMessages(library(dbscan,quietly = T))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(akima))
suppressPackageStartupMessages(library(randomForestSRC))
suppressPackageStartupMessages(library(anndata))
suppressPackageStartupMessages(library(rjson))

#### Configuration
options(rf.cores=n_threads,mc.cores=n_threads)
st_dir='/home/user/data/qijt/stdataset/h5ad_files_only_count'
st_h5ad=file.path(st_dir,dataset,section,paste0(section,".input.h5ad"))

#### Read data
cat(current_time(), "Reading data\n")
#cat(" ###", current_time(), "sc\n")
sc = read_sc(sc_h5ad = sc_h5ad,key_celltype = key_celltype)
#cat(" ###", current_time(),"st\n")
st = read_st(st_h5ad = st_h5ad)
# Restrict nPCs to be strictly less than number of cells or features
nPCs = min(dim(sc)-1,50)

#### Preprocessing
cat(current_time(), "Preprocessing\n")
# sc
#cat(" ###", current_time(), "sc\n")
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "(?i)^MT-")
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000,verbose = F)
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000,verbose = F)
sc <- sc_analysis(sc, nPCs = nPCs)

# st
#cat(" ###", current_time(),"st\n")
st[["percent.mt"]] <- PercentageFeatureSet(st, pattern = "(?i)^MT-")
st <- NormalizeData(st, normalization.method = "LogNormalize", scale.factor = 10000,verbose = F)
st <- FindVariableFeatures(st, selection.method = "vst", nfeatures = 2000,verbose = F)

#### Coembedding
#cat("####",current_time(),"Running celltrek\n")
# Rename the cells/spots with syntactically valid names
st <- RenameCells(st, new.names=make.names(Cells(st),unique = T))
sc$id = colnames(sc)
sc <- RenameCells(sc, new.names=make.names(Cells(sc),unique = T))

# Integrate sc and st
cat(current_time(), "Integrating sc and st\n")
traint_res <- traint(st_data=st, 
                     sc_data=sc, 
                     st_assay = 'Spatial',
                     sc_assay='RNA',
                     norm = 'LogNormalize',
                     nPCs = nPCs)

#### UMAP filtering
cat(current_time(), "Coembedding filtering\n")

if(knn_num>0){
  
  cat("Number of nearest neighboring cells to determine filtering cutoff: " ,knn_num, "\n")
  dist_cutoff = knn_cutoff(traint_res,reduction = 'umap',knn_num = knn_num)
  traint_res1 <- reduction_filter(traint_res,reduction = 'umap',cutoff = dist_cutoff)
  
}else{
  
  cat("Skip coembedding filtering\n")
  traint_res1 = traint_res
  
}

#### Celltrek result summary

#cat("####", current_time(), "Summarizing results\n")

# NOTE:
# (1) umap: summary of mapping results
# (2) sc_coord: mapping results, including the spot that cell was assigned, RF distance between spot and cell.

# Filtering summary

filter_summary = list()
# Summary of cells before filtering
num_raw = sort(table(traint_res@meta.data[traint_res$type == 'sc','cell_type']),decreasing = T)
# Summary of cells after filtering
num_kept = table(factor(traint_res1@meta.data[traint_res1$type == 'sc','cell_type'],levels = names(num_raw)))
# Summary of cells removed
num_filter = num_raw - num_kept
# Save summary
filter_summary$Cell_type = names(num_raw)
filter_summary$Poor_coembedding = as.vector(num_filter)
filter_summary$Passed = as.vector(num_kept)

# UMAP (before and after filtering)

umap = data.frame(traint_res@reductions$umap@cell.embeddings,
                  modality = traint_res$type,
                  cell_type = traint_res@meta.data[,'cell_type'],
                  row.names = colnames(traint_res)) %>% rownames_to_column('id')

umap1 = data.frame(traint_res1@reductions$umap@cell.embeddings,
                   modality = traint_res1$type,
                   cell_type = traint_res1@meta.data[,'cell_type'],
                   row.names = colnames(traint_res1)) %>% rownames_to_column('id')
#### Save 
## Summary

# filtering summary (before and after filtering)
write(toJSON(filter_summary),file.path(outDir,'out','json','filter_summary.filtering.json'))

# umap (before and after filtering)
umap2json(umap = umap,output = file.path(outDir,'out','json','umap.preprocessing.json'))
umap2json(umap = umap1,output = file.path(outDir,'out','json','umap.filtering.json'))

## h5ad
seurat2anndata(obj = sc,
               outFile = file.path(outDir,"sc_reduction.h5ad"),
               assay = "RNA",
               main_layer = "data",
               transfer_layers = "scale.data")
