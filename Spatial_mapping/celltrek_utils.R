read_sc <- function(sc_h5ad, key_celltype = 'cell_type') {
  
  
  sc_h5ad = anndata::read_h5ad(sc_h5ad)
  
  sc = CreateSeuratObject(counts = t(as.matrix(sc_h5ad$X)))
  
  if (! all(rownames(sc_h5ad$obs) == colnames(sc))) {
    
    stop("The cell names in count matrix and metadata are not consistent!")
    
  }
  
  if (! key_celltype %in% names(sc_h5ad$obs)) {
    
    stop("Column name of cell type are not found in metadata!")
    
  } 
  
  sc@meta.data$cell_type = sc_h5ad$obs[,key_celltype]
  
  return(sc)
  
}

read_st <- function(st_h5ad) {
  
  # To be revised
  
  st_h5ad = anndata::read_h5ad(st_h5ad)
  
  # CreateSeuratObject
  
  if (is.null(st_h5ad$raw)) {
    st = CreateSeuratObject(counts = t(as.matrix(st_h5ad$X)),project = 'st',assay = 'Spatial')
  } else{
    cts = as.matrix(st_h5ad$raw$X)
    rownames(cts) = st_h5ad$obs_names
    colnames(cts) = st_h5ad$var_names
    cts = t(cts)
    st = CreateSeuratObject(counts = cts,project = 'st',assay = 'Spatial')
  }
  
  # Add Coordinate
  
  coord = data.frame(st_h5ad$obsm$spatial,row.names = st_h5ad$obs_names)
  colnames(coord) = c('x','y')
  
  st@images$image =  new(
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "st_",
    coordinates = coord
  )
  
  # Add coordinate to meta.data
  
  if(all(rownames(st@meta.data) == rownames(st@images$image@coordinates))){
    st@meta.data$x = st@images$image@coordinates$x
    st@meta.data$y = st@images$image@coordinates$y
  } else {
    
    stop("The meta data index and coordinate index are not consistent!")
  }
  
  return(st)
  
}

sc_analysis <- function(sc,nPCs=50) {

  # Basic scRNA-seq analysis
  # NOTE:
  # (1) Including scale, PCA (based on HVGs), UMAP using seurat
    
  sc <- ScaleData(sc, verbose = F)
  sc <- RunPCA(sc, npcs = nPCs, verbose = F)
  sc <- suppressWarnings(RunUMAP(sc, reduction = "pca", dims = 1:nPCs,verbose = F))
  tryCatch({sc <- RunTSNE(sc, reduction = "pca", dims = 1:nPCs,check_duplicates = FALSE)},
           error = function(e) {})
  return(sc)
  
}

traint <- function (st_data, sc_data, st_assay='Spatial', sc_assay='RNA', norm='LogNormalize', nfeatures=2000,
                    gene_kept=NULL, nPCs=50, ...) {
  
  st_data$id <- rownames(st_data@meta.data)
  sc_data$id <- rownames(sc_data@meta.data)
  st_data$type <- 'st'
  sc_data$type <- 'sc'
  
  DefaultAssay(st_data) <- st_assay
  DefaultAssay(sc_data) <- sc_assay
  
  #cat('Finding transfer anchors... \n')
  
  ## Integration features ##
  sc_st_list <- list(st_data=st_data, sc_data=sc_data)
  sc_st_features <- Seurat::SelectIntegrationFeatures(sc_st_list, nfeatures=nfeatures,verbose = F)
  
  if (!is.null(gene_kept)) {
    sc_st_features <- union(sc_st_features, gene_kept)
  }
  
  sc_st_features <- sc_st_features[(sc_st_features %in% rownames(st_data)) &
                                     (sc_st_features %in% rownames(sc_data))]
  #cat('Using', length(sc_st_features), 'features for integration... \n')
  
  ###
  sc_st_anchors <- Seurat::FindTransferAnchors(reference = sc_data, query = st_data,
                                               reference.assay = sc_assay, query.assay = st_assay,
                                               normalization.method = norm, features = sc_st_features, reduction = 'cca', verbose = F)
  
  #cat('Data transfering... \n')
  st_data_trans <- Seurat::TransferData(anchorset = sc_st_anchors,
                                        refdata = sc_data[[sc_assay]]@data[sc_st_features, ], weight.reduction = 'cca', verbose = F)
  st_data@assays$transfer <- st_data_trans
  
  #cat('Creating new Seurat object... \n')
  sc_st_meta <- dplyr::bind_rows(st_data@meta.data, sc_data@meta.data)
  counts_temp <- cbind(data.frame(st_data[['transfer']]@data), data.frame(sc_data[[sc_assay]]@data[sc_st_features, ] %>% data.frame))
  rownames(sc_st_meta) <- sc_st_meta$id
  colnames(counts_temp) <- sc_st_meta$id
  sc_st_int <- CreateSeuratObject(counts = counts_temp, assay = 'traint', meta.data = sc_st_meta)
  sc_st_int[['traint']]@data <- sc_st_int[['traint']]@counts
  sc_st_int[['traint']]@counts <- matrix(NA, nrow = 0, ncol = 0)
  
  #cat('Scaling -> PCA -> UMAP... \n')
  sc_st_int <- ScaleData(sc_st_int, features = sc_st_features, verbose = F) %>% RunPCA(features = sc_st_features, npcs = nPCs, verbose = F)
  sc_st_int <- RunUMAP(sc_st_int, dims = 1:nPCs, verbose = F)
  sc_st_int@images <- st_data@images
  sc_st_int@images[[1]]@coordinates <- data.frame(x=sc_st_int@meta.data$x,
                                                  y=sc_st_int@meta.data$y,
                                                  row.names = rownames(sc_st_int@meta.data))
  
  return (sc_st_int)
}

celltrek_dist <- function (st_sc_int, int_assay='traint', reduction='pca', intp = T, intp_pnt=5000, intp_lin=F, nPCs=50, ntree=1000, keep_model=T) {
  DefaultAssay(st_sc_int) <- int_assay
  spot_dis <- median(dbscan::kNN(na.omit(st_sc_int@meta.data[, c('x', 'y')]), k=6)$dist)
  #cat('Distance between spots is:', spot_dis, '\n')
  
  st_sc_int_pca <- st_sc_int@reductions[[reduction]]@cell.embeddings[, 1:nPCs] %>% data.frame %>%
    mutate(id=st_sc_int$id, type=st_sc_int$type, class=st_sc_int$cell_type,
           x=st_sc_int$x, y=st_sc_int$y)
  
  #cat('Random Forest training... \n')
  ## Training on ST ##
  st_pca <- st_sc_int_pca %>% dplyr::filter(type=='st') %>% dplyr::select(-c(id:class))
  rf_train <- randomForestSRC::rfsrc(Multivar(x, y) ~ ., st_pca, block.size=5, ntree=ntree)
  
  ## Interpolation ##
  ## Uniform sampling ##
  if (intp) {
    #cat ('Interpolating...\n')
    spot_ratio <- intp_pnt/nrow(st_pca)
    st_intp_df <- apply(st_pca[, c('x', 'y')], 1, function(row_x) {
      runif_test <- runif(1)
      if (runif_test < spot_ratio%%1) {
        theta <- runif(ceiling(spot_ratio), 0, 2*pi)
        alpha <- sqrt(runif(ceiling(spot_ratio), 0, 1))
        x <- row_x[1] + (spot_dis/2)*sin(theta)*alpha
        y <- row_x[2] + (spot_dis/2)*cos(theta)*alpha
      } else {
        theta <- runif(floor(spot_ratio), 0, 2*pi)
        alpha <- sqrt(runif(floor(spot_ratio), 0, 1))
        x <- row_x[1] + (spot_dis/2)*sin(theta)*alpha
        y <- row_x[2] + (spot_dis/2)*cos(theta)*alpha
      }
      data.frame(x, y)
    }) %>% Reduce(rbind, .)
    
    suppressWarnings(st_intp_df <- apply(st_pca[, 1:nPCs], 2, function(col_x) {
      akima::interpp(x=st_pca$x, y=st_pca$y, z=col_x,
                     linear=intp_lin, xo=st_intp_df$x, yo=st_intp_df$y) %>%
        magrittr::extract2('z')
    }) %>% data.frame(., id='X', type='st_intp', st_intp_df) %>% na.omit)
    st_intp_df$id <- make.names(st_intp_df$id, unique = T)
    st_sc_int_pca <- bind_rows(st_sc_int_pca, st_intp_df)
  }
  
  #cat('Calcluate median distance in between st spots (raw and interpp mixed)...  \n')
  spot_dis_intp = median(dbscan::kNN(st_sc_int_pca[st_sc_int_pca$type!='sc', c('x', 'y')], k=4)$dist)
  
  #cat('Random Forest prediction...  \n')
  ## Testing on all ##
  rf_pred <- randomForestSRC::predict.rfsrc(rf_train, newdata=st_sc_int_pca[, c(1:nPCs)], distance='all')
  
  #cat('Making distance matrix... \n')
  rf_pred_dist <- rf_pred$distance[st_sc_int_pca$type=='sc', st_sc_int_pca$type!='sc'] %>%
    magrittr::set_rownames(st_sc_int_pca$id[st_sc_int_pca$type=='sc']) %>% 
    magrittr::set_colnames(st_sc_int_pca$id[st_sc_int_pca$type!='sc'])
  output <- list()
  output$spot_dis <- spot_dis
  output$spot_dis_intp <- spot_dis_intp
  output$celltrek_dist <- rf_pred_dist
  output$coord_df <- st_sc_int_pca %>% magrittr::set_rownames(.$id) %>% dplyr::select(-id)
  if (keep_model) {
    output$model <- rf_train
  }
  return (output)
}

celltrek_chart <- function (dist_mat, coord_df, spot_dis_intp, dist_cut=500, n_spots=10,n_cells=10,n_redundancy=1) {
  
  # Not any more return pairwise distance 
  
  #cat('Making graph... \n')
  dist_mat[dist_mat>dist_cut] <- NA
  dist_mat_dt <- data.table::data.table(Var1=rownames(dist_mat), dist_mat)
  dist_edge_list <- data.table::melt(dist_mat_dt, id=1, na.rm=T)
  colnames(dist_edge_list) <- c('Var1', 'Var2', 'value')
  
  dist_edge_list$Var1 %<>% as.character # Var1 is cells in sc
  dist_edge_list$Var2 %<>% as.character # Var2 is spots in st
  dist_edge_list$Var1_type <- 'sc'
  dist_edge_list$Var2_type <- 'non-sc'
  
  #cat('Pruning graph...\n')
  dist_edge_list_sub <- suppressMessages(dplyr::inner_join(dist_edge_list %>% group_by(Var1) %>% slice_min(order_by = value,n = n_spots,with_ties = F),
                                          dist_edge_list %>% group_by(Var2) %>% slice_min(order_by = value,n = n_cells,with_ties = F))) %>% data.frame()
  dist_edge_list_sub = dist_edge_list_sub %>% group_by(Var1) %>% slice_min(order_by = value,n = n_redundancy,with_ties = F)
  
  #cat('Spatial Charting SC data...\n')
  
  sc_coord <-  data.frame(
    id_sc_new=make.names(dist_edge_list_sub$Var1,unique = T),
    id_sc=dist_edge_list_sub$Var1,
    id_st=dist_edge_list_sub$Var2,
    dist=dist_edge_list_sub$value,
    coord_df[match(dist_edge_list_sub$Var2, rownames(coord_df)),c('x','y')],
    row.names = NULL)
  
  ## Add noise ##
  theta <- runif(nrow(dist_edge_list_sub), 0, 2*pi)
  alpha <- sqrt(runif(nrow(dist_edge_list_sub), 0, 1))
  sc_coord$x_noise <- sc_coord$x + (spot_dis_intp/2)*sin(theta)*alpha
  sc_coord$y_noise <- sc_coord$y + (spot_dis_intp/2)*cos(theta)*alpha
  
  return(sc_coord)
}

celltrek_seurat <- function(sc_coord,sc_data,traint_res,sc_assay='RNA') {
  
  #sc_data=sc
  #sc_coord=sc_coord_list_1$sc_coord
  #sc_assay='RNA'
  #traint_res # for images with raw and interpp coordinates
  
  sc_data$id <- colnames(sc_data)
  
  # Create seurat object
  
  counts = sc_data[[sc_assay]]@data[, sc_coord$id_sc] %>% magrittr::set_colnames(sc_coord$id_sc_new)
  
  metadata = sc_data@meta.data[sc_coord$id_sc, ] %>% mutate(id_new = sc_coord$id_sc_new) %>% magrittr::set_rownames(sc_coord$id_sc_new)# id & id_new
  
  sc_out <- CreateSeuratObject(counts=counts,
                               project='celltrek', assay=sc_assay,
                               meta.data=metadata)
  
  sc_out@meta.data <- cbind(sc_out@meta.data, 
                            sc_coord[match(sc_out@meta.data$id_new,sc_coord$id_sc_new),]) %>% 
    select(-c(id_sc,id_sc_new))
  
  sc_out[[sc_assay]]@counts <- matrix(nrow = 0, ncol = 0)
  
  # Add noisy coord to dimension "spatial"
  
  coord_noise = sc_out@meta.data %>% dplyr::select(c(x_noise, y_noise)) %>% magrittr::set_colnames(c('spatial_1','spatial_2')) %>% as.matrix()
  sc_out@reductions$spatial <- CreateDimReducObject(embeddings=coord_noise, assay=sc_assay, key='spatial_')
  
  # Add image
  sc_out@images <- traint_res@images
  sc_out@images[[1]]@assay <- sc_assay
  sc_out@images[[1]]@coordinates <- data.frame(x_noise=sc_out@meta.data$x_noise, y_noise=sc_out@meta.data$y_noise) %>% 
    magrittr::set_rownames(colnames(sc_out)) # add noise coordinate
  
  return(sc_out)
  
}

.regularise_df <- function(df, drop_single_values = FALSE) {
  if (ncol(df) == 0) df[["name"]] <- rownames(df)
  if (drop_single_values) {
    k_singular <- sapply(df, function(x) length(unique(x)) == 1)
    if (sum(k_singular) > 0) {
      warning(
        paste("Dropping single category variables:"),
        paste(colnames(df)[k_singular], collapse = ", ")
      )
    }
    df <- df[, !k_singular, drop = F]
    if (ncol(df) == 0) df[["name"]] <- rownames(df)
  }
  return(df)
}

seurat2anndata <- function(obj, outFile = NULL, assay = "RNA", main_layer = "data", transfer_layers = NULL, drop_single_values = FALSE) {
  if (!requireNamespace("Seurat")) {
    stop("This function requires the 'Seurat' package.")
  }
  main_layer <- match.arg(main_layer, c("data", "counts", "scale.data"))
  transfer_layers <- transfer_layers[
    transfer_layers %in% c("data", "counts", "scale.data")
  ]
  transfer_layers <- transfer_layers[transfer_layers != main_layer]
  
  X <- Seurat::GetAssayData(object = obj, assay = assay, slot = main_layer)
  
  obs <- .regularise_df(obj@meta.data, drop_single_values = drop_single_values)
  
  var <- .regularise_df(Seurat::GetAssay(obj, assay = assay)@meta.features, drop_single_values = drop_single_values)
  
  obsm <- NULL
  reductions <- names(obj@reductions)
  if (length(reductions) > 0) {
    obsm <- sapply(
      reductions,
      function(name) as.matrix(Seurat::Embeddings(obj, reduction = name)),
      simplify = FALSE
    )
    names(obsm) <- sapply(tolower(names(obj@reductions)), function(name) {if(name != 'spatial') paste0("X_", tolower(name)) else {name}}) %>% unname()
  }
  
  layers <- list()
  for (layer in transfer_layers) {
    mat <- Seurat::GetAssayData(object = obj, assay = assay, slot = layer)
    if (all(dim(mat) == dim(X))) layers[[layer]] <- Matrix::t(mat)
  }
  
  #suppressMessages(anndata <- reticulate::import("anndata", convert = FALSE))
  
  adata <- anndata::AnnData(
    X = Matrix::t(X),
    obs = obs,
    var = var,
    obsm = obsm,
    layers = layers,
    dtype = 'float64'
  )
  
  if (!is.null(outFile)) {
    anndata::write_h5ad(adata, outFile)
  }
  
}

current_time <- function() {
  # Convert local time to UTC time
  return(format(lubridate::with_tz(Sys.time(),tzone = 'UTC'),"%H:%M:%S"))
  
}

# Convert UMAP dataframe to Json

umap2json <- function(umap,output) {
  
  umap$group = paste0(umap$modality,ifelse(is.na(umap$cell_type),"",paste0('-',umap$cell_type)))
  
  umap_json = list()
  
  for (grp in unique(umap$group)) {
    
    umap_json[[grp]][['id']] = apply(umap[umap$group == grp,],1, function(x)
      unname(x['id']),
      simplify = F)  %>% unname()
    
    umap_json[[grp]][['value']] = apply(umap[umap$group == grp,],1, function(x)
      as.numeric(unname(c(x['UMAP_1'],x['UMAP_2']))),
      simplify = F) %>% unname()
    
  }
  
  write(toJSON(umap_json),output)
  
}

# Convert RF-distance to Json

rfdist2json <- function(rfdist,output) {
  
  rfdist_json = list()
  
  rfdist_json[['x']] = rfdist$x
  rfdist_json[['density']] = rfdist$density
  
  write(toJSON(rfdist_json),output)
  
}

# UMAP filtering
reduction_filter <- function(traint,reduction='umap',cutoff=1){
  st_sc_reduction <- traint@reductions[[reduction]]@cell.embeddings
  st_reduction <- st_sc_reduction[traint$type == 'st',] %>% as.data.frame()
  sc_reduction <- st_sc_reduction[traint$type == 'sc',] %>% as.data.frame()
  
  num_cells=nrow(sc_reduction)
  #cat('Total number of single cells：' , num_cells, '\n')
  #cat("Number of ST spots: ", nrow(st_reduction),'\n')
  
  sc_reduction_scoped <- data.frame(dim1=c(),dim2=c())
  
  for (row_num in 1:nrow(st_reduction)){
    #print(row_num)
    umap1 <- as.numeric(st_reduction[row_num,1])
    umap2 <- as.numeric(st_reduction[row_num,2])
    umap1_scope <- c(umap1-cutoff,umap1+cutoff)
    umap2_scope <- c(umap2-cutoff,umap2+cutoff)
    
    sc_scoped <- sc_reduction[(sc_reduction[,1] > umap1_scope[1]) & (sc_reduction[,1] < umap1_scope[2]) & (sc_reduction[,2] > umap2_scope[1]) & (sc_reduction[,2] < umap2_scope[2]),] # dim1 

    if (nrow(sc_scoped) == 0) { 
      next
    }
    
    list_outline <- c()
    for (i in c(1:nrow(sc_scoped))) {
      #print(i)
      dist = sqrt((umap1-sc_scoped[i,1])^2+(umap2-sc_scoped[i,2])^2)
      #print(num)
      if (dist >= cutoff) {
        list_outline <- c(list_outline, -i)
      } 
    }
    sc_scoped <- sc_scoped[list_outline,]
    
    #plot(sc_scoped$UMAP_1,sc_scoped$UMAP_2)
    sc_reduction <- sc_reduction[!row.names(sc_reduction) %in% row.names(sc_scoped),]
    sc_reduction_scoped <- rbind(sc_reduction_scoped, sc_scoped)
  }

  num_filter = nrow(sc_reduction)
  #cat('Number of filtered single cells: ' , num_filter,'\n')
  #cat('Proportion of filered single cells: ', num_filter/num_cells, '\n')
  
  cells_rm <- row.names(sc_reduction)
  suppressWarnings(traint_filtered <- traint[,!colnames(traint) %in% cells_rm])
  return(traint_filtered)
}

knn_cutoff <- function(traint, reduction='umap', knn_num=10){
  ## Calculate knn of nearest knn_num points
  st_sc_reduction <- traint@reductions[[reduction]]@cell.embeddings
  st_reduction <- st_sc_reduction[traint$type == 'st',] %>% as.data.frame()
  KNN <- dbscan::kNN(st_reduction, k=knn_num)
  kNN_dist <- as.data.frame(KNN$dist)
  kNN_id <- KNN$id
  
  mean_list = rowMeans(kNN_dist)
  return(mean(mean_list))
}
