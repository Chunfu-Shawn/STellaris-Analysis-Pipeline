# author: "Chunfu-Shawn"
# revise: "Juntian-Qi"
# date: "2022-11-21"

  
# install essential packages
library(ggplot2)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(viridis))
library(latex2exp)
suppressPackageStartupMessages(library(ComplexHeatmap))
library(optparse)
library(rjson)
suppressPackageStartupMessages(library(jsonlite))
library(readr)


# select cell pairs from microenvironment
select_cell_pairs_from_group = function(eucDist_group_path){
  print(" --- select cell pairs from eucDist group --- ")
  group_data = read_csv(eucDist_group_path,col_names = T, show_col_types = F)
  groups = unique(group_data[["group"]])
  groups_cell_pairs = list()
  for (group in groups){
    group_data_tmp = group_data[group_data['group'] == group,]
    for( i in 1:nrow(group_data_tmp) ){
      c1 = group_data_tmp[i,"cell_type_1"]
      c2 = group_data_tmp[i,"cell_type_2"]
      group1 = paste(group,c1,sep = '|')
      group2 = paste(group,c2,sep = '|')
      groups_cell_pairs[[group1]] = 
        c(groups_cell_pairs[[group1]], paste(c1,c2,sep = '|'))
      groups_cell_pairs[[group1]] = 
        c(groups_cell_pairs[[group1]], paste(c2,c1,sep = '|'))
      groups_cell_pairs[[group2]] = 
        c(groups_cell_pairs[[group2]], paste(c1,c2,sep = '|'))
      groups_cell_pairs[[group2]] = 
        c(groups_cell_pairs[[group2]], paste(c2,c1,sep = '|'))
    }
  }
  return(groups_cell_pairs)
}


# ggplot dotplot drawing
dotplot = function(data_ip,group=NULL,output_path,output_name){
  #adjusted parameters
  xlen = length(unique(data_ip$interacting_pair))
  ylen = length(unique(data_ip$cell_pairs))
  width = xlen/15+10
  xsize = -0.02*xlen + 9
  height = ylen/12+6
  ysize = -0.08*ylen + 10
  colorlimit = quantile(data_ip$means,probs = c(0,1))
  # draw
  p = ggplot(data=data_ip, aes(x=interacting_pair,y=cell_pairs)) +
    geom_point(aes(colour=means)) +
    scale_size_continuous(name = TeX('$-log_{10}$ p-value'), range = c(0.1,2))+ 
    scale_colour_gradientn(colours = viridis(100,option = "D"),
                           values = c(seq(0,0.3,length.out = 70),seq(0.3,1,length.out = 30)),
                           name = TeX('$log_{2}$ mean expr (molecule 1,molecule 2)')) +
    xlab('')+ylab('')+
    theme(panel.background = element_blank(),
          axis.line.x = element_line(),
          axis.line.y = element_line(),
          axis.text.x = element_text(angle = 45,hjust = 1,size=xsize),
          axis.text.y = element_text(angle = 30,hjust = 1,size=ysize),
          legend.key = element_blank(),
          legend.position = "bottom"
    )
  pdf(file.path(output_path, "out", "pdf", paste0(output_name,'.pdf')), width = width, height = height)
    print(p)
  dev.off()
  return(data_ip)
}

# dotplot save as json
dot_save_as_json = function(data_ip,group=NULL,output_path,output_name){
  # save as json
  group_inter = list()
  interacting_pairs = c(unique(arrange(data_ip,interacting_pair)$interacting_pair))
  cell_pairs = c(unique(arrange(data_ip,cell_pairs)$cell_pairs))
  group_inter[["yAxis"]] = interacting_pairs
  group_inter[["xAxis"]] = cell_pairs
  group_inter[["value"]] = apply(data_ip, 1, function(x){
    c(as.integer(which(cell_pairs==x["cell_pairs"]))-1,
      as.integer(which(interacting_pairs==x["interacting_pair"]))-1,
      x["means"])
  },simplify = F)
  return(group_inter)
}

# ligands-receptors pair save as json
lr_pair_save_as_json = function(data_ip,output_path){
  # save as json
  out_json = list()
  cell_pairs = unique(data_ip$cell_pairs)
  for ( cp in cell_pairs ){
    out_json[[cp]] = data_ip %>%
      filter(cell_pairs == cp) %>%
      apply(1, function(x){
        if(x["receptor_a"]==F){
          return(
            list(partner_l=x["partner_a"],
                 partner_r=x["partner_b"],
                 gene_l=x["gene_a"],
                 gene_r=x["gene_b"],
                 means=x["means"]
            ))
        }else{
          return(
            list(partner_r=x["partner_a"],
                 partner_l=x["partner_b"],
                 gene_r=x["gene_a"],
                 gene_l=x["gene_b"],
                 means=x["means"]
            ))
        }
      },simplify = F)
  }
  out_json = toJSON(out_json,pretty=F,auto_unbox=T)
  cat(out_json, file = paste(output_path,"/out/json",'/LR-pair_cell_pair.json',sep = ''))
}

# define function to draw interaction counts heatmap 
cells_interaction_heatmap = function(
    data_ip,eucDist_group_path,output_path,output_name="inter_count_heatmap"){
  
  eucDist_group_data = read_csv(eucDist_group_path,col_names = T,show_col_types = FALSE)
  cell_types = unique(c(eucDist_group_data[["cell_type_1"]],eucDist_group_data[["cell_type_1"]]))
  count_mat = matrix(data = 0,nrow = length(cell_types),ncol = length(cell_types))
  rownames(count_mat) = cell_types
  colnames(count_mat) = cell_types
  for (i in 1:length(rownames(data_ip))){
    cell_type = strsplit(data_ip[["cell_pairs"]][i], '\\|')[[1]]
    count_mat[cell_type[1],cell_type[2]] = count_mat[cell_type[1],cell_type[2]] +1
    if(cell_type[2]!=cell_type[1])
      count_mat[cell_type[2],cell_type[1]] = count_mat[cell_type[2],cell_type[1]] +1
  }
  
  # rescale color
  count_array = c(count_mat[lower.tri(count_mat, diag=T)]) # 取对角线
  # 判断count_array长度是否为1，或者都为一个值，导致quantile报错
  if(length(unique(count_array))==1){
    seq = quantile(c(0,1), probs = c(seq(0,0.6,length.out=10),seq(0.6,1,length.out=90)))
  }else{
    seq = quantile(count_array, probs = c(seq(0,0.6,length.out=10),seq(0.6,1,length.out=90)))
  }
  col_func = circlize::colorRamp2(seq, viridis(100,option = "D"))
  # draw
  pdf(file = paste(output_path,"/out/pdf/",output_name,'.pdf',sep = ''))
  hm = Heatmap(count_mat, name = 'Number of interaction',
               col = col_func,
               heatmap_legend_param = list(
                 direction = "horizontal"
               ),
               row_names_gp = gpar(fontsize=8),
               column_names_gp = gpar(fontsize=8),
               width = unit(10,"cm"), height = unit(10,"cm"))
  draw(hm , heatmap_legend_side="bottom")
  dev.off()
  
  # cluster
  hc = hclust(dist(count_mat),"ward.D")
  
  # save as json
  count_mat = count_mat[hc$order,hc$order]
  rownames(count_mat) = NULL
  colnames(count_mat) = NULL
  count_df = expand_grid(x=seq_len(nrow(count_mat)),y=seq_len(ncol(count_mat))) %>%
    rowwise() %>%
    mutate(count=count_mat[x,y])
  count = apply(count_df, 1, function(x){
    c(as.integer(x["x"])-1,as.integer(x["y"])-1,as.integer(x["count"]))
  },simplify = F)
  out_json =toJSON(list(cell_types=hc$labels[hc$order],count=count),pretty=F,auto_unbox=T)
  cat(out_json, file = paste(output_path,"/out/json/",output_name,'.json',sep = ''))
}


# data filter and plot
interaction_plot_save = function(means_path, pvalues_path,output_path,
                                 output_name="dot_plot", eucDist_group_path=NULL){
  means = read_tsv(means_path, col_names = T, show_col_types = FALSE)
  pvalues = read_tsv(pvalues_path, col_names = T, show_col_types = FALSE)
  print(" --- interacting ligands and receptors data filtering... --- ")
  
  # cluster by interacting pair
  clust = hclust(dist(pvalues[,-(1:11)] %>% as.matrix()))
  
  # wide table to long table and log2/10
  cell_pairs = colnames(means)
  means_all = means %>%
    mutate(row=row_number()) %>%
    pivot_longer(cols = colnames(means)[-(1:11)],names_to = "cell_pairs", values_to = "means") %>%
    select(interacting_pair,cell_pairs,means) %>% 
    mutate(interacting_pair = factor(interacting_pair, levels = unique(pvalues[clust$order,]$interacting_pair), ordered = T))
  means_all$means = log2(means_all$means)
  means_all$means[is.infinite(means_all$means)|means_all$means< -2] = -2 

  pvalues_all = pvalues %>%
    mutate(row=row_number()) %>%
    pivot_longer(cols = colnames(pvalues)[-(1:11)],names_to = "cell_pairs", values_to = "pvalue") %>%
    select(interacting_pair,cell_pairs,partner_a,partner_b,gene_a,gene_b,receptor_a,pvalue) %>% 
    mutate(interacting_pair = factor(interacting_pair, levels = unique(pvalues[clust$order,]$interacting_pair), ordered = T))


  # filter out pvalue > 0.01 and merge pvalues and means ( follow pvalue )
  pvalues_means_all = pvalues_all %>% filter(pvalue<= 0.01) %>% 
    left_join(means_all, by = c("cell_pairs","interacting_pair"))
  pvalues_means_all %>% 
    select(-pvalue) %>%
    pivot_wider(names_from = "cell_pairs", values_from = "means") %>% 
    write_csv(file = paste(output_path,'/out/table/significant_means_lt001.csv',sep = ''),
              na = "")
  
  # dot plot and save as json
  print(" --- Drawing p-value and mean expression data dotplot... --- ")
  print(" --- Saveing ligand-receptor pair as json format... --- ")
  pvalues_means_all %>% 
    dotplot(output_path=output_path,output_name=output_name) %>%
    lr_pair_save_as_json(output_path=output_path)
  
  if (!is.null(eucDist_group_path)){
    # interaction count heatmap
    print(" --- Drawing Number of interaction heatmap... --- ")
    cells_interaction_heatmap(data_ip = pvalues_means_all, 
                              eucDist_group_path=eucDist_group_path,
                              output_path=output_path)
    groups_cell_pairs = select_cell_pairs_from_group(
      eucDist_group_path=eucDist_group_path
    )
    ## for each environment
    print(" --- Drawing dotplot for each cell type in group and... --- ")
    print(" --- Saveing dotplot as json format... --- ")
    out_json = list()
    out_json[["group|cell_type"]] = names(groups_cell_pairs)
    for (group in names(groups_cell_pairs)){
      pvalues_means_group = pvalues_means_all[
        pvalues_means_all$cell_pairs %in% groups_cell_pairs[[group]],]
      out_json[[group]] = pvalues_means_group %>% 
        dotplot(output_path=output_path,output_name=paste(output_name,group,sep = "_")) %>%
        dot_save_as_json(group = group,output_path=output_path,output_name="inter_dotplot")
    }
    out_json =toJSON(out_json,pretty=F)
    cat(out_json, file = paste(output_path,"/out/json/",output_name,'.json',sep = ''))
  }
}


# operate shell parse
option_list = list(
  make_option(c("--work_dir","-w"), type = "character", default = ".", help="path of input files and output graphs")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
work_dir = opt$work_dir

# set paramters
output_path = file.path(work_dir)
means_path = file.path(work_dir,'out/table/means.txt')
pvalues_path = file.path(work_dir,'out/table/pvalues.txt')
eucDist_group_path = file.path(work_dir,'out/table/cell_types_eucDist_group.csv')


# main 
interaction_plot_save(
  means_path = means_path,
  pvalues_path = pvalues_path,
  eucDist_group_path = eucDist_group_path,
  output_path = output_path,
)
