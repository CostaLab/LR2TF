if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(dorothea)
library(dplyr)
library(Seurat)
library(tibble)
library(pheatmap)
library(tidyr)
library(viper)
library(clustree)
library(stringr)

SaveTFList <- function(scores, name, out_path){
  highly_variable_tfs_all <- scores %>%
    group_by(tf) %>%
    mutate(var = var(avg))  %>%
    ungroup() %>%
    distinct(tf)
  
  ## We prepare the data for the plot
  summarized_viper_scores_df_all <- scores %>%
    semi_join(highly_variable_tfs_all, by = "tf") %>%
    dplyr::select(-std) %>%
    spread(tf, avg) %>%
    data.frame(row.names = 1, check.names = FALSE)
  tf_scores <- t(summarized_viper_scores_df_all)
  write.csv(tf_scores, file = paste0(out_path,'/tf_scores','_',name,'.csv'))
}

GetClusTree <- function(cluster_object, out_path, name, confidence){
  png(file = paste0(out_path,'/clustree','_',name,'_',confidence,'.png'), 
      width = 833, height = 550)
  print(clustree(cluster_object, prefix = "dorothea_snn_res."))
  dev.off()
}

ClusterByTFActivity <- function(seuratobject, resolutions, res, confidence, 
                                name, single_result_path){
  seuratobject <- RunPCA(seuratobject, features = rownames(seuratobject), 
                         verbose = FALSE)
  seuratobject <- FindNeighbors(seuratobject, dims = 1:10, verbose = FALSE)
  seuratobject <- FindClusters(seuratobject, resolution = resolutions, 
                               verbose = FALSE)
  
  GetClusTree(seuratobject, single_result_path, name, confidence)
  
  Idents(object = seuratobject) <- paste0("dorothea_snn_res.",res)
  seuratobject <- RunUMAP(seuratobject, dims = 1:10, umap.method = "uwot", 
                          metric = "cosine")
  
  png(file = paste0(single_result_path,
                    '/tf_activity_cluster',name,'_',confidence,'_',res,'.png'), 
      width = 833, height = 550)
  print(DimPlot(seuratobject, reduction = "umap", label = TRUE, 
                pt.size = 0.5, split.by = "protocol" ) + NoLegend())
  dev.off()
}

PlotTFActivity <- function(viper_scores, name, confidence, res, out_path){
  ## We select the 20 most variable TFs. (20*5 populations = 100)
  highly_variable_tfs <- viper_scores %>%
    group_by(tf) %>%
    mutate(var = var(avg))  %>%
    ungroup() %>%
    top_n(100, var) %>%
    distinct(tf)
  
  summarized_viper_scores_df <- viper_scores %>%
    semi_join(highly_variable_tfs, by = "tf") %>%
    dplyr::select(-std) %>%   
    spread(tf, avg) %>%
    data.frame(row.names = 1, check.names = FALSE) 
  
  palette_length = 100
  my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)
  
  my_breaks <- c(seq(min(summarized_viper_scores_df), 0, 
                     length.out = ceiling(palette_length/2) + 1),
                 seq(max(summarized_viper_scores_df)/palette_length, 
                     max(summarized_viper_scores_df), 
                     length.out = floor(palette_length/2)))
  
  png(file = paste0(out_path,
                    '/tf_activity_top20',name,'_',confidence,'_',res,'.png'),
      width = 833, height = 550)
  viper_hmap <- pheatmap(t(summarized_viper_scores_df),fontsize = 14, 
                         fontsize_row = 10, 
                         color = my_color, breaks = my_breaks, 
                         main = "DoRothEA Top 20 Variable TFs", angle_col = 45,
                         treeheight_col = 0,  border_color = NA)
  print(viper_hmap)
  dev.off()
}

DorotheaExecution <- function(seuratobject, name, out_path, 
                              resolutions, res, confidence){
  single_result_path <- paste0(out_path,'/',name,'_',res,'_',confidence,'_')
  dir.create(single_result_path)
  
  dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
  
  regulon <- dorothea_regulon_human %>%
    dplyr::filter(confidence %in% as.list(confidence))
  write.csv(regulon, paste0(out_path, "/dorothea_regulon_", confidence,".csv"))
  
  
  seuratobject <- run_viper(seuratobject, regulon,
                            options = list(method = "scale", minsize = 4, 
                                           eset.filter = FALSE, cores = 1, 
                                           verbose = FALSE))
  
  DefaultAssay(object = seuratobject) <- "dorothea"
  seuratobject <- ScaleData(seuratobject)
  
  ClusterByTFActivity(seuratobject, resolutions, res, confidence, 
                      name, single_result_path)
  
  Idents(object = seuratobject) <- "new_annotation"
  
  seuratobject.markers <- FindAllMarkers(seuratobject, only.pos = TRUE, 
                                         min.pct = 0, logfc.threshold = 0, 
                                         verbose = FALSE)
  write.csv(seuratobject.markers, file = 
              paste0(single_result_path,'/all_specificmarker_',
                     confidence, '_',name,'.csv'))
  
  seuratobject.markers <- FindAllMarkers(seuratobject, only.pos = TRUE, 
                                         min.pct = 0.25, logfc.threshold = 0.25, 
                                         verbose = FALSE)
  write.csv(seuratobject.markers, file = paste0(single_result_path,
                                                '/filtered_specificmarker_',
                                                confidence,'_',name,'.csv'))
  
  viper_scores_df <- GetAssayData(seuratobject, slot = "scale.data", 
                                  assay = "dorothea") %>%
    data.frame(check.names = F) %>%
    t()

  CellsClusters <- data.frame(cell = names(Idents(seuratobject)), 
                              cell_type = as.character(Idents(seuratobject)),
                              check.names = F)
  
  viper_scores_clusters <- viper_scores_df  %>%
    data.frame() %>% 
    rownames_to_column("cell") %>%
    gather(tf, activity, -cell) %>%
    inner_join(CellsClusters)
  
  summarized_viper_scores <- viper_scores_clusters %>% 
    group_by(tf, cell_type) %>%
    summarise(avg = mean(activity),
              std = sd(activity))
  
  PlotTFActivity(summarized_viper_scores, name, confidence, res, 
                 single_result_path)
  
  SaveTFList(summarized_viper_scores, name, single_result_path)
  
}

#Data from CrossTalkeR paper as test data
#TODO: more generalized usage 
rds_path <- '~/dorothea_execution/20200325_step5.cca.anno.rds'
out_path <- '~/R_results/results_dorothea_ABC'
confidence_level <- 'ABC'

dir.create(out_path)

load(rds_path)
crosstalker <- UpdateSeuratObject(step5.cca)

new.cluster.ids <- c("Neural", "MSC", "Fibroblast", "Megakaryocyte", "Myeloid")
names(new.cluster.ids) <- levels(crosstalker)
crosstalker <- RenameIdents(crosstalker, new.cluster.ids)
crosstalker[["new_annotation"]] <- Idents(object = crosstalker)
crosstalker <- StashIdent(object = crosstalker, save.name = "new_annotation")

png(file = paste0(out_path,'/cluster_per_condition.png'), width = 833, 
    height = 550)
DimPlot(crosstalker, reduction = "umap", label = FALSE, pt.size = 0.5, 
        split.by = "protocol" )
dev.off()

crosstalker.list <- SplitObject(crosstalker, split.by = "protocol")

for (name in names(crosstalker.list)) {
  control <- crosstalker.list[[name]]
  
  name <- str_replace_all(name, "[,;.:-]", "_")
  
  control.averages <- AverageExpression(control)
  write.csv(control.averages[["RNA"]], file = 
              paste0(out_path,'/average_gene_expression_by_cluster_', 
                     name,'.csv'))
  
  DorotheaExecution(control, name, out_path, seq(0.4, 1.5, 0.1), 0.4, 
                    confidence_level)
}

