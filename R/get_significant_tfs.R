#' Run DoRothEA analysis to get significant transcription factors
#'
#' Description
#'
#' @param seuratobject Input Seurat Object
#' @param condition Experminet condition (e.g. disease, knockout ...)
#' @param out_path Output path to save results
#' @param tf_condition_significant condition comparison results
#' @param pval p-value to filter results
#' @param log2fc log fold change value to filter results
#' @return A data frame with transcription factor activity scores per cell type
#' @export
get_significant_tfs <- function(seuratobject, condition, out_path, tf_condition_significant, pval, log2fc) {
  single_result_path <- paste0(out_path, condition)
  dir.create(single_result_path)

  DefaultAssay(object = seuratobject) <- "tf_activities"
  seuratobject <- ScaleData(seuratobject)

  Idents(object = seuratobject) <- "tf_annotation"

  number_of_clusters = length(levels(Idents(seuratobject)))

  seuratobject.markers <- FindAllMarkers(seuratobject, only.pos = TRUE,
                                         min.pct = 0, logfc.threshold = 0,
                                         verbose = FALSE)
  write.csv(seuratobject.markers, file =
    paste0(single_result_path, '/all_specificmarker_',
           '_', condition, '.csv'))

  seuratobject.markers$tag <- sapply(seuratobject.markers$p_val_adj, function(pval) {
    if (pval < 0.001) {
      txt <- "***"
    } else if (pval < 0.01) {
      txt <- "**"
    } else if (pval < 0.05) {
      txt <- "*"
    } else {
      txt <- "ns"
    }
    return(txt)
  })

  seuratobject.markers$log_fc_tag <- sapply(seuratobject.markers$avg_log2FC, function(log_fc) {
    if (log_fc >= 1.0) {
      txt <- "***"
    } else if (log_fc > 0.5) {
      txt <- "**"
    } else if (log_fc > 0.0) {
      txt <- "*"
    } else {
      txt <- "ns"
    }
    return(txt)
  })

  tag_mapping = seuratobject.markers[c("gene", "tag", "log_fc_tag", "cluster", "avg_log2FC", "p_val_adj")]
  tag_mapping = filter(tag_mapping, p_val_adj < as.double(pval))
  tag_mapping = filter(tag_mapping, avg_log2FC > as.double(log2fc) | avg_log2FC < (0 - as.double(log2fc)))
  #tag_mapping = tag_mapping[(tag_mapping$tag == "***"),]
  #tag_mapping = tag_mapping[(tag_mapping$log_fc_tag == "***"),]
  tag_mapping = dcast(tag_mapping, gene ~ cluster, value.var = "tag")
  row.names(tag_mapping) = tag_mapping$gene
  tag_mapping$gene = NULL
  tag_mapping$log_fc_tag = NULL
  tag_mapping[is.na(tag_mapping)] <- "ns"

  viper_scores_df <- GetAssayData(seuratobject, slot = "scale.data",
                                  assay = "tf_activities") %>%
    data.frame(check.names = F) %>%
    t()

  CellsClusters <- data.frame(cell = names(Idents(seuratobject)),
                              cell_type = as.character(Idents(seuratobject)),
                              check.names = F)

  save_unfiltered_tf_scores(viper_scores_df, CellsClusters, condition, single_result_path)

  unfiltered_viper_scores = unfiltered_tf_activity_table(viper_scores_df, CellsClusters)

  col.num <- which(colnames(viper_scores_df) %in% rownames(tag_mapping))
  viper_scores_df <- viper_scores_df[, sort(c(col.num))]

  viper_scores_clusters <- viper_scores_df %>%
    data.frame() %>%
    rownames_to_column("cell") %>%
    gather(tf, activity, -cell) %>%
    inner_join(CellsClusters)

  summarized_viper_scores <- viper_scores_clusters %>%
    group_by(tf, cell_type) %>%
    summarise(avg = mean(activity),
              std = sd(activity))

  summarized_viper_scores_df <- summarized_viper_scores %>%
    dplyr::select(-std) %>%
    spread(tf, avg) %>%
    data.frame(row.names = 1, check.names = FALSE)

  summarized_viper_scores_df = t(summarized_viper_scores_df)
  write.csv(summarized_viper_scores_df, file = paste0(single_result_path, '/tf_scores', '_', condition, '.csv'))

  tf_scores_variable_table = save_variable_tf_scores(summarized_viper_scores, condition, single_result_path)

  message("Plotting top variables tf activities")
  plot_highly_variable_tfs(summarized_viper_scores, condition,
                           single_result_path, number_of_clusters)

  plot_tf_activity_compressed(summarized_viper_scores_df, condition, single_result_path)
  plot_tf_activity(summarized_viper_scores_df, tag_mapping, condition,
                   single_result_path)

  rownames(summarized_viper_scores_df) = gsub(".", "-", rownames(summarized_viper_scores_df), fixed = TRUE)
  #seuratobject.markers = seuratobject.markers[(seuratobject.markers$tag == "***"),]
  #seuratobject.markers = seuratobject.markers[(seuratobject.markers$log_fc_tag == "***"),]

  map_z_value_filtered <- function(gene, cluster) {
    if (gene %in% rownames(summarized_viper_scores_df)) {
      z_score = summarized_viper_scores_df[as.character(gene), as.character(cluster)]
      return(z_score)
    } else {
      return(NA)
    }
  }

  map_z_value <- function(gene, cluster) {
    if (gene %in% rownames(unfiltered_viper_scores)) {
      z_score = unfiltered_viper_scores[as.character(gene), as.character(cluster)]
      return(z_score)
    } else {
      return(NA)
    }
  }

  res <- list()
  seuratobject.markers$z_score = mapply(map_z_value_filtered, seuratobject.markers$gene, seuratobject.markers$cluster)
  seuratobject.markers = seuratobject.markers[!(seuratobject.markers$tag == "ns"),]
  seuratobject.markers = na.omit(seuratobject.markers)
  res[["cluster"]] = seuratobject.markers[c("gene", "tag", "cluster", "z_score")]
  write.csv(res[["cluster"]], file = paste0(single_result_path, '/significant_cluster_tf_results', '_', condition, '.csv'))

  tf_condition_significant$gene = gsub(".", "-", tf_condition_significant$gene, fixed = TRUE)
  tf_condition_significant$z_score = mapply(map_z_value, tf_condition_significant$gene, tf_condition_significant$cluster)
  tf_condition_significant = na.omit(tf_condition_significant)
  res[["condition"]] = tf_condition_significant
  write.csv(res[["condition"]], file = paste0(single_result_path, '/significant_condition_tf_results', '_', condition, '.csv'))

  return(res)
}

#' Run DoRothEA analysis to get significant transcription factors
#'
#' Description
#'
#' @param seuratobject Input Seurat Object
#' @param condition Experminet condition (e.g. disease, knockout ...)
#' @param out_path Output path to save results
#' @param tf_condition_significant condition comparison results
#' @return A data frame with transcription factor activity scores per cell type
#' @export
get_significant_tfs_single <- function(seuratobject, condition, out_path, pval, log2fc) {
  single_result_path <- paste0(out_path, condition)
  dir.create(single_result_path)

  DefaultAssay(object = seuratobject) <- "tf_activities"
  seuratobject <- ScaleData(seuratobject)

  Idents(object = seuratobject) <- "tf_annotation"

  number_of_clusters = length(levels(Idents(seuratobject)))

  seuratobject.markers <- FindAllMarkers(seuratobject, only.pos = TRUE,
                                         min.pct = 0, logfc.threshold = 0,
                                         verbose = FALSE)
  write.csv(seuratobject.markers, file =
    paste0(single_result_path, '/all_specificmarker_',
           '_', condition, '.csv'))

  seuratobject.markers$tag <- sapply(seuratobject.markers$p_val_adj, function(pval) {
    if (pval < 0.001) {
      txt <- "***"
    } else if (pval < 0.01) {
      txt <- "**"
    } else if (pval < 0.05) {
      txt <- "*"
    } else {
      txt <- "ns"
    }
    return(txt)
  })

  seuratobject.markers$log_fc_tag <- sapply(seuratobject.markers$avg_log2FC, function(log_fc) {
    if (log_fc >= 1.0) {
      txt <- "***"
    } else if (log_fc > 0.5) {
      txt <- "**"
    } else if (log_fc > 0.0) {
      txt <- "*"
    } else {
      txt <- "ns"
    }
    return(txt)
  })

  tag_mapping = seuratobject.markers[c("gene", "tag", "log_fc_tag", "cluster")]
  tag_mapping = filter(tag_mapping, p_val_adj < as.double(pval))
  tag_mapping = filter(tag_mapping, avg_log2FC > as.double(log2fc) | avg_log2FC < (0 - as.double(log2fc)))
  #tag_mapping = tag_mapping[(tag_mapping$tag == "***"),]
  #tag_mapping = tag_mapping[(tag_mapping$log_fc_tag == "***"),]
  tag_mapping = dcast(tag_mapping, gene ~ cluster, value.var = "tag")
  row.names(tag_mapping) = tag_mapping$gene
  tag_mapping$gene = NULL
  tag_mapping$log_fc_tag = NULL
  tag_mapping[is.na(tag_mapping)] <- "ns"

  viper_scores_df <- GetAssayData(seuratobject, slot = "scale.data",
                                  assay = "tf_activities") %>%
    data.frame(check.names = F) %>%
    t()

  CellsClusters <- data.frame(cell = names(Idents(seuratobject)),
                              cell_type = as.character(Idents(seuratobject)),
                              check.names = F)

  save_unfiltered_tf_scores(viper_scores_df, CellsClusters, condition, single_result_path)

  col.num <- which(colnames(viper_scores_df) %in% rownames(tag_mapping))
  viper_scores_df <- viper_scores_df[, sort(c(col.num))]

  viper_scores_clusters <- viper_scores_df %>%
    data.frame() %>%
    rownames_to_column("cell") %>%
    gather(tf, activity, -cell) %>%
    inner_join(CellsClusters)

  summarized_viper_scores <- viper_scores_clusters %>%
    group_by(tf, cell_type) %>%
    summarise(avg = mean(activity),
              std = sd(activity))

  summarized_viper_scores_df <- summarized_viper_scores %>%
    dplyr::select(-std) %>%
    spread(tf, avg) %>%
    data.frame(row.names = 1, check.names = FALSE)

  summarized_viper_scores_df = t(summarized_viper_scores_df)
  write.csv(summarized_viper_scores_df, file = paste0(single_result_path, '/tf_scores', '_', condition, '.csv'))

  tf_scores_variable_table = save_variable_tf_scores(summarized_viper_scores, condition, single_result_path)

  message("Plotting top variables tf activities")
  plot_highly_variable_tfs(summarized_viper_scores, condition,
                           single_result_path, number_of_clusters)

  plot_tf_activity_compressed(summarized_viper_scores_df, condition, single_result_path)
  plot_tf_activity(summarized_viper_scores_df, tag_mapping, condition,
                   single_result_path)

  rownames(summarized_viper_scores_df) = gsub(".", "-", rownames(summarized_viper_scores_df), fixed = TRUE)
  #seuratobject.markers = seuratobject.markers[(seuratobject.markers$tag == "***"),]
  #seuratobject.markers = seuratobject.markers[(seuratobject.markers$log_fc_tag == "***"),]

  map_z_value_filtered <- function(gene, cluster) {
    if (gene %in% rownames(summarized_viper_scores_df)) {
      z_score = summarized_viper_scores_df[as.character(gene), as.character(cluster)]
      return(z_score)
    } else {
      return(NA)
    }
  }

  res <- list()
  seuratobject.markers$z_score = mapply(map_z_value_filtered, seuratobject.markers$gene, seuratobject.markers$cluster)
  seuratobject.markers = seuratobject.markers[!(seuratobject.markers$tag == "ns"),]
  seuratobject.markers = na.omit(seuratobject.markers)
  res[["cluster"]] = seuratobject.markers[c("gene", "tag", "cluster", "z_score")]
  write.csv(res[["cluster"]], file = paste0(single_result_path, '/significant_cluster_tf_results', '_', condition, '.csv'))

  return(res)
}