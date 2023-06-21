#' Run TF activity analysis
#'
#' Description
#'
#' @param seuratobject Input Seurat Object
#' @param tf_activities matrix with TF activities for each cell in the scRMA-seq data
#' @param arguments_list named R list with custom options for the analysis
#' @export
dorothea_tf_prediction <- function(seuratobject, tf_activities = NA, arguments_list) {

  if (typeof(seuratobject) == "character") {
    seuratobject <- readRDS(seuratobject)
  }

  Idents(object = seuratobject) <- arguments_list$celltype
  dir.create(arguments_list$out_path)

  if (is.na(tf_activities)[[1]]) {
    seuratobject <- decoupleR_viper_analysis(seuratobject,)
  } else {
    if (typeof(tf_activities) == "character") {
      tf_activities <- t(read.csv(tf_activities, header = TRUE, row.names = 1))
    }

    tf_acticities <- CreateAssayObject(data = tf_activities)
    seuratobject[["tf_acticities"]] <- tf_acticities
  }

  Idents(object = seuratobject) <- arguments_list$condition
  if (length(arguments_list$comparison_list) > 0 & length(levels(Idents(seuratobject))) < 2) {
    arguments_list$comparison_list <- NA
    print("Only one condition was found in the data, although a list of comparisons was provided. The analyses are performed only for the present condition!")
  }
  Idents(object = seuratobject) <- arguments_list$celltype

  if (is.na(arguments_list$comparison_list)[[1]]) {
    seuratobject[['tf_annotation']] <- Idents(object = seuratobject)
    result_list <- list()

    Idents(object = seuratobject) <- arguments_list$condition
    name <- levels(seuratobject)
    name <- str_replace_all(name, "[,;.:-]", "_")

    seuratobject.averages <- AverageExpression(seuratobject, group.by = arguments_list$celltype, assays = "RNA")
    write.csv(seuratobject.averages[["RNA"]], file =
      paste0(arguments_list$out_path, 'average_gene_expression_by_cluster_',
             name, '.csv'))

    tf_activity_scores = get_significant_tfs_single(seuratobject, name, arguments_list$out_path)
    result_list[[name]] = tf_activity_scores
    result_list[[paste0(name, "_average_expression")]] = seuratobject.averages[["RNA"]]

    saveRDS(result_list, file = paste0(arguments_list$out_path, "result_list.RDS"))
    saveRDS(seuratobject, file = paste0(out_arguments_list$out_pathpath, "result_seurat_object.RDS"))
    return(result_list)
  } else {
    out_path_compared <- paste0(arguments_list$out_path, "compared")
    dir.create(out_path_compared)
    compared_significant_tfs <- condition_comparison_significant(seuratobject, out_path_compared, arguments_list$celltype, arguments_list$condition, arguments_list$comparison_list)

    plot_condition_tf_activities(compared_significant_tfs, out_path_compared)
    plot_condition_tf_activities_compressed(compared_significant_tfs, out_path_compared)

    seuratobject[['tf_annotation']] <- Idents(object = seuratobject)
    seuratobject_list <- SplitObject(seuratobject, split.by = arguments_list$condition)

    result_list <- list()
    for (name in names(seuratobject_list)) {
      sub_object <- seuratobject_list[[name]]

      compared_tfs = data.frame(
        gene = character(),
        tag = character(),
        cluster = character()
      )

      for (result_name in names(compared_significant_tfs)) {
        if (grepl(name, result_name, fixed = TRUE)) {
          tf_condition_significant = compared_significant_tfs[[result_name]]
          #tf_condition_significant = tf_condition_significant[(tf_condition_significant$tag == "***"),]
          tf_condition_significant = filter(tf_condition_significant, FDR < as.double(arguments_list$pval))
          tf_condition_significant = filter(tf_condition_significant, logFC > as.double(arguments_list$logfc) | logFC < (0 - as.double(arguments_list$logfc)))
          tf_condition_significant = tf_condition_significant[c("tf", "tag", "CellType")]
          tf_condition_significant = tf_condition_significant %>% rename(gene = tf, cluster = CellType)
          compared_tfs = rbind(compared_tfs, tf_condition_significant)
        }
      }

      compared_tfs = compared_tfs[!duplicated(rownames(compared_tfs)),]

      name <- str_replace_all(name, "[,;.:-]", "_")

      sub_object.averages <- AverageExpression(sub_object, group.by = arguments_list$celltype, assays = "RNA")
      write.csv(sub_object.averages[["RNA"]], file =
        paste0(arguments_list$out_path, 'average_gene_expression_by_cluster_',
               name, '.csv'))

      tf_activity_scores = get_significant_tfs(sub_object, name, arguments_list$out_path, compared_tfs, pval = arguments_list$pval, log2fc = arguments_list$logfc)
      result_list[[name]] = tf_activity_scores
      result_list[[paste0(name, "_average_expression")]] = sub_object.averages[["RNA"]]
    }
    saveRDS(result_list, file = paste0(arguments_list$out_path, "result_list.RDS"))
    saveRDS(seuratobject, file = paste0(arguments_list$out_path, "result_seurat_object.RDS"))
    return(result_list)
  }

}

