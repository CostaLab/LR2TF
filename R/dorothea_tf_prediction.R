#' Run complete DoRothEA analysis
#'
#' Description
#'
#' @param seuratobject Input Seurat Object
#' @param out_path Output path to save results
#' @param confidence_level Curation confidence level to filter DoRothEA regulon (default "ABC")
#' @param organism Organism of sample origin
#' @param condition_Ident Meta data field name that contains condition annotation of cells
#' @param celltype_Ident Meta data field name that contains cell type annotation of cells
#' @param comparison_list List of condition comparisons
#' @param pvalue p-value for filtering
#' @param log2fc log2fc value for filtering
#' @export
dorothea_tf_prediction <- function(seuratobject, out_path, confidence_level = c("A", "B", "C"),
                                   organism = "human", condition_Ident, celltype_Ident, comparison_list = NA, pvalue = 0.05, log2fc = 0.0) {

  Idents(object = seuratobject) <- celltype_Ident

  out_path <- paste0(out_path, "/")
  dir.create(out_path)
  seuratobject = dorothea_base_execution(seuratobject, out_path, confidence_level, organism)

  Idents(object = seuratobject) <- condition_Ident
  if(length(comparison_list)> 0 & length(levels(Idents(seuratobject)))<2){
    comparison_list = NA
    print("Only one condition was found in the data, although a list of comparisons was provided. The analyses are performed only for the present condition!")
  }
  Idents(object = seuratobject) <- celltype_Ident

  if (is.na(comparison_list)[[1]]) {
    seuratobject[['doro_annotation']] <- Idents(object = seuratobject)
    result_list <- list()

    Idents(object = seuratobject) <- condition_Ident
    name <- levels(seuratobject)
    name <- str_replace_all(name, "[,;.:-]", "_")

    seuratobject.averages <- AverageExpression(seuratobject, group.by = celltype_annotation, assays = "RNA")
    write.csv(seuratobject.averages[["RNA"]], file =
      paste0(out_path, '/average_gene_expression_by_cluster_',
             name, '.csv'))

    tf_activity_scores = get_significant_tfs_single(seuratobject, name, out_path)
    result_list[[name]] = tf_activity_scores
    result_list[[paste0(name, "_average_expression")]] = seuratobject.averages[["RNA"]]

    saveRDS(result_list, file = paste0(out_path, "/result_list.RDS"))
    return(result_list)
  } else {
    out_path_compared <- paste0(out_path, "/compared")
    dir.create(out_path_compared)
    compared_significant_tfs <- condition_comparison_significant(seuratobject, out_path_compared, celltype_Ident, condition_Ident, comparison_list)

    plot_condition_tf_activities(compared_significant_tfs, out_path_compared)
    plot_condition_tf_activities_compressed(compared_significant_tfs, out_path_compared)

    seuratobject[['doro_annotation']] <- Idents(object = seuratobject)
    seuratobject_list <- SplitObject(seuratobject, split.by = condition_Ident)

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
          tf_condition_significant = filter(tf_condition_significant, FDR < as.double(pvalue))
          tf_condition_significant = filter(tf_condition_significant, logFC > as.double(log2fc) | logFC < (0 - as.double(log2fc)))
          tf_condition_significant = tf_condition_significant[c("tf", "tag", "CellType")]
          tf_condition_significant = tf_condition_significant %>% rename(gene = tf, cluster = CellType)
          compared_tfs = rbind(compared_tfs, tf_condition_significant)
        }
      }

      compared_tfs = compared_tfs[!duplicated(compared_tfs$gene),]

      name <- str_replace_all(name, "[,;.:-]", "_")

      sub_object.averages <- AverageExpression(sub_object, group.by = celltype_Ident, assays = "RNA")
      write.csv(sub_object.averages[["RNA"]], file =
        paste0(out_path, '/average_gene_expression_by_cluster_',
               name, '.csv'))

      tf_activity_scores = get_significant_tfs(sub_object, name, out_path, compared_tfs, pval = pvalue, log2fc = log2fc)
      result_list[[name]] = tf_activity_scores
      result_list[[paste0(name, "_average_expression")]] = sub_object.averages[["RNA"]]
    }
    saveRDS(result_list, file = paste0(out_path, "/result_list.RDS"))
    return(result_list)
  }

}

