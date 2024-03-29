#' Generate cluster and condition heatmap with r effect size only for significant genes
#'
#' Description
#'
#' @param seuratobject Input Seurat Object
#' @param out_path Output path to save results
#' @param celltype_annotation meta data field with celltype annotations
#' @param condition_annotation meta data field with condition annotation
#' @param comparison_list list of wished comparisons
#' @import glue
#' @import maditr
#' @export
condition_comparison_significant <- function(seuratobject, out_path, celltype_annotation, condition_annotation, comparison_list) {

  DefaultAssay(object = seuratobject) <- "tf_activities"
  seuratobject <- ScaleData(seuratobject)

  Idents(object = seuratobject) <- condition_annotation
  seuratobject[['tf_condition']] <- Idents(object = seuratobject)
  Idents(object = seuratobject) <- celltype_annotation
  seuratobject[['tf_annotation']] <- Idents(object = seuratobject)

  vs_df_list <- list()

  for (vs in comparison_list) {
    vs1 <- vs[1]
    vs2 <- vs[2]

    message("vs: ", vs1, " ", vs2, " ", date(), "\n")
    Idents(seuratobject) <- condition_annotation

    pws <- rownames(seuratobject@assays$tf_activities)

    ###
    res <- list()
    for (i in levels(seuratobject@meta.data$tf_annotation)) {
      a_sub <- subset(seuratobject, cells = rownames(seuratobject@meta.data)[seuratobject@meta.data$tf_annotation == i & (seuratobject@meta.data$tf_condition %in% vs)])
      if (length(unique(a_sub@meta.data$tf_condition)) == 2) {
        condition_table <- as.data.frame(a_sub@meta.data$tf_condition)
        names(condition_table)[1] <- "condition"
        metadata_counts <- condition_table %>%
          group_by(condition) %>%
          summarise(total_count = n())
        if (all(metadata_counts$total_count > 10)) {
          g <- as.character(a_sub@meta.data$tf_condition)
          g <- factor(g, levels = c(vs1, vs2)) ###############################################
          res[[i]] <- scran::findMarkers(as.matrix(a_sub@assays$tf_activities@scale.data), g)[[1]]
          res[[i]] <- as.data.frame(res[[i]])
          r <- sapply(pws, function(pw) rcompanion::wilcoxonR(as.vector(a_sub@assays$tf_activities@scale.data[pw,]), g))
          nms <- sapply(stringr::str_split(names(r), "\\."), function(x)x[1])
          names(r) <- nms
          res[[i]][nms, "r"] <- r
          res[[i]] <- res[[i]][nms,]
        }
      }
    }

    for (cl in names(res)) {
      res[[cl]]$tf <- rownames(res[[cl]])
      res[[cl]]$CellType <- cl
      colnames(res[[cl]]) <- c("Top", "p.value", "FDR", "summary.logFC", "logFC", "r", "tf", "CellType")

    }
    res_df <- do.call("rbind", res)
    res_df <- na.omit(res_df)
    res_df$tag <- sapply(res_df$FDR, function(pval) {
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

    #significant_res <- res_df[res_df$tag == '***',]
    #significant_genes <- unique(significant_res$tf)

    #end_res <- filter(res_df, res_df$tf %in% significant_genes)
    end_res <- res_df

    vs_df_list[[glue("{vs1} vs {vs2}")]] <- end_res
    write.csv(res_df, paste0(out_path, "/all_tfs_", glue("{vs1}_vs_{vs2}", ".csv")))
  }

  saveRDS(vs_df_list, file = paste0(out_path, "/comparison_dfs.RDS"))
  return(vs_df_list)
}