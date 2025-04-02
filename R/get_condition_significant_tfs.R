#' Generate cluster and condition heatmap with r effect size only for significant genes
#'
#' Description
#'
#' @param seuratobject Input Seurat Object
#' @param out_path Output path to save results
#' @param comparison_list list of wished comparisons
#' @param num_cell_filter minimum number of cells in each cell type
#' @import glue
#' @import maditr
#' @import Seurat
#' @import scran
#' @import rcompanion
#' @export
condition_comparison_significant <- function(seuratobject, out_path, comparison_list, num_cell_filter = 0, test_type = "binom") {
  comparison_df_list <- list()

  for (comparison in comparison_list) {
    cond1 <- comparison[1]
    cond2 <- comparison[2]

    message("vs: ", cond1, " ", cond2, " ", date(), "\n")
    Idents(seuratobject) <- "tf_condition"

    all_tf_list <- rownames(seuratobject@assays$tf_activities)

    res <- list()
    comparison_sub <- subset(seuratobject, cells = rownames(seuratobject@meta.data)[seuratobject@meta.data$tf_condition %in% comparison])
    if (length(unique(comparison_sub@meta.data$tf_condition)) == 2) {
      for (i in unique(comparison_sub@meta.data$tf_annotation)) {
        a_sub <- subset(comparison_sub, cells = rownames(comparison_sub@meta.data)[comparison_sub@meta.data$tf_annotation == i])
        if (length(unique(a_sub@meta.data$tf_condition)) == 2) {
          condition_table <- as.data.frame(a_sub@meta.data$tf_condition)
          names(condition_table)[1] <- "condition"
          metadata_counts <- condition_table %>%
            group_by(condition) %>%
            summarise(total_count = n())
          if (all(metadata_counts$total_count > num_cell_filter)) {
            g <- as.character(a_sub@meta.data$tf_condition)
            g <- factor(g, levels = c(cond1, cond2))
            res[[i]] <- scran::findMarkers(as.matrix(a_sub@assays$tf_activities@scale.data), g, test.type = test_type)[[1]]
            res[[i]] <- as.data.frame(res[[i]])
            r <- sapply(all_tf_list, function(single_tf) rcompanion::wilcoxonR(as.vector(a_sub@assays$tf_activities@scale.data[single_tf, ]), g))
            nms <- sapply(stringr::str_split(names(r), "\\."), function(x) x[1])
            names(r) <- nms
            res[[i]][nms, "r"] <- r
            res[[i]] <- res[[i]][nms, ]
          }
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

    end_res <- res_df

    comparison_df_list[[glue("{cond1} vs {cond2}")]] <- end_res
    write.csv(res_df, paste0(out_path, "/all_tfs_", glue("{cond1}_vs_{cond2}", ".csv")))
  }

  saveRDS(comparison_df_list, file = paste0(out_path, "/comparison_dfs.RDS"))
  return(comparison_df_list)
}



#' Generate cluster and condition heatmap with r effect size only for significant genes
#'
#' Description
#'
#' @param seuratobject Input Seurat Object
#' @param out_path Output path to save results
#' @param comparison_list list of wished comparisons
#' @param num_cell_filter minimum number of cells in each cell type
#' @import glue
#' @import maditr
#' @import Seurat
#' @import scran
#' @import rcompanion
#' @export
condition_comparison_significant_Seurat <- function(seuratobject, out_path, comparison_list, num_cell_filter = 0, test_type = "binom") {
  comparison_df_list <- list()

  for (comparison in comparison_list) {
    cond1 <- comparison[1]
    cond2 <- comparison[2]

    message("vs: ", cond1, " ", cond2, " ", date(), "\n")
    Idents(seuratobject) <- "tf_condition"

    all_tf_list <- rownames(seuratobject@assays$tf_activities)

    res <- list()
    comparison_sub <- subset(seuratobject, cells = rownames(seuratobject@meta.data)[seuratobject@meta.data$tf_condition %in% comparison])
    if (length(unique(comparison_sub@meta.data$tf_condition)) == 2) {
      seuratobject_list <- SplitObject(comparison_sub, split.by = "tf_annotation")
      results_condition <- list()
      r_values <- list()
      for (name in names(seuratobject_list)) {
        sub_object <- seuratobject_list[[name]]
        if (length(unique(sub_object@meta.data$tf_condition)) == 2) {
          condition_table <- as.data.frame(sub_object@meta.data$tf_condition)
          names(condition_table)[1] <- "condition"
          metadata_counts <- condition_table %>%
            group_by(condition) %>%
            summarise(total_count = n())
          if (all(metadata_counts$total_count > num_cell_filter)) {
            g <- as.character(sub_object@meta.data$tf_condition)
            g <- factor(g, levels = c(cond1, cond2))

            condition_tfs <- FindMarkers(sub_object, ident.1 = cond1, ident.2 = cond2, only.pos = FALSE, min.pct = 0.1)
            condition_tfs$CellType <- name
            condition_tfs$tf <- rownames(condition_tfs)
            colnames(condition_tfs) <- c("p.value", "logFC", "pct.1", "pct.2", "FDR", "CellType", "tf")
            results_condition[[name]] <- condition_tfs

            r <- sapply(all_tf_list, function(single_tf) rcompanion::wilcoxonR(as.vector(sub_object@assays$tf_activities@scale.data[single_tf, ]), g))
            nms <- sapply(stringr::str_split(names(r), "\\."), function(x) x[1])
            names(r) <- nms
            r_df <- as.data.frame(r)
            r_df$CellType <- name
            r_df$tf <- rownames(r_df)
            r_values[[name]] <- r_df
          }
        }
      }
    
      res_df <- do.call("rbind", results_condition)
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

      r_df <- do.call("rbind", r_values)

      end_res <- merge(r_df, res_df, by = 0, all = TRUE)
      end_res <- end_res[, !grepl("\\.y$", colnames(end_res))]
      colnames(end_res) <- gsub("\\.x$", "", colnames(end_res))
      end_res$tag[is.na(end_res$tag)] <- "ns"
      end_res[is.na(end_res)] <- 0.0
      rownames(end_res) <- end_res$Row.names
      end_res$Row.names <- NULL

      comparison_df_list[[glue("{cond1} vs {cond2}")]] <- end_res
      write.csv(res_df, paste0(out_path, "/all_tfs_", glue("{cond1}_vs_{cond2}", ".csv")))
    }
    }
  saveRDS(comparison_df_list, file = paste0(out_path, "/comparison_dfs.RDS"))
  return(comparison_df_list)
}