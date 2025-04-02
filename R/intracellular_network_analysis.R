#' Run IntraTalker analysis
#'
#' Perform intracellular network anaysis based on conditions and clusters using decoupleR TF activity results and scRNA-seq data.
#'
#' @param seuratobject Input Seurat Object
#' @param tf_activities matrix with TF activities for each cell in the scRMA-seq data
#' @param arguments_list named R list with custom options for the analysis
#' @return TF object containing cluster and condition specific TFs as well as CrossTalkeR input
#' @import Seurat
#' @import stringr
#' @import dplyr
#' @export
IntraTalker_analysis <- function(seuratobject, tf_activities = NA, arguments_list) {
  if (is.character(seuratobject)) {
    seuratobject <- readRDS(seuratobject)
  }

  arguments_list <- validate_input_arguments(arguments_list)

  Idents(object = seuratobject) <- arguments_list$celltype
  dir.create(arguments_list$out_path)
  tf_path <- paste0(arguments_list$out_path, "TF_results/")
  dir.create(tf_path)

  if (!is.na(tf_activities)[[1]]) {
    if (is.character(tf_activities)) {
      tf_activities <- t(read.csv(tf_activities, header = TRUE, row.names = 1))
    }
    ## TODO check if the matrix is in the right format
    tf_activities <- CreateAssayObject(data = tf_activities)
    seuratobject[["tf_activities"]] <- tf_activities
  }

  Idents(object = seuratobject) <- arguments_list$condition
  if (length(arguments_list$comparison_list) > 0 & length(levels(Idents(seuratobject))) < 2) {
    arguments_list$comparison_list <- NA
    print("Only one condition was found in the data, although a list of comparisons was provided. The analyses are performed only for the present condition!")
  }

  DefaultAssay(object = seuratobject) <- "tf_activities"
  seuratobject <- ScaleData(seuratobject)

  Idents(object = seuratobject) <- arguments_list$condition
  seuratobject[["tf_condition"]] <- Idents(object = seuratobject)

  Idents(object = seuratobject) <- arguments_list$celltype
  seuratobject[["tf_annotation"]] <- Idents(object = seuratobject)
  #annotation_df <- data.frame(cell_type = seuratobject@meta.data[[arguments_list$celltype]], row.names = rownames(seuratobject@meta.data))
  #annotation_df$cell_type <- gsub("_", "-", annotation_df$cell_type)
  #seuratobject@meta.data$tf_annotation <- annotation_df$cell_type

  if (is.na(arguments_list$comparison_list)[[1]]) {
    result_list <- list()
    gene_expression_list <- list()
    CTR_cluster_list <- list()
    intranet_cluster_list <- list()

    Idents(object = seuratobject) <- "tf_condition"
    seuratobject_list <- SplitObject(seuratobject,
      split.by = "tf_condition"
    )
    for (name in names(seuratobject_list)) {
      sub_object <- seuratobject_list[[name]]

      name <- str_replace_all(name, "[,;.:-]", "_")

      sub_object.averages <- AverageExpression(sub_object,
        group.by = "tf_annotation",
        assays = "RNA"
      )
      # write.csv(seuratobject.averages[["RNA"]], file =
      #   paste0(arguments_list$out_path, 'average_gene_expression_by_cluster_',
      #          name, '.csv'))

      #colnames(sub_object.averages[["RNA"]]) <- unique(sub_object[[arguments_list$celltype]][[1]])

      tf_activity_scores <- get_significant_tfs(sub_object,
        name,
        tf_path,
        pval = arguments_list$pval,
        log2fc = arguments_list$logfc
      )
      result_list[[name]] <- tf_activity_scores
      gene_expression_list[[paste0(name, "_average_expression")]] <- sub_object.averages[["RNA"]]

      CTR_cluster_list[[name]] <- generate_CrossTalkeR_input(
        tf_activity_scores[["cluster"]],
        gene_expression_list[[paste0(name, "_average_expression")]],
        arguments_list$reg,
        arguments_list$organism
      )

      intranet_cluster_list[[name]] <- generate_intracellular_network(
        tf_activity_scores[["cluster"]],
        gene_expression_list[[paste0(name, "_average_expression")]],
        arguments_list$reg,
        arguments_list$organism
      )
    }
    tf <- new("TFObj",
      tf_activities_condition = list(),
      tf_activities_cluster = result_list,
      average_gene_expression = gene_expression_list,
      regulon = arguments_list$reg,
      CTR_input_condition = list(),
      CTR_input_cluster = CTR_cluster_list,
      intracellular_network_condition = list(),
      intracellular_network_cluster = intranet_cluster_list
    )

    saveRDS(tf, file = paste0(tf_path, "result_TF_object.RDS"))
    # saveRDS(seuratobject, file = paste0(tf_path, "TF_seurat_object.RDS"))
    return(tf)
  } else {
    out_path_compared <- paste0(tf_path, "compared")
    dir.create(out_path_compared)

    if(arguments_list$Seurat) {
      compared_significant_tfs <- condition_comparison_significant_Seurat(
      seuratobject,
      out_path_compared,
      arguments_list$comparison_list,
      arguments_list$num_cell_filter,
      arguments_list$test_type)
    } else {
      compared_significant_tfs <- condition_comparison_significant(
      seuratobject,
      out_path_compared,
      arguments_list$comparison_list,
      arguments_list$num_cell_filter,
      arguments_list$test_type)
    }
    
    plot_condition_tf_activities(compared_significant_tfs, out_path_compared)
    plot_condition_tf_activities_compressed(compared_significant_tfs, out_path_compared)

    seuratobject_list <- SplitObject(seuratobject, split.by = arguments_list$condition)

    result_condition_list <- list()
    result_cluster_list <- list()
    gene_expression_list <- list()
    CTR_condition_list <- list()
    CTR_cluster_list <- list()
    intranet_condition_list <- list()
    intranet_cluster_list <- list()
    for (name in names(seuratobject_list)) {
      if (name %in% unique(unlist(arguments_list$comparison_list))) {
        sub_object <- seuratobject_list[[name]]

        compared_tfs <- data.frame(
          gene = character(),
          tag = character(),
          cluster = character()
        )

        for (result_name in names(compared_significant_tfs)) {
          if (grepl(name, result_name, fixed = TRUE)) {
            tf_condition_significant <- compared_significant_tfs[[result_name]]
            tf_condition_significant <- filter(
              tf_condition_significant,
              FDR < as.double(arguments_list$pval)
            )
            if(arguments_list$test_type == "wilcox") {
              tf_condition_significant <- filter(
                tf_condition_significant,
                logFC > (0.5 + as.double(arguments_list$logfc)) | logFC < (0.5 - as.double(arguments_list$logfc))
              )
            } else {
              tf_condition_significant <- filter(
              tf_condition_significant,
              logFC > as.double(arguments_list$logfc) | logFC < (0 - as.double(arguments_list$logfc))
            )
            }
            tf_condition_significant <- tf_condition_significant[c("tf", "tag", "CellType")]
            tf_condition_significant <- tf_condition_significant %>%
              rename(
                gene = tf,
                cluster = CellType
              )
            compared_tfs <- rbind(compared_tfs, tf_condition_significant)
          }
        }

        compared_tfs <- compared_tfs[!duplicated(rownames(compared_tfs)), ]

        name <- str_replace_all(name, "[,;.:-]", "_")

        sub_object.averages <- AverageExpression(sub_object,
          group.by = "tf_annotation",
          assays = "RNA"
        )
        # write.csv(sub_object.averages[["RNA"]], file =
        #   paste0(arguments_list$out_path, 'average_gene_expression_by_cluster_',
        #          name, '.csv'))
        #colnames(sub_object.averages[["RNA"]]) <- unique(sub_object[[arguments_list$celltype]][[1]])
        tf_activity_scores <- get_significant_tfs(sub_object,
          name,
          tf_path,
          pval = arguments_list$pval,
          log2fc = arguments_list$logfc,
          compared_tfs,
          TRUE
        )
        result_condition_list[[name]] <- tf_activity_scores[["condition"]]
        result_cluster_list[[name]] <- tf_activity_scores[["cluster"]]
        gene_expression_list[[paste0(name, "_average_expression")]] <- sub_object.averages[["RNA"]]


        CTR_condition_list[[name]] <- generate_CrossTalkeR_input(
          tf_activity_scores[["condition"]],
          gene_expression_list[[paste0(name, "_average_expression")]],
          arguments_list$reg,
          arguments_list$organism
        )
        CTR_cluster_list[[name]] <- generate_CrossTalkeR_input(
          tf_activity_scores[["cluster"]],
          gene_expression_list[[paste0(name, "_average_expression")]],
          arguments_list$reg,
          arguments_list$organism
        )

        intranet_condition_list[[name]] <- generate_intracellular_network(
          tf_activity_scores[["condition"]],
          gene_expression_list[[paste0(name, "_average_expression")]],
          arguments_list$reg,
          arguments_list$organism
        )
        intranet_cluster_list[[name]] <- generate_intracellular_network(
          tf_activity_scores[["cluster"]],
          gene_expression_list[[paste0(name, "_average_expression")]],
          arguments_list$reg,
          arguments_list$organism
        )
      }
    }

    tf <- new("TFObj",
      tf_activities_condition = result_condition_list,
      tf_activities_cluster = result_cluster_list,
      average_gene_expression = gene_expression_list,
      regulon = arguments_list$reg,
      CTR_input_condition = CTR_condition_list,
      CTR_input_cluster = CTR_cluster_list,
      intracellular_network_condition = intranet_condition_list,
      intracellular_network_cluster = intranet_cluster_list
    )

    saveRDS(tf, file = paste0(tf_path, "result_TF_object.RDS"))
    # saveRDS(seuratobject, file = paste0(tf_path, "TF_seurat_object.RDS"))

    return(tf)
  }
}
