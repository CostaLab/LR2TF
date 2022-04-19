
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
condition_comparison_significant <- function(seuratobject, out_path, celltype_annotation, condition_annotation, comparison_list){

  DefaultAssay(object = seuratobject) <- "dorothea"
  seuratobject <- ScaleData(seuratobject)

  Idents(object = seuratobject) <- condition_annotation
  seuratobject[['doro_condition']] <- Idents(object=seuratobject)
  Idents(object = seuratobject) <- celltype_annotation
  seuratobject[['doro_annotation']] <- Idents(object=seuratobject)

  vs_df_list <- list()

  for(vs in comparison_list){
    vs1 <- vs[1]
    vs2 <- vs[2]

    message("vs: ",vs1, " ", vs2, " ", date(), "\n")
    Idents(seuratobject) <- condition_annotation

    pws <- rownames(seuratobject@assays$dorothea)

    ###
    res <- list()
    for(i in levels(seuratobject@meta.data$doro_annotation)){
      a_sub <- subset(seuratobject, cells=rownames(seuratobject@meta.data)[seuratobject@meta.data$doro_annotation==i & (seuratobject@meta.data$doro_condition %in% vs)])
      g <- as.character(a_sub@meta.data$doro_condition)
      g <- factor(g, levels=c(vs1, vs2)) ###############################################
      res[[i]] <- scran::findMarkers(as.matrix(a_sub@assays$dorothea@scale.data), g)[[1]]
      res[[i]] <- as.data.frame(res[[i]])
      r <- sapply(pws, function(pw) rcompanion::wilcoxonR(as.vector(a_sub@assays$dorothea@scale.data[pw,]), g))
      nms <- sapply(stringr::str_split(names(r), "\\."), function(x)x[1])
      names(r) <- nms
      res[[i]][nms, "r"] <- r
      res[[i]] <- res[[i]][nms, ]

    }

    for (cl in names(res)) {
      res[[cl]]$tf <- rownames(res[[cl]])
      res[[cl]]$CellType <- cl
      colnames(res[[cl]]) <-  c("Top","p.value","FDR", "summary.logFC","logFC","r","tf","CellType")

    }
    res_df <- do.call("rbind", res)
    res_df$tag <- sapply(res_df$FDR, function(pval) {
      if(pval< 0.001) {
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

    significant_res <- res_df[res_df$tag == '***',]
    significant_genes <- unique(significant_res$tf)

    end_res <- filter(res_df, res_df$tf %in% significant_genes)

    vs_df_list[[glue("{vs1} vs {vs2}")]] <- end_res
    write.csv(res_df,paste0(out_path,"/all_tfs_",glue("{vs1}_vs_{vs2}", ".csv")))
  }

  return(vs_df_list)
}