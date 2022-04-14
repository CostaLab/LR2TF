#' Plots transcription factor activities compared between cluster and conditions
#'
#' Description
#'
#' @param tf_activity_tables Table with tf activities
#' @param out_path Output path to save results
#' @import qpdf
plot_condition_tf_activities <-
  function(tf_activity_tables, out_path) {

    tmp_out_path = paste0(out_path, "/tmp")
    dir.create(tmp_out_path)

    for (nm in names(tf_activity_tables)) {
      nm_df <- data.frame(tf_activity_tables[[nm]])

      plot_width = (((length(unique(nm_df$CellType))) * 15) / 25.4) + 5
      plot_height = (((length(unique(nm_df$tf))) * 4) / 25.4) + 5

      pdf(paste0(tmp_out_path, "/", chartr(" ", "_", nm), ".pdf"), height = plot_height, width = plot_width)

      tag_mapping = nm_df[c("tf", "tag", "CellType")]
      tag_mapping = dcast(tag_mapping, tf ~ CellType, value.var = "tag")
      row.names(tag_mapping) = tag_mapping$tf

      nm_df_short = nm_df[c("r", "tf", "CellType")]
      nm_df_clust = tapply(nm_df_short$r, list(tf = nm_df_short$tf, CellType = nm_df_short$CellType), mean)
      nm_df_clust = data.frame(nm_df_clust)
      nm_df_clust = as.matrix(nm_df_clust)

      fh = function(x) fastcluster::hclust(dist(x))
      p <- Heatmap(nm_df_clust, name = "r", cluster_columns = fh, width = ncol(nm_df_clust) * unit(15, "mm"), row_title = "Transcription Factor", column_title = "Cell Type",
                   height = nrow(nm_df_clust) * unit(4, "mm"), cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(tag_mapping[as.character(rownames(nm_df_clust)[i]), as.character(colnames(nm_df_clust)[j])], x, y, gp = gpar(fontsize = 10))
        })
      draw(p, column_title = glue("{nm}  r effect size"))

      dev.off()
    }

    file_list = list.files(path = tmp_out_path, pattern = "*.pdf", full.names = TRUE)
    qpdf::pdf_combine(input = file_list,
                      output = paste0(out_path, "/cluster_condition_activity_difference.pdf"))

    unlink(tmp_out_path, recursive = TRUE)
  }

#' Plots transcription factor activities compared between cluster and conditions
#'
#' Description
#'
#' @param tf_activity_tables Table with tf activities
#' @param out_path Output path to save results
#' @import qpdf
plot_condition_tf_activities_compressed <-
  function(tf_activity_tables, out_path) {
    tmp_out_path = paste0(out_path, "/tmp")
    dir.create(tmp_out_path)

    for (nm in names(tf_activity_tables)) {
      nm_df <- data.frame(tf_activity_tables[[nm]])

      plot_width = (((length(unique(nm_df$CellType))) * 15) / 25.4) + 5
      plot_height = (((length(unique(nm_df$tf))) * 4) / 25.4) + 5

      pdf(paste0(tmp_out_path, "/", chartr(" ", "_", nm), ".pdf"), height = plot_height, width = plot_width)

      nm_df_short = nm_df[c("r", "tf", "CellType")]
      nm_df_clust = tapply(nm_df_short$r, list(tf = nm_df_short$tf, CellType = nm_df_short$CellType), mean)
      nm_df_clust = data.frame(nm_df_clust)
      #nm_df_clust$tf <- rownames(nm_df_clust)
      nm_df_clust = as.matrix(nm_df_clust)

      fh = function(x) fastcluster::hclust(dist(x))
      p <- Heatmap(as.matrix(nm_df_clust), name = "r", cluster_rows = fh, width = ncol(nm_df_clust) * unit(15, "mm"),
                   height = nrow(nm_df_clust) * unit(0.1, "mm"), show_row_names = FALSE, row_title = "Transcription Factor", column_title = "Cell Type")
      draw(p, column_title = glue("{nm}  r effect size"))
      dev.off()
    }

    file_list = list.files(path = tmp_out_path, pattern = "*.pdf", full.names = TRUE)
    qpdf::pdf_combine(input = file_list,
                      output = paste0(out_path, "/cluster_condition_activity_difference.pdf"))

    unlink(tmp_out_path, recursive = TRUE)
  }


#' Plot heatmap of highly variable transcription factors
#'
#' This function plots a heatmap with the transcription factor activity per
#' cell type. Only the most variable transcription factors are included, which
#' is determined by the statistical variance of the activity scores.
#'
#' @param tf_scores data frame with transcription factor activity scores per cell type
#' @param condition Sample condition for file naming(e.g. control, disease ...)
#' @param out_path Output path to save results
#' @param clusters Number of clusters
#' @import dplyr
#' @import tibble
#' @import pheatmap
#' @import tidyr
#' @import stringr
#' @export
plot_highly_variable_tfs <-
  function(tf_scores, condition, out_path, clusters) {
    sort_factor = clusters * 20

    highly_variable_tfs <- tf_scores %>%
      group_by(tf) %>%
      mutate(var = var(avg))  %>%
      ungroup() %>%
      top_n(sort_factor, var) %>%
      distinct(tf)

    summarized_scores_df <- tf_scores %>%
      semi_join(highly_variable_tfs, by = "tf") %>%
      dplyr::select(-std) %>%
      spread(tf, avg) %>%
      data.frame(row.names = 1, check.names = FALSE)

    summarized_scores_df = as.matrix(t(summarized_scores_df))

    plot_width = ((ncol(summarized_scores_df) * 15) / 25.4) + 5
    plot_height = ((nrow(summarized_scores_df) * 4) / 25.4) + 5

    pdf(
      file = paste0(out_path,
                    '/tf_activity_top20_variable_', condition, '.pdf'),
      width = plot_width,
      height = plot_height
    )

    fh = function(x)
      fastcluster::hclust(dist(x))
    p <-
      Heatmap(
        summarized_scores_df,
        name = "z-score",
        cluster_columns = fh,
        width = ncol(summarized_scores_df) * unit(15, "mm"),
        height = nrow(summarized_scores_df) * unit(4, "mm"),
        row_title = "Transcription Factor",
        column_title = "Cell Type"
      )
    draw(p)

    #print(viper_hmap)
    dev.off()
  }


#' Plot heatmap of all transcription factor activities compressed version
#'
#' This function plots a heatmap with the transcription factor activities per
#' cell type.
#'
#' @param tf_scores data frame with transcription factor activity scores per cell type
#' @param condition Sample condition for file naming(e.g. control, disease ...)
#' @param tag_mapping Labels for heatmap
#' @param out_path Output path to save results
#' @import dplyr
#' @import tibble
#' @import pheatmap
#' @import tidyr
#' @import stringr
#' @export
plot_tf_activity_compressed <-
  function(tf_scores, condition, out_path) {
    plot_width = ((ncol(tf_scores) * 15) / 25.4) + 5
    plot_height = ((nrow(tf_scores) * 0.4) / 25.4) + 5

    tf_scores = as.matrix(tf_scores)
    pdf(
      file = paste0(out_path, '/tf_activity_compressed_', condition, '.pdf'),
      height = plot_height,
      width = plot_width
    )

    fh = function(x)
      fastcluster::hclust(dist(x))
    p <-
      Heatmap(
        tf_scores,
        name = "z-score",
        cluster_rows = fh,
        width = ncol(tf_scores) * unit(15, "mm"),
        height = nrow(tf_scores) * unit(0.4, "mm"),
        row_title = "Transcription Factor",
        column_title = "Cell Type",
        show_row_names = FALSE
      )
    draw(p)
    dev.off()
  }