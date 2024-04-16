#' Generate CrossTalkeR input from significant tf table
#'
#' This function loads the transcription factor activity table/data frame
#' for multiple cell types, and generates the CrossTalkeR input table.
#' The returned input table contains receptor-transcription factor and
#' transcription factor-ligand interactions based on the OmniPath database and
#' the DoRothEA regulon.
#'
#' @param tf_activities Data Frame with transcription factor activities by cell type
#' @param confidence_level Confidence Level for the DoRothEA regulon used to get the transcription factor activities
#' @param gene_activities Table with average gene expression levels
#' @return A data frame with CrossTalkeR input
#' @import dplyr
#' @import tibble
#' @export
generate_CrossTalkeR_input <-
  function(tf_activities,
           gene_expression,
           regulon = NA) {
    if (any(tf_activities$z_score > 0)) {
      regulon <- regulon %>%
        rename(tf = source)

      ligands <- ligands_human

      R2TF <- aggregate(RTF_DB_2$receptor ~ RTF_DB_2$tf, FUN = c)
      colnames(R2TF) <- c('tf', 'receptors')
      R2TF <- R2TF %>%
        remove_rownames %>%
        tibble::column_to_rownames(var = 'tf')

      sorted_regulon <-
        aggregate(regulon$target ~ regulon$tf, FUN = c)
      colnames(sorted_regulon) <- c('tf', 'targets')
      sorted_regulon <- sorted_regulon %>%
        remove_rownames %>%
        tibble::column_to_rownames(var = 'tf')

      tf_activities <- tf_activities %>%
        filter(z_score > 0)

      output_df <- create_empty_CTR_dataframe()

      for (row in 1:nrow(tf_activities)) {

        r_tf <- create_empty_CTR_dataframe()
        tf_l <- create_empty_CTR_dataframe()

        #if (tf_activities[row, "z_score"] > 0) {
        tf <- as.character(tf_activities[row, "gene"])
        targets <- sorted_regulon[tf,][1]
        receptors <- R2TF[tf,][1]
        tf_ligands <- intersect(targets[[1]], ligands)
        if (length(tf_ligands) > 0) {
          for (ligand in tf_ligands) {
            expressed <- FALSE
            if (ligand %in% rownames(gene_expression)) {
              ex_value <- gene_expression[ligand, tf_activities[row, "cluster"]]
              if (ex_value != 0) {
                expressed <- TRUE
              }
            }

            if (expressed == TRUE) {
              df <- add_entry_to_CTR_dataframe(tf_activities[row, "cluster"],
                                               tf_activities[row, "cluster"],
                                               tf_activities[row, "gene"],
                                               ligand,
                                               'Transcription Factor',
                                               'Ligand',
                                               tf_activities[row, "z_score"])
              tf_l <- rbind(tf_l, df)
            }
          }
        }
        if (length(receptors[[1]]) > 0) {
          for (receptor in receptors[[1]]) {
            df <- add_entry_to_CTR_dataframe(tf_activities[row, "cluster"],
                                             tf_activities[row, "cluster"],
                                             receptor,
                                             tf_activities[row, "gene"],
                                             'Receptor',
                                             'Transcription Factor',
                                             tf_activities[row, "z_score"])
            r_tf <- rbind(r_tf, df)
          }
        }
        #}

        r_tf$gene_A <- gsub("_", "+", r_tf$gene_A, fixed = TRUE)
        r_tf$gene_B <- gsub("_", "+", r_tf$gene_B, fixed = TRUE)
        tf_l$gene_A <- gsub("_", "+", tf_l$gene_A, fixed = TRUE)
        tf_l$gene_B <- gsub("_", "+", tf_l$gene_B, fixed = TRUE)

        output_df <- rbind(output_df, r_tf)
        output_df <- rbind(output_df, tf_l)
      }
    } else {
      output_df = NA
    }

    return(output_df)
  }


#' Generate CrossTalkeR input from significant tf table for mouse data
#'
#' This function loads the transcription factor activity table/data frame
#' for multiple cell types, and generates the CrossTalkeR input table.
#' The returned input table contains receptor-transcription factor and
#' transcription factor-ligand interactions based on the OmniPath database and
#' the DoRothEA regulon.
#'
#' @param tf_activities Data Frame with transcription factor activities by cell type
#' @param confidence_level Confidence Level for the DoRothEA regulon used to get the transcription factor activities
#' @param gene_activities Table with average gene expression levels
#' @return A data frame with CrossTalkeR input
#' @import dplyr
#' @import tibble
#' @export
generate_CrossTalkeR_input_mouse <-
  function(tf_activities,
           gene_expression,
           regulon = NA) {
    if (any(tf_activities$z_score > 0)) {
      ligands <- converted_ligands

      R2TF <- aggregate(RTF_DB_mouse$receptor ~ RTF_DB_mouse$tf, FUN = c)
      colnames(R2TF) <- c('tf', 'receptors')
      R2TF <- R2TF %>%
        remove_rownames %>%
        tibble::column_to_rownames(var = 'tf')

      regulon <- regulon %>%
        rename(tf = source)
      sorted_regulon <-
        aggregate(regulon$target ~ regulon$tf, FUN = c)
      colnames(sorted_regulon) <- c('tf', 'targets')
      sorted_regulon <- sorted_regulon %>%
        remove_rownames %>%
        tibble::column_to_rownames(var = 'tf')

      tf_activities <- tf_activities %>%
        filter(z_score > 0)

      output_df <- create_empty_CTR_dataframe()

      for (row in 1:nrow(tf_activities)) {

        r_tf <- create_empty_CTR_dataframe()

        tf_l <- create_empty_CTR_dataframe()

        #if (tf_activities[row, "z_score"] > 0) {
        tf <- as.character(tf_activities[row,]["gene"])
        targets <- sorted_regulon[tf,][1]
        receptors <- R2TF[tf,][1]
        tf_ligands <- intersect(targets[[1]], ligands)
        if (length(tf_ligands) > 0) {
          for (ligand in tf_ligands) {
            expressed <- FALSE
            translations <- ligand
            if (length(translations) > 0) {
              for (l in translations) {
                if (l %in% rownames(gene_expression)) {
                  ex_value <- gene_expression[l, tf_activities[row, "cluster"]]
                  if (ex_value != 0) {
                    expressed <- TRUE
                  }
                }
              }
            }

            if (expressed == TRUE) {
              df <- add_entry_to_CTR_dataframe(tf_activities[row, "cluster"],
                                               tf_activities[row, "cluster"],
                                               tf_activities[row, "gene"],
                                               ligand,
                                               'Transcription Factor',
                                               'Ligand',
                                               tf_activities[row, "z_score"])
              tf_l <- rbind(tf_l, df)
            }
          }
        }
        if (length(receptors[[1]]) > 0) {
          for (receptor in receptors[[1]]) {
            df <- add_entry_to_CTR_dataframe(tf_activities[row, "cluster"],
                                             tf_activities[row, "cluster"],
                                             receptor,
                                             tf_activities[row, "gene"],
                                             'Receptor',
                                             'Transcription Factor',
                                             tf_activities[row, "z_score"])
            r_tf <- rbind(r_tf, df)
          }
        }
        #}

        r_tf$gene_A <- gsub("_", "+", r_tf$gene_A, fixed = TRUE)
        r_tf$gene_B <- gsub("_", "+", r_tf$gene_B, fixed = TRUE)
        tf_l$gene_A <- gsub("_", "+", tf_l$gene_A, fixed = TRUE)
        tf_l$gene_B <- gsub("_", "+", tf_l$gene_B, fixed = TRUE)

        output_df <- rbind(output_df, r_tf)
        output_df <- rbind(output_df, tf_l)
      }
    } else {
      output_df = NA
    }
    return(output_df)
  }

#' Generate connections in intracellular network
#'
#' This function loads the transcription factor activity table/data frame
#' for multiple cell types, and generates a table containing all detected intracellular connections.
#'
#' @param tf_activities Data Frame with transcription factor activities by cell type
#' @param confidence_level Confidence Level for the DoRothEA regulon used to get the transcription factor activities
#' @param gene_activities Table with average gene expression levels
#' @return A data frame with CrossTalkeR input
#' @import dplyr
#' @import tibble
#' @export
generate_intracellular_network <-
  function(tf_activities,
           gene_expression,
           regulon,
           organism = "human") {
    if (dim(tf_activities)[1] > 0) {
      if (any(tf_activities$z_score > 0)) {
        if (organism == "human") {
          R2TF <- aggregate(RTF_DB_2$receptor ~ RTF_DB_2$tf, FUN = c)
          colnames(R2TF) <- c('tf', 'receptors')
          R2TF <- R2TF %>%
            tibble::remove_rownames() %>%
            tibble::column_to_rownames(var = 'tf')
        } else {
          R2TF <- aggregate(RTF_DB_mouse$receptor ~ RTF_DB_mouse$tf, FUN = c)
          colnames(R2TF) <- c('tf', 'receptors')
          R2TF <- R2TF %>%
            tibble::remove_rownames() %>%
            tibble::column_to_rownames(var = 'tf')
        }

        regulon <- regulon %>%
          rename(tf = source)
        sorted_regulon <-
          aggregate(regulon$target ~ regulon$tf, FUN = c)
        colnames(sorted_regulon) <- c('tf', 'targets')
        sorted_regulon <- sorted_regulon %>%
          tibble::remove_rownames() %>%
          tibble::column_to_rownames(var = 'tf')

        #recept_regulon <- create_empty_Regulon_dataframe()

        tf_activities <- tf_activities %>%
          filter(z_score > 0)

        RTF_df <- data.frame(
          Receptor = character(),
          TF = character()
        )

        TFTG_df <- data.frame(
          celltype = character(),
          TF = character(),
          Target_Gene = character(),
          TF_Score = numeric()
        )

        for (row in 1:nrow(tf_activities)) {
          #if (tf_activities[row, "z_score"] > 0) {
          tf <- as.character(tf_activities[row, "gene"])
          targets <- sorted_regulon[tf,][[1]]
          receptors <- R2TF[tf,][1]
          if (length(targets) > 0) {
            if (length(receptors[[1]]) > 0) {
              for (target in targets) {
                expressed <- FALSE
                if (target %in% rownames(gene_expression)) {
                  ex_value <- gene_expression[target, tf_activities[row, "cluster"]]
                  if (ex_value != 0) {
                    expressed <- TRUE
                  }
                }
                if (expressed == TRUE) {
                  TFTG_tmp <-
                    data.frame(
                      tf_activities[row, "cluster"],
                      tf_activities[row, "gene"],
                      target,
                      tf_activities[row, "z_score"]
                    )
                  names(TFTG_tmp) <-
                    c(
                      'celltype',
                      'TF',
                      'Target_Gene',
                      'TF_Score'
                    )
                  TFTG_df <- rbind(TFTG_df, TFTG_tmp)
                }
              }
              for (receptor in receptors[[1]]) {
                RTF_tmp <-
                  data.frame(
                    receptor,
                    tf_activities[row, "gene"]
                  )
                names(RTF_tmp) <-
                  c(
                    'Receptor',
                    'TF'
                  )
                RTF_df <- rbind(RTF_df, RTF_tmp)
              }
            }
          }
          #}
        }
        recept_regulon <- merge(x = RTF_df, y = TFTG_df,
                                by = "TF", all = TRUE)

      } else {
        recept_regulon = NA
      }
      return(recept_regulon)
    }
  }


