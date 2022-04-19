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
#' @import OmnipathR
#' @import dorothea
#' @export
generate_CrossTalkeR_input_significant_table <-
  function(tf_activities,
           confidence_level = c("A", "B", "C"),
           gene_expression) {

    dorothea_regulon_human <-
      get(data("dorothea_hs", package = "dorothea"))
    regulon <- dorothea_regulon_human %>%
      dplyr::filter(confidence %in% confidence_level)

    intercell <-
      import_intercell_network(
        receiver_param = list(categories = c('receptor')),
        transmitter_param = list(categories = c('ligand'))
      )
    ligands = list(intercell$source_genesymbol)

    posttranslational <-
      import_post_translational_interactions(organism = '9606',
                                             genesymbol = TRUE)
    posttranslational = aggregate(posttranslational$source_genesymbol ~ posttranslational$target_genesymbol,
                                  FUN = c)
    colnames(posttranslational) <- c('tf', 'receptors')
    posttranslational = posttranslational %>%
      remove_rownames %>%
      tibble::column_to_rownames(var = 'tf')

    sorted_regulon <-
      aggregate(regulon$target ~ regulon$tf, FUN = c)
    colnames(sorted_regulon) <- c('tf', 'targets')
    sorted_regulon = sorted_regulon %>%
      remove_rownames %>%
      tibble::column_to_rownames(var = 'tf')

    output_df = data.frame(
      source = character(),
      target = character(),
      gene_A = character(),
      gene_B = character(),
      type_gene_A = character(),
      type_gene_B = character(),
      MeanLR = numeric()
    )

    for (row in 1:nrow(tf_activities)) {

      r_tf = data.frame(
        source = character(),
        target = character(),
        gene_A = character(),
        gene_B = character(),
        type_gene_A = character(),
        type_gene_B = character(),
        MeanLR = double()
      )

      tf_l = data.frame(
        source = character(),
        target = character(),
        gene_A = character(),
        gene_B = character(),
        type_gene_A = character(),
        type_gene_B = character(),
        MeanLR = double()
      )

      if (tf_activities[row, "z_score"] > 0) {
        tf = as.character(tf_activities[row, "gene"])
        targets = sorted_regulon[tf,][1]
        receptors = posttranslational[tf,][1]
        tf_ligands = intersect(targets[[1]], ligands[[1]])
        if (length(tf_ligands) > 0) {
          for (ligand in tf_ligands) {
            expressed = FALSE
            if (ligand %in% rownames(gene_expression)) {
              ex_value = gene_expression[ligand, tf_activities[row, "cluster"]]
              if (ex_value != 0) {
                expressed = TRUE
              }
            }

            if (expressed == TRUE) {
              df <-
                data.frame(
                  tf_activities[row, "cluster"],
                  tf_activities[row, "cluster"],
                  tf_activities[row, "gene"],
                  ligand,
                  'Transcription Factor',
                  'Ligand',
                  tf_activities[row, "z_score"]
                )
              names(df) <-
                c(
                  'source',
                  'target',
                  'gene_A',
                  'gene_B',
                  'type_gene_A',
                  'type_gene_B',
                  "MeanLR"
                )
              tf_l = rbind(tf_l, df)
            }
          }
        }
        if (length(receptors[[1]]) > 0) {
          for (receptor in receptors[[1]]) {
            df <-
              data.frame(
                tf_activities[row, "cluster"],
                tf_activities[row, "cluster"],
                receptor,
                tf_activities[row, "gene"],
                'Receptor',
                'Transcription Factor',
                tf_activities[row, "z_score"]
              )
            names(df) <-
              c(
                'source',
                'target',
                'gene_A',
                'gene_B',
                'type_gene_A',
                'type_gene_B',
                "MeanLR"
              )
            r_tf = rbind(r_tf, df)
          }
        }
      }

      r_tf$gene_A = gsub("_", "+", r_tf$gene_A, fixed = TRUE)
      r_tf$gene_B = gsub("_", "+", r_tf$gene_B, fixed = TRUE)
      tf_l$gene_A = gsub("_", "+", tf_l$gene_A, fixed = TRUE)
      tf_l$gene_B = gsub("_", "+", tf_l$gene_B, fixed = TRUE)

      output_df = rbind(output_df, r_tf)
      output_df = rbind(output_df, tf_l)
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
#' @import OmnipathR
#' @import dorothea
#' @export

generate_CrossTalkeR_input_mouse_significant_table <-
  function(tf_activities,
           confidence_level = c("A", "B", "C"),
           gene_expression) {

    dorothea_regulon_human <-
      get(data("dorothea_hs", package = "dorothea"))
    regulon <- dorothea_regulon_human %>%
      dplyr::filter(confidence %in% confidence_level)

    intercell <-
      import_intercell_network(
        receiver_param = list(categories = c('receptor')),
        transmitter_param = list(categories = c('ligand'))
      )
    ligands = list(intercell$source_genesymbol)

    posttranslational <-
      import_post_translational_interactions(organism = '9606',
                                             genesymbol = TRUE)
    posttranslational = aggregate(posttranslational$source_genesymbol ~ posttranslational$target_genesymbol,
                                  FUN = c)
    colnames(posttranslational) <- c('tf', 'receptors')
    posttranslational = posttranslational %>%
      remove_rownames %>%
      tibble::column_to_rownames(var =
                                   'tf')

    sorted_regulon <-
      aggregate(regulon$target ~ regulon$tf, FUN = c)
    colnames(sorted_regulon) <- c('tf', 'targets')
    sorted_regulon = sorted_regulon %>%
      remove_rownames %>%
      tibble::column_to_rownames(var =
                                   'tf')

    output_df = data.frame(
      source = character(),
      target = character(),
      gene_A = character(),
      gene_B = character(),
      type_gene_A = character(),
      type_gene_B = character(),
      MeanLR = numeric()
    )

    for (row in 1:nrow(tf_activities)) {

      r_tf = data.frame(
        source = character(),
        target = character(),
        gene_A = character(),
        gene_B = character(),
        type_gene_A = character(),
        type_gene_B = character(),
        MeanLR = double()
      )

      tf_l = data.frame(
        source = character(),
        target = character(),
        gene_A = character(),
        gene_B = character(),
        type_gene_A = character(),
        type_gene_B = character(),
        MeanLR = double()
      )

      if (tf_activities[row, "z_score"] > 0) {
        tf = as.character(tf_activities[row,]["gene"])
        targets = sorted_regulon[tf,][1]
        receptors = posttranslational[tf,][1]
        tf_ligands = intersect(targets[[1]], ligands[[1]])
        if (length(tf_ligands) > 0) {
          for (ligand in tf_ligands) {
            expressed = FALSE
            translations = translate_to_mouse_ligands(ligand)
            if (length(translations) > 0) {
              for (l in translations) {
                if (l %in% rownames(gene_expression)) {
                  ex_value = gene_expression[l, tf_activities[row, "cluster"]]
                  if (ex_value != 0) {
                    expressed = TRUE
                  }
                }
              }
            }

            if (expressed == TRUE) {
              df <-
                data.frame(
                  tf_activities[row, "cluster"],
                  tf_activities[row, "cluster"],
                  tf_activities[row, "gene"],
                  ligand,
                  'Transcription Factor',
                  'Ligand',
                  tf_activities[row, "z_score"]
                )
              names(df) <-
                c(
                  'source',
                  'target',
                  'gene_A',
                  'gene_B',
                  'type_gene_A',
                  'type_gene_B',
                  "MeanLR"
                )
              tf_l = rbind(tf_l, df)
            }
          }
        }
        if (length(receptors[[1]]) > 0) {
          for (receptor in receptors[[1]]) {
            df <-
              data.frame(
                tf_activities[row, "cluster"],
                tf_activities[row, "cluster"],
                receptor,
                tf_activities[row, "gene"],
                'Receptor',
                'Transcription Factor',
                tf_activities[row, "z_score"]
              )
            names(df) <-
              c(
                'source',
                'target',
                'gene_A',
                'gene_B',
                'type_gene_A',
                'type_gene_B',
                "MeanLR"
              )
            r_tf = rbind(r_tf, df)
          }
        }
      }

      r_tf$gene_A = gsub("_", "+", r_tf$gene_A, fixed = TRUE)
      r_tf$gene_B = gsub("_", "+", r_tf$gene_B, fixed = TRUE)
      tf_l$gene_A = gsub("_", "+", tf_l$gene_A, fixed = TRUE)
      tf_l$gene_B = gsub("_", "+", tf_l$gene_B, fixed = TRUE)

      output_df = rbind(output_df, r_tf)
      output_df = rbind(output_df, tf_l)
    }

    return(output_df)
  }