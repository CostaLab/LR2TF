#' Saves most variable transcription factor activity scores per cell type into a table
#'
#' This function saves the transcription factor activity scores per cell type
#' into a csv table.
#'
#' @param tf_scores data frame with transcription factor activity scores per cell type
#' @param condition Sample condition for file naming(e.g. control, disease ...)
#' @param out_path Output path to save results
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import stringr
#' @export
save_variable_tf_scores <- function(tf_scores, condition, out_path) {
  highly_variable_tfs_all <- tf_scores %>%
    group_by(tf) %>%
    mutate(var = var(avg)) %>%
    ungroup() %>%
    distinct(tf)

  summarized_viper_scores_df_all <- tf_scores %>%
    semi_join(highly_variable_tfs_all, by = "tf") %>%
    dplyr::select(-std) %>%
    spread(tf, avg) %>%
    data.frame(row.names = 1, check.names = FALSE)
  tf_scores <- t(summarized_viper_scores_df_all)
  write.csv(tf_scores, file = paste0(out_path, "/variable_tf_scores", "_", condition, ".csv"))

  return(tf_scores)
}

#' Adds gene type to the name of the gene
#'
#' Description
#'
#' @param df dataframe with all interactions
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @export
add_node_type <- function(df) {
  df <- df %>%
    mutate(gene_A = ifelse(type_gene_A == "Ligand", paste0(gene_A, "|L"), gene_A))
  df <- df %>%
    mutate(gene_A = ifelse(type_gene_A == "Receptor", paste0(gene_A, "|R"), gene_A))
  df <- df %>%
    mutate(gene_A = ifelse(type_gene_A == "Transcription Factor", paste0(gene_A, "|TF"), gene_A))
  df <- df %>%
    mutate(gene_B = ifelse(type_gene_B == "Ligand", paste0(gene_B, "|L"), gene_B))
  df <- df %>%
    mutate(gene_B = ifelse(type_gene_B == "Receptor", paste0(gene_B, "|R"), gene_B))
  df <- df %>%
    mutate(gene_B = ifelse(type_gene_B == "Transcription Factor", paste0(gene_B, "|TF"), gene_B))
}

#' Combining Ligand-Receptor interaction prediction with Transcription Factor interaction predictions
#'
#' Description
#'
#' @param tf_table table with tf interactions
#' @param LR_prediction path to or dataframe with ligand-receptor interaction prediction
#' @param out_path path to save results
#' @param condition sample condition of data
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @export
combine_LR_and_TF <- function(tf_table, LR_prediction, out_path, condition, add_node_type = FALSE) {
  if (!is.data.frame(LR_prediction)) {
    lr_table <- read.csv(LR_prediction)
    row.names(lr_table) <- lr_table$X
    lr_table$X <- NULL
  } else {
    lr_table <- LR_prediction
  }

  intra_connections <- tf_table[NULL, ]
  for (celltype in unique(append(lr_table$source, lr_table$target))) {
    lr_filtered_ligands <- lr_table[lr_table$source == celltype, ]
    lr_filtered_receptors <- lr_table[lr_table$target == celltype, ]
    lr_ligands <- unique(lr_filtered_ligands$gene_A)
    lr_receptors <- unique(lr_filtered_receptors$gene_B)
    tf_table_receptors <- tf_table[tf_table$target == celltype & tf_table$type_gene_A == "Receptor", ]
    tf_table_ligands <- tf_table[tf_table$source == celltype & tf_table$type_gene_B == "Ligand", ]
    tf_receptor_interactions <- tf_table_receptors %>%
      filter(gene_A %in% lr_receptors)
    tf_ligand_interactions <- tf_table_ligands %>%
      filter(gene_B %in% lr_ligands)
    intra_connections <- rbind(intra_connections, tf_receptor_interactions, tf_ligand_interactions)
  }
  intra_connections$all_pair <- paste0(
    intra_connections$source, "/",
    intra_connections$gene_A, "/",
    intra_connections$target, "/",
    intra_connections$gene_B
  )
  intra_connections <- intra_connections[!duplicated(intra_connections$all_pair), ]
  intra_connections$all_pair <- NULL
  complete_interactions <- rbind(intra_connections, lr_table)
  if (add_node_type) {
    complete_interactions <- add_node_type(complete_interactions)
  }
  write.csv(complete_interactions, paste0(out_path, "CrossTalkeR_input_", condition, ".csv"), row.names = FALSE)
  return(complete_interactions)
}

#' Combining Ligand-Receptor interaction prediction with Transcription Factor interaction predictions considering receptor complexes
#'
#' Description
#'
#' @param tf_table table with tf interactions
#' @param LR_prediction path to or dataframe with ligand-receptor interaction prediction
#' @param out_path path to save results
#' @param condition sample condition of data
#' @return complete_interactions table with all interactions
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @export
combine_LR_and_TF_complexes <- function(tf_table, LR_prediction, out_path, condition, add_node_type = FALSE) {
  if (!is.data.frame(LR_prediction)) {
    lr_table <- read.csv(LR_prediction)
    row.names(lr_table) <- lr_table$X
    lr_table$X <- NULL
  } else {
    lr_table <- LR_prediction
  }

  intra_connections <- tf_table[NULL, ]
  for (celltype in unique(append(lr_table$source, lr_table$target))) {
    lr_filtered_ligands <- lr_table[lr_table$source == celltype, ]
    lr_filtered_receptors <- lr_table[lr_table$target == celltype, ]
    lr_ligands <- unique(lr_filtered_ligands$gene_A)
    lr_receptors <- unique(lr_filtered_receptors$gene_B)

    contains_complex <- grepl("_", lr_receptors)
    R_with_complex <- lr_receptors[contains_complex]
    R_without_complex <- lr_receptors[!contains_complex]

    tf_table_receptors <- tf_table[tf_table$target == celltype & tf_table$type_gene_A == "Receptor", ]
    tf_receptor_interactions <- tf_table_receptors %>%
      filter(gene_A %in% R_without_complex)

    complex_df <- tf_table[0, ]
    if (length(R_with_complex) > 0) {
      for (i in 1:length(R_with_complex)) {
        complex <- R_with_complex[i]
        receptors <- strsplit(complex, "_")[[1]]
        R_TF_with_complex <- tf_table_receptors[tf_table_receptors$gene_A %in% receptors, ]
        if (nrow(R_TF_with_complex) == 0) {
          next
        }
        unique(R_TF_with_complex)
        R_TF_with_complex$gene_A <- complex
        complex_df <- rbind(complex_df, R_TF_with_complex)
      }
      complex_df <- complex_df[!duplicated(complex_df), ]
    }
    tf_receptor_interactions <- rbind(tf_receptor_interactions, complex_df)

    tf_table_ligands <- tf_table[tf_table$source == celltype & tf_table$type_gene_B == "Ligand", ]
    tf_ligand_interactions <- tf_table_ligands %>%
      filter(gene_B %in% lr_ligands)
    intra_connections <- rbind(intra_connections, tf_receptor_interactions, tf_ligand_interactions)
  }
  intra_connections$all_pair <- paste0(
    intra_connections$source, "/",
    intra_connections$gene_A, "/",
    intra_connections$target, "/",
    intra_connections$gene_B
  )
  intra_connections <- intra_connections[!duplicated(intra_connections$all_pair), ]
  intra_connections$all_pair <- NULL
  complete_interactions <- rbind(intra_connections, lr_table)
  if (add_node_type) {
    complete_interactions <- add_node_type(complete_interactions)
  }
  write.csv(complete_interactions, paste0(out_path, "CrossTalkeR_input_", condition, ".csv"), row.names = FALSE)
  return(complete_interactions)
}

#' Combining Ligand-Receptor interaction prediction with Transcription Factor interaction predictions
#'
#' Description
#'
#' @param tf_table table with tf interactions
#' @param LR_path path to ligand receptor interaction prediction
#' @param out_path path to save results
#' @param condition sample condition of data
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @export
combine_LR_and_TF_unfiltered <- function(tf_table, LR_path, out_path, condition) {
  if (!is.data.frame(LR_prediction)) {
    lr_table <- read.csv(LR_prediction)
    row.names(lr_table) <- lr_table$X
    lr_table$X <- NULL
  } else {
    lr_table <- LR_prediction
  }

  complete_interactions <- rbind(tf_table, lr_table)
  complete_interactions <- add_node_type(complete_interactions)

  write.csv(complete_interactions, paste0(out_path, "CrossTalkeR_input_", condition, ".csv"), row.names = FALSE)
}

#' Create an empty dataframe with CrossTalkeR input table format
#'
#' @export
create_empty_CTR_dataframe <- function() {
  empty_df <- data.frame(
    source = character(),
    target = character(),
    gene_A = character(),
    gene_B = character(),
    type_gene_A = character(),
    type_gene_B = character(),
    MeanLR = numeric()
  )
  return(empty_df)
}

#' Add entry to a dataframe in the CrossTalkeR input format
#'
#' @export
add_entry_to_CTR_dataframe <- function(source, target, gene_A, gene_B, type_gene_A, type_gene_B, MeanLR) {
  df <-
    data.frame(
      source,
      target,
      gene_A,
      gene_B,
      type_gene_A,
      type_gene_B,
      MeanLR
    )
  names(df) <-
    c(
      "source",
      "target",
      "gene_A",
      "gene_B",
      "type_gene_A",
      "type_gene_B",
      "MeanLR"
    )

  return(df)
}

#' Create an empty dataframe with CrossTalkeR input table format
#'
#' @export
create_empty_Regulon_dataframe <- function() {
  empty_df <- data.frame(
    celltype = character(),
    Receptor = character(),
    TF = character(),
    Target_Gene = character(),
    TF_Score = numeric()
  )
  return(empty_df)
}

#' Add entry to a dataframe in the CrossTalkeR input format
#'
#' @export
add_entry_to_Regulon_dataframe <- function(celltype, Receptor, TF, Target_Gene, TF_Score) {
  df <-
    data.frame(
      celltype,
      Receptor,
      TF,
      Target_Gene,
      TF_Score
    )
  names(df) <-
    c(
      "celltype",
      "Receptor",
      "TF",
      "Target_Gene",
      "TF_Score"
    )

  return(df)
}

#' Convert Seurat object to anndata object and save anndata object file
#'
#' @param seuratobject Input Seurat Object
#' @param out_path Output path to save results
#' @import sceasy
#' @export
convert_seurat_to_anndata <- function(seuratobject, out_path) {
  sceasy::convertFormat(seuratobject, from = "seurat", to = "anndata", outFile = paste0(out_path, "anndata_object.h5ad"))
}

#' Check arguments passed by User for validity
#'
#' @param arguments_list list of user defined arguments
#' @return val_arguments validated list of arguments
#' @export
validate_input_arguments <- function(arguments_list) {
  if (is.null(arguments_list$out_path)) {
    print("Please provide an output path")
  } else {
    if (substring(arguments_list$out_path, length(arguments_list$out_path) - 1, length(arguments_list$out_path)) != "/") {
      arguments_list$out_path <- paste0(arguments_list$out_path, "/")
    }
  }
  if (is.null(arguments_list$celltype)) {
    print("Please provide the name of the metadata field containing cell type annotations")
  }
  if (is.null(arguments_list$condition)) {
    print("Please provide the name of the metadata field containing condition annotations")
  }
  if (is.null(arguments_list$organism)) {
    arguments_list$organism <- "human"
  }
  if (is.null(arguments_list$comparison_list)) {
    arguments_list$comparison_list <- NA
  }
  if (is.null(arguments_list$logfc)) {
    arguments_list$logfc <- 0.0
  }
  if (is.null(arguments_list$pval)) {
    arguments_list$pval <- 0.05
  }
  if (is.null(arguments_list$num_cell_filter)) {
    arguments_list$num_cell_filter <- 0
  }
  if (is.null(arguments_list$reg)) {
    arguments_list$reg <- load_dorothea_regulon(arguments_list$organism)
  } else {
    if (typeof(arguments_list$reg) == "character") {
      arguments_list$reg <- read.csv(arguments_list$reg, header = TRUE)
    }
    if (!all(c("source", "target", "weight") %in% names(arguments_list$reg))) {
      stop("Not all necessary columns found in regulon table! Please make sure that the regulon has the columns source, target and weight!")
    }
  }
  if (is.null(arguments_list$plot)) {
    arguments_list$plot <- TRUE
  } else {
    if (!is.boolean(arguments_list$plot)) {
      stop("Plot argument must be a boolean value!")
    }
  }
  if (is.null(arguments_list$test_type)) {
    arguments_list$test_type <- "t"
  } else {
    if (!any(arguments_list$test_type %in% c("binom", "t", "wilcox"))) {
      arguments_list$test_type <- "binom"
    }
  }
  if (is.null(arguments_list$Seurat)) {
    arguments_list$Seurat <- FALSE
  }
  return(arguments_list)
}
