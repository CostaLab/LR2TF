#' Translates human gene symbols to homologous mouse gene symbols
#'
#' This function translates human gene symbols to homologous mouse gene symbols
#'
#' @param df data frame with the transcription factor activities and the gene names as rownames
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @export
translate_to_mouse <- function(df){
  regulon_data = LR2TF::mouse_regulon_homologous

  to_translate = rownames(df)
  translated <- regulon_data[regulon_data$'9606' %in% to_translate ,]

  df <- rownames_to_column(df, '9606')

  result_df <- inner_join(df, translated)
  rownames(result_df) <- result_df[,'10090']
  result_df$gene <- result_df[,'10090']
  result_df$'9606' = NULL
  result_df$'10090' = NULL

  return(result_df)
}

#' Translates mouse gene symbols to homologous human gene symbols
#'
#' This function translates mouse gene symbols to homologous human gene symbols
#'
#' @param df data frame with the transcription factor activities and the gene names as rownames
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @export
translate_to_human_old <- function(df){
  data("mouse_regulon_homologous", envir=environment())
  regulon_data = LR2TF::mouse_regulon_homologous

  to_translate = rownames(df)
  translated <- regulon_data[regulon_data$'10090' %in% to_translate ,]

  df <- rownames_to_column(df, '10090')

  result_df <- inner_join(df, translated)
  result_df <- result_df[!duplicated(result_df$'9606'),]
  rownames(result_df) <- result_df[,'9606']
  result_df$gene <- result_df[,'9606']
  result_df$'9606' = NULL
  result_df$'10090' = NULL

  return(result_df)
}

#' Translates mouse gene symbols to homologous human gene symbols
#'
#' This function translates mouse gene symbols to homologous human gene symbols
#'
#' @param df data frame with the transcription factor activities and the gene names as rownames
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @export
translate_to_human <- function(df){
  data("mouse_regulon_homologous", envir=environment())
  regulon_data = LR2TF::mouse_regulon_homologous

  to_translate = df$gene
  translated <- regulon_data[regulon_data$'10090' %in% to_translate ,]

  colnames(df)[colnames(df) == 'gene'] <- '10090'

  result_df <- inner_join(df, translated)
  result_df <- result_df[!duplicated(result_df$'9606'),]
  result_df$gene <- result_df[,'9606']
  result_df$'9606' = NULL
  result_df$'10090' = NULL

  rownames(result_df) <- paste(result_df$cluster, result_df$gene, sep = ".")

  return(result_df)
}


#' Translates mouse gene symbols to homologous human gene symbols for ligands
#'
#' This function translates mouse gene symbols to homologous human gene symbols
#'
#' @param ligand_list list with ligand genesymbols to be translated
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @export
translate_to_human_ligands <- function(ligand_list){
  data("ligands_homologous", envir=environment())
  ligand_data = LR2TF::ligands_homologous

  translated <- ligand_data[ligand_data$'10090' %in% ligand_list ,]
  result_list = translated$'9606'

  return(result_list)
}

#' Translates human gene symbols to homologous mouse gene symbols for ligands
#'
#' This function translates human gene symbols to homologous mouse gene symbols
#'
#' @param ligand_list list with ligand genesymbols to be translated
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @export
translate_to_mouse_ligands <- function(ligand_list){
  data("ligands_homologous", envir=environment())
  ligand_data = LR2TF::ligands_homologous

  translated <- ligand_data[ligand_data$'9606' %in% ligand_list ,]
  result_list = translated$'10090'

  return(result_list)
}
