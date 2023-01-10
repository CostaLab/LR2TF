#' Mouse-Human homologous genes table transcription factors
#'
#' A data frame containing the mouse-human homologous genes of the dorothea regulon
#'
#' @format A data frame with 1026 rows and 4 variables:
#' \itemize{
#'   \item 10090: Mouse gene symbol
#'   \item 9606: Human gene symbol
#'   \item 10090_ID: Mouse NCBI gene ID
#'   \item 9606_ID: Human NCBI gene ID
#' }
"mouse_regulon_homologous"

#' Mouse-Human homologous genes table ligands
#'
#' A data frame containing the mouse-human homologous genes of the ligands in the Omnipath database
#'
#' @format A data frame with 1026 rows and 4 variables:
#' \itemize{
#'   \item 10090: Mouse gene symbol
#'   \item 9606: Human gene symbol
#' }
"ligands_homologous"

#' Example dataset
#'
#' Example dataset for tutorial execution of the package
#'
#' @export bone_marrow_stromal_cell_example
"bone_marrow_stromal_cell_example"

#' Dataframe with pre-computed receptor to transcription faactor connections
#'
#' Dataframe with pre-computed receptor to transcription faactor connections from Omnipath
#'
"RTF_DB"