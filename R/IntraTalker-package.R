#' @keywords internal
"_PACKAGE"

#' @importFrom grDevices dev.off pdf
#' @importFrom methods new
#' @importFrom stats aggregate dist na.omit sd var
#' @importFrom utils read.csv write.csv
NULL

# Column names used for non-standard evaluation in dplyr/maditr pipelines, and
# lazy-loaded package datasets referenced by name.
utils::globalVariables(c(
  "CellType", "FDR", "activity", "avg", "avg_log2FC", "cell", "cell_type",
  "condition", "gene_A", "gene_B", "ligands_human", "ligands_mouse", "logFC",
  "p_val_adj", "std", "tf", "type_gene_A", "type_gene_B", "var", "z_score"
))
