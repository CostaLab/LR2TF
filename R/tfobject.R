#'
#'Object structure for the results of the TF activity analyses
#'
#'@slot tf_activities All Cell Cell Interaction Networks
#'@slot average_gene_expression All tables from single condition
#'@slot regulon  Max meanLR from all
#'@slot CTR_input All Celltype in the experiment
#'@slot intracellular_network  Cell Cell Interaction Plots
tfobject <- setClass("TFObj",
                     slots = list(tf_activities_condition = "list",
                                  tf_activities_cluster = "list",
                                  average_gene_expression = "list",
                                  regulon = "data.frame",
                                  CTR_input_condition = "list",
                                  CTR_input_cluster = "list",
                                  intracellular_network_condition = "list",
                                  intracellular_network_cluster = "list"
                     )
)