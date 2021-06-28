import dorothea_based_r2tf2l as dba
import network_utils as nu

if __name__ == '__main__':

    data_dir = "~/R_results/results_complete_ABC"
    dorothea_regulon_file = data_dir + "/dorothea_regulon_ABC.csv"
    conditions = ["Control", "PMF_MF2"]
    out_dir = "~/R_results/interaction_analysis_results_ABC"
    dba.get_ligand_targets(data_dir=data_dir, dorothea_regulon_file=dorothea_regulon_file, conditions=conditions,
                           out_dir=out_dir)

    Control_lr_path = "~/R_results/ControlLR.csv"
    PMF_lr_path = "~/R_results/PMFLR.csv"

    Control_rtf_path = out_dir + "/Control_r2tf_network.csv"
    PMF_rtf_path = out_dir + "/PMF_MF2_r2tf_network.csv"

    Control_tfl_path = out_dir + "/Control_tf2ligand_network.csv"
    PMF_tfl_path = out_dir + "/PMF_MF2_tf2ligand_network.csv"

    nu.network_from_tables(Control_lr_path, Control_rtf_path, Control_tfl_path, out_dir, "Control")
    nu.network_from_tables(PMF_lr_path, PMF_rtf_path, PMF_tfl_path, out_dir, "PMF")

    nu.convert_table_to_crosstalker_input(Control_lr_path, Control_rtf_path, Control_tfl_path, out_dir, "Control")
    nu.convert_table_to_crosstalker_input(PMF_lr_path, PMF_rtf_path, PMF_tfl_path, out_dir, "PMF")
