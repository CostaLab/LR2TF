import os
import omnipath as om
import pandas as pd
from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML


def get_ligands(target_list, lr):
    ligand_interactions = lr[lr['genesymbol_intercell_source'].isin(target_list)]
    ligand_list = list(set(ligand_interactions['genesymbol_intercell_source'].tolist()))
    return ligand_list


def read_input_filenames(conditions, data_path):
    files_dict = {}
    for condition in conditions:
        condition_files = {}
        for path, subdirs, files in os.walk(data_path):
            for name in files:
                if "cluster_per" in name and "png" in name:
                    condition_files["cluster_comparison_image"] = os.path.join(path, name)
                if condition in name:
                    if "average_gene_expression" in name:
                        condition_files["average_gene_expression"] = os.path.join(path, name)
                    elif "tf_scores" in name:
                        condition_files["tf_activity_scores"] = os.path.join(path, name)
                    elif "all_specificmarker" in name:
                        condition_files["all_specific_markers"] = os.path.join(path, name)
                    elif "filtered" in name:
                        condition_files["filtered_specific_markers"] = os.path.join(path, name)
                    elif "top20" in name and "png" in name:
                        condition_files["heatmap_image"] = os.path.join(path, name)
                    elif "tf_activity_cluster" in name and "png" in name:
                        condition_files["cluster_image"] = os.path.join(path, name)

        files_dict[condition] = condition_files
    print("File names loaded")
    return files_dict


def get_ligand_targets(data_dir, dorothea_regulon_file, out_dir, conditions):
    os.makedirs(out_dir)
    condition_mapping = {"control": "Control", "condition": "PMF,MF2"}

    lr = om.interactions.import_intercell_network(transmitter_params={'categories': 'ligand'},
                                                  receiver_params={'categories': 'receptor'})
    r2tf = om.interactions.PostTranslational.get(organism='human', genesymbol=True)  # ppi

    files_per_condition = read_input_filenames(conditions, data_dir)

    overview_data_per_condition = {}
    cell_data_per_condition = {}

    for condition in conditions:
        file_dict = files_per_condition[condition]
        r2tf_element_list = []

        gene_expression = pd.read_csv(file_dict["average_gene_expression"], sep=',', header=0)
        gene_expression = gene_expression.rename(columns={gene_expression.columns.tolist()[0]: 'gene'})
        tf_activity = pd.read_csv(file_dict["tf_activity_scores"], sep=',', header=0, index_col=0)
        dorothea_regulon = pd.read_csv(dorothea_regulon_file, sep=',', header=0)
        specific_marker = pd.read_csv(file_dict["all_specific_markers"], sep=',', header=0, index_col=0)
        cell_types = tf_activity.columns.tolist()

        grouped_regulon = dorothea_regulon.groupby('tf')['target'].apply(list).reset_index(name='targets')
        grouped_regulon = grouped_regulon.set_index('tf')
        regulon_dict = grouped_regulon['targets'].to_dict()

        # specific_marker_dict = specific_marker.set_index('gene').transpose().to_dict(orient='dict')
        gene_expression_dict = gene_expression.set_index('gene').transpose().to_dict(orient='dict')
        tf_activity_dict = tf_activity.transpose().to_dict(orient='dict')
        number_cell_types = len(cell_types)
        number_active_tfs = len(tf_activity_dict.keys())
        all_possible_targets = []
        all_possible_ligands = []
        all_expressed_targets = []
        all_expressed_ligands = []

        tf_csv_data = []
        cell_data = []
        tfs_for_receptors = []
        for cell in cell_types:

            cell_fix_name = ""

            if cell == "Myeloid_Granulocyte":
                cell_fix_name = "Myeloid"
            else:
                cell_fix_name = cell

            cell_possible_targets = []
            cell_possible_ligands = []
            cell_expressed_targets = []
            cell_expressed_ligands = []
            cell_number_active_tfs = 0
            tf_data = {}

            cell_specific_markers = specific_marker[specific_marker['cluster'].isin([cell])]
            specific_marker_dict = cell_specific_markers.set_index('gene').transpose().to_dict(orient='dict')
            for key in tf_activity_dict:
                if tf_activity_dict[key][cell] > 0:
                    cell_number_active_tfs = cell_number_active_tfs + 1

                    if "NKX" in key:
                        key = key.replace('.', '-')
                    possible_targets = regulon_dict[key]
                    all_possible_targets.extend(possible_targets)
                    cell_possible_targets.extend(possible_targets)

                    ligands = get_ligands(possible_targets, lr)
                    all_possible_ligands.extend(ligands)
                    cell_possible_ligands.extend(ligands)

                    tf_targets = gene_expression[gene_expression['gene'].isin(possible_targets)]
                    expressed_tf_targets = tf_targets.query("(" + cell + " < 0) or (" + cell + " > 0)")
                    expressed_tf_targets_list = expressed_tf_targets['gene'].tolist()
                    expressed_ligands = get_ligands(expressed_tf_targets_list, lr)

                    if key in specific_marker_dict:
                        p_value = str(specific_marker_dict[key]['p_val'])
                    else:
                        continue
                    all_expressed_targets.extend(expressed_tf_targets_list)
                    cell_expressed_targets.extend(expressed_tf_targets_list)
                    all_expressed_ligands.extend(expressed_ligands)
                    cell_expressed_ligands.extend(expressed_ligands)

                    if "NKX" in key:
                        key = key.replace('-', '.')

                    for ligand in expressed_ligands:
                        tf_csv_data.append((key, p_value, tf_activity_dict[key][cell], ligand, cell_fix_name))

                    tf_data[key] = {  # "possible_targets": len(possible_targets),
                        # "possible_ligands": len(ligands),
                        "cell_type": cell_fix_name,
                        "p-value": p_value,
                        "z-score": tf_activity_dict[key][cell],
                        "num_expressed_targets": len(expressed_tf_targets_list),
                        "num_expressed_ligands": len(expressed_ligands),

                        "ligands": ', '.join(expressed_ligands)}

                    if key not in tfs_for_receptors:
                        tfs_for_receptors.append(key)

                    receptors_df = r2tf[r2tf['target_genesymbol'].isin([key])]
                    receptors_df = receptors_df.loc[:, ['source_genesymbol', 'target_genesymbol', 'n_references']]
                    receptors_list = receptors_df['source_genesymbol'].tolist()
                    for receptor in receptors_list:
                        tuple_table = (cell_fix_name, receptor, cell_fix_name, key, tf_activity_dict[key][cell])
                        r2tf_element_list.append(tuple_table)

            tf_dataframe = pd.DataFrame.from_dict(tf_data).T
            tf_dataframe_csv = tf_dataframe
            tf_dataframe = tf_dataframe.sort_values('num_expressed_ligands', ascending=False).head(10)
            table_file_name = out_dir + "/" + condition.replace(',', '_') + "_" + cell_fix_name + "_tf2ligands.csv"
            tf_dataframe_csv.to_csv(
                out_dir + "/" + condition.replace(',', '_') + "_" + cell_fix_name + "_tf2ligands.csv")
            cell_number_possible_targets = len(set(cell_possible_targets))
            cell_number_possible_ligands = len(set(cell_possible_ligands))
            cell_number_expressed_targets = len(set(cell_expressed_targets))
            cell_number_expressed_ligands = len(set(cell_expressed_ligands))
            cell_data.append({"cell_type": cell_fix_name, "active_tfs": cell_number_active_tfs,
                              "possible_targets": cell_number_possible_targets,
                              "possible_ligands": cell_number_possible_ligands,
                              "expressed_targets": cell_number_expressed_targets,
                              "expressed_ligands": cell_number_expressed_ligands,
                              "file_name": table_file_name,
                              "tf_data": tf_dataframe.to_html(classes='mystyle')})

        # receptors_df = r2tf[r2tf['target_genesymbol'].isin(tfs_for_receptors)]
        # receptors_df = receptors_df.loc[:, ['source_genesymbol', 'target_genesymbol', 'n_references']]
        # receptors_df.to_csv(condition_mapping[condition].replace(',', '_') + "_receptor2tf.csv")
        # csv_dataframe = pd.DataFrame(tf_csv_data, columns=['tf', 'p_value', 'z_score', 'ligand', 'cell_type'])
        r_df = pd.DataFrame(r2tf_element_list, columns=['sender', 'receptor', 'receiver', 'tf', 'z_score'])
        r_df.to_csv(out_dir + "/" + condition.replace(',', '_') + "_r2tf_network.csv")

        tfl_df = pd.DataFrame(tf_csv_data, columns=['tf', 'p_value', 'z_score', 'ligand', 'cell_type'])
        tfl_df.to_csv(out_dir + "/" + condition.replace(',', '_') + "_tf2ligand_network.csv")

        total_number_possible_targets = len(set(all_possible_targets))
        total_number_possible_ligands = len(set(all_possible_ligands))
        total_number_expressed_targets = len(set(all_expressed_targets))
        total_number_expressed_ligands = len(set(all_expressed_ligands))

        overview_data = [number_cell_types, number_active_tfs,
                         total_number_possible_targets, total_number_possible_ligands,
                         total_number_expressed_targets, total_number_expressed_ligands, cell_data[-1]["cell_type"],
                         file_dict["cluster_comparison_image"], file_dict["cluster_image"], file_dict["heatmap_image"]]

        overview_data_per_condition[condition] = overview_data
        cell_data_per_condition[condition] = cell_data

    ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
    env = Environment(loader=FileSystemLoader(ROOT_DIR))
    template = env.get_template("/html/dorothea_report.html")
    template_vars = {"title": "DoRothEA report (ABCDE)",
                     "conditions": conditions,
                     "overview_data": overview_data_per_condition,
                     "cell_data": cell_data_per_condition,
                     "name_mapping": condition_mapping
                     }
    html_out = template.render(template_vars)

    with open(out_dir + "/" + 'test.html', 'w') as f:
        f.write(html_out)

    HTML(string=html_out, base_url=ROOT_DIR).write_pdf(out_dir + "/" + "dorothea_r2tf2l_report.pdf",
                                                       stylesheets=[ROOT_DIR + "/html/style.css"])

    print("finish")
