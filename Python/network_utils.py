import pandas as pd
import networkx as nx


def network_from_tables(lr_table, rtf_table, tfl_table, out_dir, condition):

    lr_df = pd.read_csv(lr_table, sep=',', header=0)
    rtf_df = pd.read_csv(rtf_table, sep=',', header=0, index_col=0)
    tfl_df = pd.read_csv(tfl_table, sep=',', header=0, index_col=0)

    tab1 = pd.DataFrame()
    tab1['sender'] = tfl_df['cell_type']
    tab1['tf'] = tfl_df['tf']
    tab1['receiver'] = tfl_df['cell_type']
    tab1['ligand'] = tfl_df['ligand']
    tab1['z_score'] = tfl_df['z_score']

    graph_tfl_dir = nx.from_pandas_edgelist(tab1, source="tf", target="ligand",
                                            edge_attr=["sender", "receiver", "z_score"], create_using=nx.DiGraph())
    graph_rtf_dir = nx.from_pandas_edgelist(rtf_df, source="receptor", target="tf",
                                            edge_attr=["sender", "receiver", "z_score"], create_using=nx.DiGraph())
    graph_lr_dir = nx.from_pandas_edgelist(lr_df, source="Ligand", target="Receptor",
                                           edge_attr=["Ligand.Cluster", "Receptor.Cluster", "MeanLR"],
                                           create_using=nx.DiGraph())

    graph_lr_rtf_dir = nx.compose(graph_lr_dir, graph_rtf_dir)
    graph_lr_rtf_tfl_dir = nx.compose(graph_lr_rtf_dir, graph_tfl_dir)

    nx.write_graphml(graph_lr_rtf_tfl_dir,
                     out_dir + "/" + condition + "_l2r2tf2l_network.graphml")


def convert_table_to_crosstalker_input(lr_table, rtf_table, tfl_table, out_dir, condition):

    rtf_df = pd.read_csv(rtf_table, sep=',', header=0, index_col=0)
    tfl_df = pd.read_csv(tfl_table, sep=',', header=0, index_col=0)
    lr_df = pd.read_csv(lr_table, sep=',', header=0, index_col=0)

    rtf_list = []
    for i in range(rtf_df.shape[0]):
        curr = rtf_df.iloc[i, :]
        rtf_list.append((curr["sender"], curr["receiver"], curr["receptor"], "tf-" + curr["tf"], curr["z_score"]))
    result_rtf_df = pd.DataFrame(rtf_list,
                                 columns=["Ligand.Cluster", "Receptor.Cluster", "Ligand", "Receptor", "MeanLR"])
    result_rtf_df.to_csv(out_dir + "/" + condition + "_r2tf_crosstalker_input.csv")

    tfl_list = []
    for i in range(tfl_df.shape[0]):
        curr = tfl_df.iloc[i, :]
        tfl_list.append((curr["cell_type"], curr["cell_type"], "tf-" + curr["tf"], curr["ligand"], curr["z_score"]))
    result_tfl_df = pd.DataFrame(tfl_list,
                                 columns=["Ligand.Cluster", "Receptor.Cluster", "Ligand", "Receptor", "MeanLR"])
    result_tfl_df.to_csv(out_dir + "/" + condition + "_tf2l_crosstalker_input.csv")

    combined_df = lr_df.append(rtf_df)
    combined_df = combined_df.append(tfl_df)
    combined_df.to_csv(out_dir + "/" + condition + "_combined_crosstalker_input.csv")

    print("finish")
