import scanpy as sc
import liana as li
import pandas as pd
import os

def run_liana(ann_object_path, condition_field, cluster_field):
    data = sc.read_h5ad(ann_object_path)
    data.raw = data

    for i in set(data.obs[condition_field]):
        print(i)
        lr = li.method.cellphonedb(data[data.obs[condition_field] == i],
                                   groupby=cluster_field,
                                   expr_prop=0.0,
                                   verbose=True,
                                   resource_name='consensus', inplace=False)
        # lr.to_csv(f"{i}_lr_liana_consensus.csv")
        lr.to_csv(f"{i}_lr_liana_consensus_unfiltered.csv")


    for i in os.listdir():
        if i.endswith('lr_liana_consensus_unfiltered.csv'):
            evfull = pd.read_csv(i)
            evfull = evfull.loc[:, ['ligand', 'receptor', 'source', 'target', 'lr_means', 'cellphone_pvals']]
            evfull['type_gene_A'] = 'Ligand'
            evfull['type_gene_B'] = 'Receptor'
            evfull['gene_A'] = evfull['ligand']
            evfull['gene_B'] = evfull['receptor']
            evfull['MeanLR'] = evfull['lr_means']
            evfull = evfull.loc[list(evfull.cellphone_pvals.to_numpy() <= 0.05), :]
            evfull = evfull.loc[:, ['source', 'target', 'type_gene_A', 'type_gene_B', 'gene_A', 'gene_B', 'MeanLR']]
            k = i[0:i.find('_lr_')]
            evfull.to_csv(f'{k}_lr_ready.csv')
