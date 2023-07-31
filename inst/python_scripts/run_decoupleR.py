import scanpy as sc
import decoupler as dc
import pandas as pd

def run_decoupler_viper(ann_object_path, reg_path, out_path):
    ann_data = sc.read_h5ad(ann_object_path)
    reg = pd.read_csv(reg_path)

    dc.run_viper(mat=ann_data, net=reg, source='source', target='target', weight='weight', verbose=True, use_raw=False)

    estimates = ann_data.obsm['viper_estimate']
    estimates.to_csv(
        out_path + "decoupler_results.csv")
    print("Finished calculating TF activities")

def run_decoupler_ulm(ann_object_path, reg_path, out_path):
    ann_data = sc.read_h5ad(ann_object_path)
    reg = pd.read_csv(reg_path)

    dc.run_ulm(mat=ann_data, net=reg, source='source', target='target', weight='weight', verbose=True, use_raw=False)

    estimates = ann_data.obsm['ulm_estimate']
    estimates.to_csv(
        out_path + "decoupler_results.csv")
    print("Finished calculating TF activities")

def run_decoupler_mlm(ann_object_path, reg_path, out_path):
    ann_data = sc.read_h5ad(ann_object_path)
    reg = pd.read_csv(reg_path)

    dc.run_mlm(mat=ann_data, net=reg, source='source', target='target', weight='weight', verbose=True, use_raw=False)

    estimates = ann_data.obsm['mlm_estimate']
    estimates.to_csv(
        out_path + "decoupler_results.csv")
    print("Finished calculating TF activities")

def run_decoupler_wmean(ann_object_path, reg_path, out_path):
    ann_data = sc.read_h5ad(ann_object_path)
    reg = pd.read_csv(reg_path)

    dc.run_wmean(mat=ann_data, net=reg, source='source', target='target', weight='weight', verbose=True, use_raw=False)

    estimates = ann_data.obsm['wmean_estimate']
    estimates.to_csv(
        out_path + "decoupler_results.csv")
    print("Finished calculating TF activities")
